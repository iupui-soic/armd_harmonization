"""CSV staging against synthetic files that reproduce release quirks."""

import duckdb
import pytest

from armd_common import convert
from armd_common.config import SiteConfig

HEADER = "anon_id,pat_enc_csn_id_coded,order_proc_id_coded,organism,code,days"


@pytest.fixture()
def site(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    (tmp_path / "raw").mkdir()
    return SiteConfig(site="test", raw_subdir="raw", parquet_subdir="pq", tables={"cohort": "c.csv"})


def write_csv(cfg, extra_lines, n=200):
    lines = [HEADER] + [f"p{i},e{i},{i},ESCHERICHIA COLI,{i},{i}" for i in range(n)] + extra_lines
    cfg.csv_path("cohort").write_text("\r\n".join(lines) + "\r\n")


def read(cfg, sql):
    return duckdb.connect().execute(sql.format(src=f"read_parquet('{cfg.parquet_glob('cohort')}')")).fetchall()


def test_quoted_comma_after_sniff_sample(site, monkeypatch):
    monkeypatch.setattr(convert, "SAMPLE_SIZE", 2048)
    write_csv(site, [
        'p900000,e900000,900000,"GRAM NEGATIVE RODS - NONFERM, UNIDENTIFIED",1,Null',
        'p900001,e900001,900001,"QUOTED NAME",2,7',
    ], n=60_000)
    try:
        leaked = duckdb.connect().execute(
            f"SELECT COUNT(*) FROM read_csv('{site.csv_path('cohort')}', header=true, sample_size=2048) WHERE organism LIKE '%\"%'"
        ).fetchone()[0]
        sniffing_alone_breaks = leaked > 0
    except duckdb.Error:
        sniffing_alone_breaks = True
    if not sniffing_alone_breaks:
        pytest.skip("this DuckDB version sniffed the late quotes, so the file does not reproduce the release quirk")

    result = convert.convert_table(duckdb.connect(), site, "cohort")
    assert result["rows"] == 60_002
    assert result["csv_options"]["quote"] == '"'
    assert dict(read(site, "SELECT anon_id, organism FROM {src} WHERE anon_id IN ('p900000', 'p900001')")) == {
        "p900000": "GRAM NEGATIVE RODS - NONFERM, UNIDENTIFIED",
        "p900001": "QUOTED NAME",
    }
    assert read(site, "SELECT COUNT(*) FROM {src} WHERE days IS NULL") == [(1,)]
    types = dict(result["columns"])
    assert types["anon_id"] == types["pat_enc_csn_id_coded"] == types["order_proc_id_coded"] == "VARCHAR"


def test_column_retried_as_varchar(site, monkeypatch):
    write_csv(site, ["p200,e200,200,ESCHERICHIA COLI,ABC,3"])
    original = convert._read_csv_sql

    def force_code_bigint(path, varchar_columns):
        sql = original(path, varchar_columns)
        return sql if "code" in varchar_columns else sql.replace("types = {", "types = {'code': 'BIGINT', ")

    monkeypatch.setattr(convert, "_read_csv_sql", force_code_bigint)
    result = convert.convert_table(duckdb.connect(), site, "cohort")
    assert result["retried_as_varchar"] == ["code"]
    assert dict(result["columns"])["code"] == "VARCHAR"
    assert read(site, "SELECT code FROM {src} WHERE anon_id = 'p200'") == [("ABC",)]


def test_conversion_error_message_is_sanitized(site):
    write_csv(site, ["p200,e200,200,ESCHERICHIA COLI,SECRETVALUE,3"])
    with pytest.raises(duckdb.Error) as err:
        duckdb.connect().execute(f"SELECT * FROM read_csv('{site.csv_path('cohort')}', header=true, types={{'code': 'BIGINT'}})").fetchall()
    assert convert._COLUMN_RE.search(str(err.value)).group(1) == "code"
    assert "SECRETVALUE" in str(err.value)
    assert "SECRETVALUE" not in convert._safe_message(err.value)


def test_structural_error_hides_row_content(site):
    write_csv(site, ["p200,e200,200,SECRET ORGANISM,1,2,unexpected"])
    with pytest.raises(convert.ConversionError) as err:
        convert.convert_table(duckdb.connect(), site, "cohort")
    assert "SECRET" not in str(err.value)
    assert not (site.parquet_dir / "cohort").exists()
