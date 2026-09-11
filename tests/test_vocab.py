import duckdb
import pytest

from armd_common import organisms, vocab


@pytest.mark.parametrize(
    "raw, group, idsa, ambiguous",
    [
        ("ESCHERICHIA COLI", "ecoli", True, False),
        ("ESCHERICHIA COLI (CARBAPENEM RESISTANT)", "ecoli", True, False),
        ("ESCHERICHIA HERMANII", "other_enterobacterales", False, False),
        ("KLEBSIELLA PNEUMONIAE", "klebsiella", True, False),
        ("KLEBSIELLA PNEUMONIAE SSP. OZAENAE", "klebsiella", True, False),
        ("KLEBSIELLA OXYTOCA", "klebsiella", True, False),
        ("ZZZENTEROBACTER AEROGENES", "klebsiella", False, False),
        ("KLEBSIELLA AEROGENES", "klebsiella", False, False),
        ("PROTEUS MIRABILIS", "proteus", True, False),
        ("ZZZPROTEUS VULGARIS", "proteus", False, False),
        ("ENTEROBACTER CLOACAE COMPLEX (CARBAPENEM RESISTANT)", "enterobacter", False, False),
        ("CITROBACTER KOSERI (DIVERSUS)", "other_enterobacterales", False, False),
        ("SERRATIA MARCESCENS", "other_enterobacterales", False, False),
        ("MORGANELLA MORGANII", "other_enterobacterales", False, False),
        ("SALMONELLA TYPHI/PARATYPHI", "other_enterobacterales", False, False),
        ("LACTOSE FERMENTING GRAM NEGATIVE RODS", "other_enterobacterales", False, True),
        ("ENTERIC GRAM NEGATIVE RODS", "other_enterobacterales", False, True),
        ("NON LACTOSE FERMENTING GRAM NEGATIVE RODS", "non_enterobacterales", False, True),
        ("NON FERMENTING GRAM NEGATIVE RODS", "non_enterobacterales", False, False),
        ("GRAM NEGATIVE RODS - NONFERM, UNIDENTIFIED", "non_enterobacterales", False, False),
        ("GRAM NEGATIVE ORGANISM", "non_enterobacterales", False, True),
        ("MIXED GENITAL FLORA", "non_enterobacterales", False, False),
        ("CULTURE POSITIVE GRAM NEGATIVE COCCI", "non_enterobacterales", False, False),
        ("PSEUDOMONAS AERUGINOSA (MUCOID)", "non_enterobacterales", False, False),
        ("STAPH AUREUS {MRSA}", "non_enterobacterales", False, False),
        ("ENTEROCOCCUS FAECALIS", "non_enterobacterales", False, False),
        ("CANDIDA ALBICANS", "non_enterobacterales", False, False),
    ],
)
def test_classify(raw, group, idsa, ambiguous):
    info = organisms.classify(raw)
    assert (info.group, info.idsa_esbl_species, info.ambiguous) == (group, idsa, ambiguous)
    assert info.is_enterobacterales == (group in organisms.ENTEROBACTERALES_GROUPS)


def test_classify_missing():
    assert organisms.classify(None).group is None


def test_normalize_name():
    assert organisms.normalize_name("STAPH AUREUS(COLONY VARIANT - SMALL COLONY OR OTHER MORPHOTYPE)") == "STAPHYLOCOCCUS AUREUS"
    assert organisms.normalize_name("STREPTOCOCCUS DYSGALACTIAE SUBSP. EQUISIMILIS (GROUP C/G STR]") == "STREPTOCOCCUS DYSGALACTIAE SUBSP. EQUISIMILIS"


def test_register_lookup():
    con = duckdb.connect()
    table = organisms.register_lookup(con, ["ESCHERICHIA COLI", None, "CANDIDA ALBICANS", "ESCHERICHIA COLI"])
    assert con.execute(f"SELECT COUNT(*), COUNT(*) FILTER (WHERE is_enterobacterales) FROM {table}").fetchone() == (2, 1)


@pytest.mark.parametrize(
    "name, code",
    [
        ("Ceftriaxone", "CRO"),
        ("Piperacillin/Tazobactam", "TZP"),
        ("Trimethoprim/Sulfamethoxazole", "SXT"),
        ("Amoxicillin/Clavulanic Acid", "AMC"),
        ("Ceftazidime/Avibactam", "CAZ_AVI"),
        ("Ticacarcillin/Clavulanic Acid", "TIM"),
        ("Imipenem/Ebactam", "IPM_REL"),
        ("Cephalexin/Cephalothin", "CEF"),
        ("ceftolozane_tazobactam", "TOL_TAZ"),
        ("Ampicillin/Sulbactam", "SAM"),
        ("Pip/tazobactam", "TZP"),
        ("Amox/clavulanate", "AMC"),
        ("Augmentin", "AMC"),
        ("Amp/sulbactam", "SAM"),
        ("Trimeth/sulf", "SXT"),
        ("Trimeth-sulfamethoxazole", "SXT"),
        ("Ceftazidime+clavulanic Acid", "CAZ_CLV"),
        ("Cefotaxime+clavulanic Acid", "CTX_CLV"),
        ("Mcim", "MCIM"),
        ("Polymixin e", "CST"),
        ("5 fluorocytosine", "5FC"),
        ("Synercid", "QUD"),
        ("Gentamicin screen", None),
        ("Not A Drug", None),
    ],
)
def test_antibiotic_code(name, code):
    assert vocab.antibiotic_code(name) == code


def test_unmapped_by_design_names_stay_unmapped():
    assert all(vocab.antibiotic_code(name.replace("_", " ")) is None for name in vocab.UNMAPPED_BY_DESIGN)


def test_sentinel_agents_and_aliases_resolve_to_known_codes():
    assert set(vocab.SENTINEL_AGENTS) <= set(vocab.CODES)
    assert set(vocab.ALIASES.values()) <= set(vocab.CODES)


@pytest.mark.parametrize(
    "raw, value",
    [
        ("Susceptible", "S"),
        ("Intermediate", "I"),
        ("Resistant", "R"),
        ("Susceptible dose dependent", "I"),
        ("Susceptible dose-dependent", "I"),
        ("Non-susceptible", "R"),
        ("Inconclusive", "other"),
        ("Synergism", "other"),
        ("something new", "other"),
        (None, None),
    ],
)
def test_sir(raw, value):
    assert vocab.sir(raw) == value


@pytest.mark.parametrize("raw, band", [("18-24 years", "18-24"), ("85-89 years", "85-89"), ("above 90", "90+"), ("90+", "90+"), ("10-17 years", None), (None, None)])
def test_age_band(raw, band):
    assert vocab.age_band(raw) == band


def test_sex():
    assert [vocab.sex(v) for v in ("Female", "Male", "Unknown", "X", None)] == ["f", "m", "unknown", "unknown", "unknown"]


def test_sql_case_matches_python():
    con = duckdb.connect()
    expr = vocab.sql_case("v", vocab.SIR, default="'other'")
    for raw in ("Susceptible", " RESISTANT ", "Susceptible dose dependent", "Inconclusive", "weird"):
        assert con.execute(f"SELECT {expr} FROM (SELECT ? AS v)", [raw]).fetchone()[0] == vocab.sir(raw)
