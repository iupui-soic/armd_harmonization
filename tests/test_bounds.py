"""E3 bounds by hand."""

import pytest

from armd_common import bounds


def test_manski_bounds_by_hand():
    # 1,000 cultured events: 40 TP, 60 FP, 880 TN, 20 FN. u = 0.5 adds 1,000 uncultured events.
    b = bounds.manski_bounds(40, 60, 880, 20, 0.5)
    assert b["prevalence"] == pytest.approx((60 / 2000, 1060 / 2000))
    assert b["ppv"] == pytest.approx((40 / 1100, 1040 / 1100))
    assert b["npv"] == pytest.approx((880 / 1900, 1880 / 1900))


def test_no_uncultured_events_gives_point_estimates():
    b = bounds.manski_bounds(40, 60, 880, 20, 0.0)
    assert b["ppv"] == pytest.approx((0.4, 0.4)) and b["prevalence"] == pytest.approx((0.06, 0.06))
    with pytest.raises(ValueError):
        bounds.manski_bounds(1, 1, 1, 1, 1.0)


def test_tipping_relative_risk_zeroes_population_net_benefit():
    # A weak model: NB at t = 0.20 is barely positive among cultured events.
    tp, fp, tn, fn, u, t = 20, 60, 900, 20, 0.5, 0.20
    out = bounds.tipping_relative_risk(tp, fp, tn, fn, u, t)
    rr, w, n = out["relative_risk"], t / (1 - t), 1000
    assert rr == pytest.approx(0.6)
    assert out["net_benefit_cultured"] == pytest.approx((tp - fp * w) / n)
    nb_uncultured = rr * tp / n * (1 + w) - (tp + fp) / n * w
    assert (1 - u) * out["net_benefit_cultured"] + u * nb_uncultured == pytest.approx(0, abs=1e-12)
    # Net benefit is positive among cultured events, so it only turns negative if the uncultured
    # carry less risk.
    assert 0 < rr < 1


def test_tipping_point_outside_the_feasible_range():
    # A useful model with few uncultured events: even no resistance among them keeps NB positive.
    out = bounds.tipping_relative_risk(40, 60, 880, 20, 0.25, 0.10)
    assert out["net_benefit_cultured"] > 0 and out["relative_risk"] is None
