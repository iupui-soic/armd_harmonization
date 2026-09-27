"""Bounds for the population estimand (E3).

ARMD holds cultured events only. If a share u of empiric prescribing events has no culture, the
population contains N_u = N_c * u / (1 - u) uncultured events for every N_c cultured ones, with
unknown resistance and unknown model output.

manski_bounds gives worst-case (Manski) bounds on population prevalence, PPV and NPV from the
confusion counts of the cultured events, leaving the uncultured events unrestricted: each may be
flagged or not and resistant or not.

tipping_relative_risk adds assumptions to find where the population net benefit at threshold t
changes sign. Uncultured events are assumed to be flagged at the cultured rate, with the model's
sensitivity unchanged and their prevalence scaled by a relative risk RR. Population net benefit is
then linear in RR, and the function returns the RR at which it crosses zero, or None when no RR in
the feasible range does.
"""

from __future__ import annotations


def uncultured_events(n_cultured: float, u: float) -> float:
    if not 0 <= u < 1:
        raise ValueError("u must be in [0, 1)")
    return n_cultured * u / (1 - u)


def manski_bounds(tp: float, fp: float, tn: float, fn: float, u: float) -> dict:
    """Worst-case bounds on population prevalence, PPV and NPV. Each value is (lower, upper)."""
    n_c = tp + fp + tn + fn
    n_u = uncultured_events(n_c, u)
    n = n_c + n_u
    flagged, unflagged = tp + fp, tn + fn
    return {
        "uncultured_share": u,
        "prevalence": ((tp + fn) / n, (tp + fn + n_u) / n),
        # Every uncultured event flagged and susceptible, or flagged and resistant.
        "ppv": (tp / (flagged + n_u), (tp + n_u) / (flagged + n_u)),
        # Every uncultured event unflagged and resistant, or unflagged and susceptible.
        "npv": (tn / (unflagged + n_u), (tn + n_u) / (unflagged + n_u)),
    }


def tipping_relative_risk(tp: float, fp: float, tn: float, fn: float, u: float, t: float) -> dict:
    """The relative risk of resistance among uncultured events at which population NB(t) is zero."""
    n_c = tp + fp + tn + fn
    w = t / (1 - t)
    tp_rate, flag_rate, prevalence = tp / n_c, (tp + fp) / n_c, (tp + fn) / n_c
    nb_cultured = (tp - fp * w) / n_c
    out = {"uncultured_share": u, "threshold": t, "net_benefit_cultured": nb_cultured, "relative_risk": None}
    if u == 0 or tp_rate == 0:
        return out
    # NB_pop(RR) = (1 - u) NB_c + u [RR tp_rate (1 + w) - flag_rate w], which is zero at:
    rr = (u * flag_rate * w - (1 - u) * nb_cultured) / (u * tp_rate * (1 + w))
    # The uncultured true positives cannot exceed the flagged events, nor the prevalence exceed 1.
    ceiling = min(flag_rate / tp_rate, 1 / prevalence)
    if 0 <= rr <= ceiling:
        out["relative_risk"] = rr
    out["relative_risk_ceiling"] = ceiling
    return out
