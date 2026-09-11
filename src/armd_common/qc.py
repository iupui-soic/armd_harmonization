"""QC and profiling helpers.

Every count that leaves the package for a report, log or shared file passes through cell(), so
counts from 1 to MIN_CELL - 1 never appear. Zero is not a small cell and is shown.
"""

from __future__ import annotations

from .config import MIN_CELL

SMALL = f"<{MIN_CELL}"


def cell(n):
    if n is None:
        return None
    n = int(n)
    return SMALL if 0 < n < MIN_CELL else n


def rate(numerator, denominator, digits: int = 4):
    """A proportion, or a marker when the numerator or denominator is a small cell."""
    if not denominator:
        return None
    if 0 < numerator < MIN_CELL or 0 < denominator < MIN_CELL:
        return SMALL
    return round(numerator / denominator, digits)
