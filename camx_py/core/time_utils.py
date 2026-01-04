"""Time utilities for CAMx-style YYJJJ + HHMM timestamps."""
from __future__ import annotations

from typing import Tuple


def add_minutes(date_yyjjj: int, time_hhmm: float, minutes: float) -> Tuple[int, float]:
    """Advance YYJJJ/HHMM by a number of minutes."""
    hour = int(time_hhmm // 100)
    minute = int(round(time_hhmm % 100))
    total = hour * 60 + minute + int(round(minutes))
    new_hour = (total // 60) % 24
    new_minute = total % 60
    day_incr = total // (24 * 60)
    new_date = date_yyjjj
    for _ in range(day_incr):
        new_date = _add_day(new_date)
    return new_date, float(new_hour * 100 + new_minute)


def is_before_or_equal(
    date_a: int, time_a: float, date_b: int, time_b: float
) -> bool:
    """Return True if A <= B using YYJJJ + HHMM ordering."""
    if date_a < date_b:
        return True
    if date_a > date_b:
        return False
    return time_a <= time_b


def is_after_or_equal(
    date_a: int, time_a: float, date_b: int, time_b: float
) -> bool:
    """Return True if A >= B using YYJJJ + HHMM ordering."""
    if date_a > date_b:
        return True
    if date_a < date_b:
        return False
    return time_a >= time_b


def _add_day(date_yyjjj: int) -> int:
    """Increment YYJJJ by one day with simple leap-year handling."""
    year = date_yyjjj // 1000
    jday = date_yyjjj % 1000
    leap = (year % 4 == 0)
    jday += 1
    if jday > 365:
        if leap:
            if jday == 367:
                jday = 1
                year += 1
        else:
            jday = 1
            year += 1
    return year * 1000 + jday
