"""Time-rate of change helpers."""
from __future__ import annotations

import numpy as np


def timrates(current: np.ndarray, nxt: np.ndarray, dt_minutes: float) -> np.ndarray:
    """Compute time rates (nxt - current) / dt_seconds."""
    dt_seconds = dt_minutes * 60.0
    if dt_seconds <= 0:
        raise ValueError("dt_minutes must be positive")
    return (nxt - current) / dt_seconds

