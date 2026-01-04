"""Advection stubs (CPU/GPU backends to be implemented)."""
from __future__ import annotations

import numpy as np


def advect_ppm(conc: np.ndarray, wind_u: np.ndarray, wind_v: np.ndarray, dt: float, dx: float, dy: float) -> np.ndarray:
    """Placeholder for PPM advection; returns updated concentration array."""
    raise NotImplementedError("Advection kernel not yet ported.")

