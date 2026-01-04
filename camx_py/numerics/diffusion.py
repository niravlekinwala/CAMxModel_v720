"""Diffusion stubs (CPU/GPU backends to be implemented)."""
from __future__ import annotations

import numpy as np


def vertical_diffusion(conc: np.ndarray, k_profile: np.ndarray, dt: float, dz: np.ndarray) -> np.ndarray:
    """Placeholder for vertical diffusion; returns updated concentration array."""
    raise NotImplementedError("Diffusion kernel not yet ported.")

