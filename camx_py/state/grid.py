"""Grid and projection metadata for CAMx Python port."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Tuple


Projection = Literal["latlon", "utm", "lambert", "rotated_polar", "polar", "mercator"]


@dataclass
class GridSpec:
    """Spatial grid definition matching CAMx ncf_set_global semantics."""

    ncol: int
    nrow: int
    nlays: int
    xorig_km: float
    yorig_km: float
    delx_km: float
    dely_km: float
    projection: Projection
    polelon: float
    polelat: float
    tlat1: float
    tlat2: float
    mesh_factor: int = 1
    nested_bounds: Tuple[int, int, int, int] | None = None  # inst1, inst2, jnst1, jnst2
    staggered_winds: bool = False
    latlon_coords: Tuple = ()  # Optional precomputed (lon, lat) arrays

    def inner_dims(self) -> Tuple[int, int]:
        """Return COL/ROW dimension sizes following CAMx nested grid trimming."""
        if self.mesh_factor == 1 and not self.nested_bounds:
            return self.ncol, self.nrow
        # nested grids drop the halo cells
        return self.ncol - 2, self.nrow - 2
