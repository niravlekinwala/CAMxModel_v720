"""In-memory model state scaffolding."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np

from camx_py.state.grid import GridSpec


@dataclass
class GridState:
    grid: GridSpec
    tsurf: np.ndarray
    tsurf_next: np.ndarray
    tsurf_rate: np.ndarray
    snow: np.ndarray
    snowage: np.ndarray
    snowrat: np.ndarray
    height: np.ndarray
    height_next: np.ndarray
    height_rate: np.ndarray
    press: np.ndarray
    press_next: np.ndarray
    press_rate: np.ndarray
    windu: np.ndarray
    windu_next: np.ndarray
    windu_rate: np.ndarray
    windv: np.ndarray
    windv_next: np.ndarray
    windv_rate: np.ndarray
    tempk: np.ndarray
    temp_next: np.ndarray
    temp_rate: np.ndarray
    water: np.ndarray
    water_next: np.ndarray
    water_rate: np.ndarray
    rkv: np.ndarray
    rkv_next: np.ndarray
    rkv_rate: np.ndarray
    cldwtr: np.ndarray
    ranwtr: np.ndarray
    snowtr: np.ndarray
    gplwtr: np.ndarray
    cldod: np.ndarray
    cldph: np.ndarray
    cigfrc: np.ndarray
    cigtim: np.ndarray
    cigwtr: np.ndarray
    cigpcp: np.ndarray
    cigent: np.ndarray
    cigdet: np.ndarray
    conc: np.ndarray

    @classmethod
    def allocate(cls, grid: GridSpec, nspec: int) -> "GridState":
        ncol, nrow, nlay = grid.ncol, grid.nrow, grid.nlays
        zeros2d = lambda: np.zeros((ncol, nrow), dtype=np.float32)
        zeros3d = lambda: np.zeros((ncol, nrow, nlay), dtype=np.float32)
        return cls(
            grid=grid,
            tsurf=zeros2d(),
            tsurf_next=zeros2d(),
            tsurf_rate=zeros2d(),
            snow=zeros2d(),
            snowage=zeros2d(),
            snowrat=zeros2d(),
            height=zeros3d(),
            height_next=zeros3d(),
            height_rate=zeros3d(),
            press=zeros3d(),
            press_next=zeros3d(),
            press_rate=zeros3d(),
            windu=zeros3d(),
            windu_next=zeros3d(),
            windu_rate=zeros3d(),
            windv=zeros3d(),
            windv_next=zeros3d(),
            windv_rate=zeros3d(),
            tempk=zeros3d(),
            temp_next=zeros3d(),
            temp_rate=zeros3d(),
            water=zeros3d(),
            water_next=zeros3d(),
            water_rate=zeros3d(),
            rkv=zeros3d(),
            rkv_next=zeros3d(),
            rkv_rate=zeros3d(),
            cldwtr=zeros3d(),
            ranwtr=zeros3d(),
            snowtr=zeros3d(),
            gplwtr=zeros3d(),
            cldod=zeros3d(),
            cldph=zeros3d(),
            cigfrc=zeros2d(),
            cigtim=zeros2d(),
            cigwtr=zeros3d(),
            cigpcp=zeros3d(),
            cigent=zeros3d(),
            cigdet=zeros3d(),
            conc=np.zeros((ncol, nrow, nlay, nspec), dtype=np.float32),
        )


@dataclass
class ModelState:
    grids: List[GridState]

    @classmethod
    def allocate(cls, grids: List[GridSpec], nspec: int) -> "ModelState":
        return cls([GridState.allocate(grid, nspec) for grid in grids])

