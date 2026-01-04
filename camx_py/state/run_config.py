"""Run-level configuration derived from CAMx namelist."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from .grid import GridSpec
from .species import SpeciesDef


@dataclass
class RunConfig:
    """Minimal runtime configuration to drive NetCDF I/O and loop setup."""

    begin_date: int  # YYJJJ
    begin_time: float  # HHMM
    end_date: int
    end_time: float
    start_year: int
    end_year: int
    dtout_minutes: float
    dtinp_minutes: float
    grids: List[GridSpec]
    species: List[SpeciesDef]
    run_message: str = ""
    version_string: str = "CAMx v7.20"
    time_zone: int = 0
    probing_tool: str = ""
    chem_solver: str = ""
    advection_solver: str = ""
    drydep_model: str = ""
    wetdep: bool = False
    acm2: bool = False
    cig_model: bool = False
    surface_model: bool = False
    inline_ix_emiss: bool = False
    bidi_nh3_drydep: bool = False
    strat_o3_profile: bool = False
    super_stepping: bool = False
    gridded_emiss: bool = False
    point_emiss: bool = False
    ignore_emiss_dates: bool = False
    output_3d: bool = True
    pig_sample_grid: bool = False
    pig_sample_bckgnd: bool = False
    met2d_paths: List[str] = field(default_factory=list)
    met3d_paths: List[str] = field(default_factory=list)
    cloud_paths: List[str] = field(default_factory=list)
    kv_paths: List[str] = field(default_factory=list)
    ic_paths: List[str] = field(default_factory=list)
    bnd_paths: List[str] = field(default_factory=list)
    top_paths: List[str] = field(default_factory=list)
