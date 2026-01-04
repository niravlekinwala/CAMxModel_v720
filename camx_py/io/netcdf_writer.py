"""NetCDF writer aligned with CAMx v7.20 Fortran implementation.

This mirrors the Fortran routines:
  - ncf_createfile.F90
  - ncf_set_global.f / ncf_wrt_global.f
  - ncf_wrt_dim.f
  - ncf_wrt_vars_base.f / ncf_set_vars_base.f
  - ncf_wrt_vars_species.F90
  - ncf_wrt_data_grid.f
  - ncf_wrt_data_tstep.f
  - ncf_wrt_data_species.f
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
from netCDF4 import Dataset

from camx_py.state.grid import GridSpec
from camx_py.state.run_config import RunConfig
from camx_py.state.species import SpeciesDef


@dataclass
class Chunking:
    """Chunking and compression settings."""

    compress: bool = False
    deflate_level: int = 0
    chunk_x: int = 64
    chunk_y: int = 64


def _now_ioapi() -> Tuple[int, int]:
    """Return current date/time in IOAPI YYJJJ, HHMMSS format."""
    now = _dt.datetime.utcnow()
    y = now.year % 100
    jday = now.timetuple().tm_yday
    date = y * 1000 + jday
    time = now.hour * 10000 + now.minute * 100 + now.second
    return date, time


class CamxNetCDFWriter:
    """Produces CAMx-compatible NetCDF files."""

    def __init__(self, cfg: RunConfig, chunking: Chunking | None = None):
        self.cfg = cfg
        self.chunking = chunking or Chunking()

    # ---------------------- Public API ---------------------- #
    def create_file(
        self,
        path: Path | str,
        grid: GridSpec,
        is_sample: bool = False,
        var_dim_is_layers: bool = True,
    ) -> Dataset:
        """Create an empty file with dimensions, vars, and attributes defined."""
        nc = Dataset(path, "w", format="NETCDF4_CLASSIC" if self.chunking.compress else "NETCDF3_CLASSIC")
        dims = self._write_dimensions(nc, grid, len(self.cfg.species))
        self._set_global_attrs(nc, grid, is_sample, len(self.cfg.species))
        self._define_base_vars(nc, grid, dims)
        self._define_species_vars(nc, grid, dims, var_dim_is_layers)
        return nc

    def write_grid_data(self, nc: Dataset, grid: GridSpec, lon: np.ndarray, lat: np.ndarray, topo: np.ndarray) -> None:
        """Write static grid variables X, Y, layer, longitude, latitude, topo."""
        col = nc.dimensions["COL"].size
        row = nc.dimensions["ROW"].size
        # X/Y in km unless latlon (degrees)
        factor = 1.0 if grid.projection == "latlon" else 1000.0
        x = (grid.xorig_km * factor) + (np.arange(grid.ncol) + 0.5) * grid.delx_km * factor
        y = (grid.yorig_km * factor) + (np.arange(grid.nrow) + 0.5) * grid.dely_km * factor
        if col != grid.ncol:
            x = x[1:-1]
        if row != grid.nrow:
            y = y[1:-1]
        nc.variables["X"][:] = x
        nc.variables["Y"][:] = y
        nc.variables["layer"][:] = np.arange(1, grid.nlays + 1, dtype=np.int32)
        # trim halo for nested grids
        lon_trim = lon
        lat_trim = lat
        topo_trim = topo
        if col != grid.ncol or row != grid.nrow:
            lon_trim = lon[1:-1, 1:-1]
            lat_trim = lat[1:-1, 1:-1]
            topo_trim = topo[1:-1, 1:-1]
        nc.variables["longitude"][:, :] = lon_trim.astype(np.float64)
        nc.variables["latitude"][:, :] = lat_trim.astype(np.float64)
        nc.variables["topo"][:, :] = topo_trim.astype(np.float32)

    def write_tflags(self, nc: Dataset, tstep_index: int, start_flag: Tuple[int, int], end_flag: Tuple[int, int]) -> None:
        """Write TFLAG/ETFLAG for a timestep."""
        nvars = nc.dimensions["VAR"].size
        nc.variables["TFLAG"][0:2, 0:nvars, tstep_index] = np.array([[start_flag[0]] * nvars, [start_flag[1]] * nvars], dtype=np.int32)
        nc.variables["ETFLAG"][0:2, 0:nvars, tstep_index] = np.array([[end_flag[0]] * nvars, [end_flag[1]] * nvars], dtype=np.int32)

    def write_species(
        self,
        nc: Dataset,
        grid: GridSpec,
        tstep_index: int,
        height: np.ndarray,
        species_fields: Sequence[np.ndarray],
        is_output: Sequence[bool],
        num_dims: int,
    ) -> None:
        """Write z and species fields for one timestep.

        height: shape (ncol, nrow, nlays_out)
        species_fields: list of arrays matching num_dims (3 or 4)
        """
        col = nc.dimensions["COL"].size
        row = nc.dimensions["ROW"].size
        # z first
        z_var = nc.variables["z"]
        z_data = height
        if col != grid.ncol or row != grid.nrow:
            z_data = height[1:-1, 1:-1, :]
        z_var[:, :, :, tstep_index] = z_data.astype(np.float32)

        for spec_def, field, write_flag in zip(self.cfg.species, species_fields, is_output):
            if not write_flag:
                continue
            name = spec_def.name.strip()
            var = nc.variables[name]
            data = field
            if col != grid.ncol or row != grid.nrow:
                data = field[1:-1, 1:-1, ...]
            if num_dims == 4:
                var[:, :, :, tstep_index] = data.astype(np.float32)
            else:
                var[:, :, tstep_index] = data.astype(np.float32)

    # ---------------------- Private helpers ---------------------- #
    def _write_dimensions(self, nc: Dataset, grid: GridSpec, nspcs: int) -> Tuple[int, int, int, int, int, int]:
        ncf_col = grid.ncol if grid.mesh_factor == 1 and not grid.nested_bounds else grid.ncol - 2
        ncf_row = grid.nrow if grid.mesh_factor == 1 and not grid.nested_bounds else grid.nrow - 2
        nc.createDimension("TSTEP", None)
        nc.createDimension("DATE-TIME", 2)
        nc.createDimension("LAY", grid.nlays)
        nc.createDimension("VAR", nspcs + 1)  # z + species
        nc.createDimension("ROW", ncf_row)
        nc.createDimension("COL", ncf_col)
        return ncf_col, ncf_row, grid.nlays, nspcs + 1, 2, None

    def _set_global_attrs(self, nc: Dataset, grid: GridSpec, is_sample: bool, nspcs: int) -> None:
        cdate, ctime = _now_ioapi()
        ncf_nrows = grid.nrow if grid.mesh_factor == 1 and not grid.nested_bounds else grid.nrow - 2
        ncf_ncols = grid.ncol if grid.mesh_factor == 1 and not grid.nested_bounds else grid.ncol - 2
        nc.setncattr("FTYPE", 1)
        nc.setncattr("CDATE", cdate)
        nc.setncattr("CTIME", ctime)
        nc.setncattr("WDATE", cdate)
        nc.setncattr("WTIME", ctime)
        nc.setncattr("SDATE", self.cfg.begin_date)
        nc.setncattr("STIME", int(self.cfg.begin_time * 100))
        if self.cfg.dtout_minutes < 60.0:
            tstep = int(self.cfg.dtout_minutes * 100)
        else:
            tstep = int(self.cfg.dtout_minutes / 60.0 * 10000)
        nc.setncattr("TSTEP", tstep)
        nc.setncattr("NTHIK", 1)
        nc.setncattr("NCOLS", ncf_ncols)
        nc.setncattr("NROWS", ncf_nrows)
        nc.setncattr("NLAYS", grid.nlays)
        nc.setncattr("NVARS", nspcs + 1)
        # projection attributes
        proj_map = {
            "latlon": (0, 1),
            "utm": (1, 5),
            "lambert": (2, 2),
            "rotated_polar": (3, 4),
            "polar": (4, 6),
            "mercator": (5, 7),
        }
        cproj, gdtyp = proj_map.get(grid.projection, (0, 1))
        nc.setncattr("GDTYP", gdtyp)
        nc.setncattr("P_ALP", float(min(grid.tlat1, grid.tlat2) if grid.tlat1 > 0 else max(grid.tlat1, grid.tlat2)))
        nc.setncattr("P_BET", float(max(grid.tlat1, grid.tlat2) if grid.tlat1 > 0 else min(grid.tlat1, grid.tlat2)))
        nc.setncattr("P_GAM", float(grid.polelon))
        nc.setncattr("XCENT", float(grid.polelon))
        nc.setncattr("YCENT", float(grid.polelat))
        nc.setncattr("XORIG", float(grid.xorig_km * 1000.0 if grid.projection != "latlon" else grid.xorig_km))
        nc.setncattr("YORIG", float(grid.yorig_km * 1000.0 if grid.projection != "latlon" else grid.yorig_km))
        nc.setncattr("XCELL", float(grid.delx_km * 1000.0 if grid.projection != "latlon" else grid.delx_km))
        nc.setncattr("YCELL", float(grid.dely_km * 1000.0 if grid.projection != "latlon" else grid.dely_km))
        nc.setncattr("VGTYP", 6)
        nc.setncattr("VGTOP", 10000.0)
        nc.setncattr("VGLVLS", np.zeros(grid.nlays + 1, dtype=np.int32))
        nc.setncattr("GDNAM", self.cfg.version_string)
        nc.setncattr("UPNAM", self.cfg.version_string)
        nc.setncattr("IOAPI_VERSION", "IOAPI-CAMx")
        nc.setncattr("EXEC_ID", "CAMx")
        # VAR-LIST
        var_list = ["z"] + [s.name.strip() for s in self.cfg.species]
        var_concat = "".join([v.ljust(16) for v in var_list])
        nc.setncattr("VAR-LIST", var_concat)
        nc.setncattr("FILEDESC", self.cfg.run_message or self.cfg.version_string)
        nc.setncattr("HISTORY", f"Generated by {self.cfg.version_string}")
        nc.setncattr("IUTM", 0)
        nc.setncattr("ISTAG", 1 if grid.staggered_winds else 0)
        nc.setncattr("CPROJ", cproj)
        nc.setncattr("NSTEPS", 0)
        nc.setncattr("CAMx_NAME", self.cfg.version_string)
        nc.setncattr("NOTE", self.cfg.run_message)
        nc.setncattr("ITZON", self.cfg.time_zone)
        nc.setncattr("UPDSC", self.cfg.version_string)
        nc.setncattr("GRID_ID", 1)
        inst1, inst2, jnst1, jnst2 = 1, grid.ncol, 1, grid.nrow
        if grid.nested_bounds:
            inst1, inst2, jnst1, jnst2 = grid.nested_bounds
        nc.setncattr("I_GRID_START", inst1)
        nc.setncattr("I_GRID_END", inst2)
        nc.setncattr("J_GRID_START", jnst1)
        nc.setncattr("J_GRID_END", jnst2)
        nc.setncattr("GRID_MESH_FACTOR", grid.mesh_factor)
        nc.setncattr("FLEXI_NEST", 0)
        # model flags
        nc.setncattr("ADVECTION", self.cfg.advection_solver)
        nc.setncattr("CHEM_SOLVER", self.cfg.chem_solver)
        nc.setncattr("PIG", "")
        nc.setncattr("PROBING_TOOL", self.cfg.probing_tool)
        nc.setncattr("CHEMISTRY", "")
        nc.setncattr("TOTAL_SPECIES", len(self.cfg.species))
        nc.setncattr("RADICAL_SPECIES", 0)
        nc.setncattr("GAS_SPECIES", 0)
        nc.setncattr("PM_SPECIES", 0)
        nc.setncattr("REACTIONS", 0)
        nc.setncattr("DRYDEP", self.cfg.drydep_model)
        nc.setncattr("WETDEP", 1 if self.cfg.wetdep else 0)
        nc.setncattr("ACM2", 1 if self.cfg.acm2 else 0)
        nc.setncattr("CIG_MODEL", 1 if self.cfg.cig_model else 0)
        nc.setncattr("SURFACE_MODEL", 1 if self.cfg.surface_model else 0)
        nc.setncattr("INLINE_IX_EMISS", 1 if self.cfg.inline_ix_emiss else 0)
        nc.setncattr("BIDI_NH3_DRYDEP", 1 if self.cfg.bidi_nh3_drydep else 0)
        nc.setncattr("STRAT_O3_PROFILE", 1 if self.cfg.strat_o3_profile else 0)
        nc.setncattr("SUPER_STEPPING", 1 if self.cfg.super_stepping else 0)
        nc.setncattr("GRIDDED_EMISS", 1 if self.cfg.gridded_emiss else 0)
        nc.setncattr("POINT_EMISS", 1 if self.cfg.point_emiss else 0)
        nc.setncattr("IGNORE_EMISS_DATES", 1 if self.cfg.ignore_emiss_dates else 0)
        nc.setncattr("OUTPUT_3D", 1 if self.cfg.output_3d else 0)
        nc.setncattr("PIG_SAMPLE_GRID", 1 if self.cfg.pig_sample_grid else 0)
        if is_sample:
            nc.setncattr("PIG_SAMPLE_GRID_ID", 1)
            nc.setncattr("I_SAMPLE_START", inst1)
            nc.setncattr("I_SAMPLE_END", inst2)
            nc.setncattr("J_SAMPLE_START", jnst1)
            nc.setncattr("J_SAMPLE_END", jnst2)
        nc.setncattr("PIG_SAMPLE_BCKGND", 1 if self.cfg.pig_sample_bckgnd else 0)

    def _define_base_vars(self, nc: Dataset, grid: GridSpec, dims: Tuple[int, int, int, int, int, int]) -> None:
        col_dim = nc.dimensions["COL"]
        row_dim = nc.dimensions["ROW"]
        lay_dim = nc.dimensions["LAY"]
        tstep_dim = nc.dimensions["TSTEP"]
        date_time_dim = nc.dimensions["DATE-TIME"]
        var_dim = nc.dimensions["VAR"]

        nc.createVariable("X", "f8", ("COL",)).setncatts(
            {"units": "degrees" if grid.projection == "latlon" else "km",
             "long_name": "X coordinate",
             "var_desc": "Longitude degrees east" if grid.projection == "latlon" else "X cartesian distance from projection origin"}
        )
        nc.createVariable("Y", "f8", ("ROW",)).setncatts(
            {"units": "degrees" if grid.projection == "latlon" else "km",
             "long_name": "Y coordinate",
             "var_desc": "Longitude degrees north" if grid.projection == "latlon" else "Y cartesian distance from projection origin"}
        )
        nc.createVariable("layer", "i4", ("LAY",)).setncatts(
            {"units": "Layer index", "long_name": "Model layer", "var_desc": "Model layer"}
        )
        nc.createVariable("TFLAG", "i4", ("DATE-TIME", "VAR", "TSTEP")).setncatts(
            {"units": "YYYYDDD,HHMMSS", "long_name": "Start time flag", "var_desc": "Timestep start date and time"}
        )
        nc.createVariable("ETFLAG", "i4", ("DATE-TIME", "VAR", "TSTEP")).setncatts(
            {"units": "YYYYDDD,HHMMSS", "long_name": "End time flag", "var_desc": "Timestep end date and time"}
        )
        nc.createVariable("longitude", "f8", ("COL", "ROW")).setncatts(
            {"units": "Degrees east", "long_name": "Longitude", "var_desc": "Longitude degrees east", "coordinates": "latitude longitude"}
        )
        nc.createVariable("latitude", "f8", ("COL", "ROW")).setncatts(
            {"units": "Degrees north", "long_name": "Latitude", "var_desc": "Latitude degrees north", "coordinates": "latitude longitude"}
        )
        nc.createVariable("topo", "f4", ("COL", "ROW")).setncatts(
            {"units": "m MSL", "long_name": "topographic elevation", "var_desc": "topographic elevation m above sea level", "coordinates": "latitude longitude"}
        )
        nc.createVariable("z", "f4", ("COL", "ROW", "LAY", "TSTEP")).setncatts(
            {"units": "m", "long_name": "Layer height", "var_desc": "Layer interface heights AGL", "coordinates": "latitude longitude"}
        )

    def _define_species_vars(self, nc: Dataset, grid: GridSpec, dims: Tuple[int, int, int, int, int, int], var_dim_is_layers: bool) -> None:
        col_dim = nc.dimensions["COL"]
        row_dim = nc.dimensions["ROW"]
        lay_dim = nc.dimensions["LAY"]
        tstep_dim = nc.dimensions["TSTEP"]
        for spec in self.cfg.species:
            dims_tuple = ("COL", "ROW", "LAY", "TSTEP") if var_dim_is_layers else ("COL", "ROW", "TSTEP")
            var = nc.createVariable(spec.name.strip(), "f4", dims_tuple, zlib=self.chunking.compress, complevel=self.chunking.deflate_level if self.chunking.compress else 0)
            if self.chunking.compress:
                chunk_x = min(self.chunking.chunk_x, col_dim.size)
                chunk_y = min(self.chunking.chunk_y, row_dim.size)
                chunk_z = max(1, lay_dim.size // 2) if "LAY" in dims_tuple else 1
                var.set_chunksizes((chunk_x, chunk_y, chunk_z, 1) if dims_tuple == ("COL", "ROW", "LAY", "TSTEP") else (chunk_x, chunk_y, 1))
            var.setncatts(
                {
                    "long_name": spec.long_name,
                    "units": spec.units,
                    "var_desc": spec.description,
                    "coordinates": spec.coordinates,
                }
            )

