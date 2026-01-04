"""Top-level CAMx Python driver scaffolding."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from camx_py.config import NamelistParser, map_namelist_to_config
from camx_py.io.netcdf_reader import NetCDFReader
from camx_py.io.netcdf_writer import CamxNetCDFWriter
from camx_py.core.time_utils import add_minutes, is_after_or_equal, is_before_or_equal
from camx_py.numerics.timrates import timrates
from camx_py.state.model_state import ModelState
from camx_py.state.run_config import RunConfig

@dataclass
class DriverSettings:
    camx_in_path: str = "CAMx.in"
    output_root: str = "CAMxPy"
    use_gpu: bool = False
    optimized: bool = False


def parse_config(camx_in_path: str = "CAMx.in") -> RunConfig:
    """Parse CAMx.in and return a RunConfig."""
    parser = NamelistParser()
    groups = parser.parse(camx_in_path)
    return map_namelist_to_config(groups)


def run_model(settings: DriverSettings) -> None:
    """Top-level driver: parse config, orchestrate I/O, run time loop."""
    camx_in = Path(settings.camx_in_path)
    if not camx_in.exists():
        raise FileNotFoundError(f"CAMx.in not found: {camx_in}")

    # Parse configuration
    config = parse_config(str(camx_in))

    # Wire I/O helpers
    reader = NetCDFReader()
    writer = CamxNetCDFWriter(config, chunking=None)

    # Allocate in-memory state
    state = ModelState.allocate(config.grids, len(config.species))

    # Open input datasets if provided
    met2d_ds = [reader.open(p).dataset for p in config.met2d_paths]
    met3d_ds = [reader.open(p).dataset for p in config.met3d_paths]
    cloud_ds = [reader.open(p).dataset for p in config.cloud_paths]
    kv_ds = [reader.open(p).dataset for p in config.kv_paths]
    ic_ds = [reader.open(p).dataset for p in config.ic_paths]
    bnd_ds = [reader.open(p).dataset for p in config.bnd_paths]
    top_ds = [reader.open(p).dataset for p in config.top_paths]

    # Open output datasets
    out_datasets = []
    output_root = Path(settings.output_root)
    if output_root.suffix:
        output_root = output_root.with_suffix("")
    if output_root.exists() and output_root.is_dir():
        base_dir = output_root
        base_prefix = "camx_py_output"
    else:
        base_dir = output_root.parent if output_root.parent != Path("") else Path(".")
        base_prefix = output_root.name
    base_dir.mkdir(parents=True, exist_ok=True)
    for idx, grid in enumerate(config.grids):
        out_path = base_dir / f"{base_prefix}.grid{idx+1}.nc"
        ds = writer.create_file(out_path, grid, is_sample=False, var_dim_is_layers=config.output_3d)
        # Placeholder lat/lon/topo until landuse/geo inputs are wired
        lon = np.zeros((grid.ncol, grid.nrow), dtype=np.float32)
        lat = np.zeros((grid.ncol, grid.nrow), dtype=np.float32)
        topo = np.zeros((grid.ncol, grid.nrow), dtype=np.float32)
        writer.write_grid_data(ds, grid, lon=lon, lat=lat, topo=topo)
        out_datasets.append(ds)

    # Initialize concentrations from IC if provided (master grid only for now)
    if ic_ds:
        ic_data = reader.read_initial_conditions(
            ic_ds[0],
            tuple(s.name for s in config.species),
            config.begin_date,
            config.begin_time,
        )
        grid_state = state.grids[0]
        for sp_idx, sp in enumerate(config.species):
            key = sp.name.strip()
            if key in ic_data:
                grid_state.conc[:, :, :, sp_idx] = ic_data[key]

    def _apply_boundaries(grid_state, bnd_data):
        conc = grid_state.conc
        for sp_idx, sp in enumerate(config.species):
            key = sp.name.strip()
            if key not in bnd_data:
                continue
            arr = bnd_data[key]
            conc[0, 1:-1, :, sp_idx] = arr[0, 1:-1, :]
            conc[-1, 1:-1, :, sp_idx] = arr[-1, 1:-1, :]
            conc[1:-1, 0, :, sp_idx] = arr[1:-1, 0, :]
            conc[1:-1, -1, :, sp_idx] = arr[1:-1, -1, :]

    def _apply_top(grid_state, top_data):
        conc = grid_state.conc
        top_layer = conc.shape[2] - 1
        for sp_idx, sp in enumerate(config.species):
            key = sp.name.strip()
            if key not in top_data:
                continue
            arr = top_data[key]
            conc[:, :, top_layer, sp_idx] = arr[:, :, 0]

    date = config.begin_date
    time_hhmm = config.begin_time
    dtinp = config.dtinp_minutes
    dtout = config.dtout_minutes
    step_minutes = min(dtinp, dtout)

    next_met_date, next_met_time = date, time_hhmm
    next_out_date, next_out_time = date, time_hhmm

    def _update_field(curr, nxt, rate):
        rate[:] = timrates(curr, nxt, dtinp)
        curr[:] = nxt

    def _flag(date_val: int, time_val: float) -> tuple[int, int]:
        century = config.start_year // 100 if config.start_year else 0
        return century * 100000 + date_val, int(time_val * 100)

    out_step_index = 0

    while is_before_or_equal(date, time_hhmm, config.end_date, config.end_time):
        if is_after_or_equal(date, time_hhmm, next_met_date, next_met_time):
            for idx, grid_state in enumerate(state.grids):
                if idx < len(met2d_ds):
                    met2d = reader.read_met_2d(met2d_ds[idx], date, time_hhmm, dtinp)
                    grid_state.tsurf_next[:] = met2d["sfctemperature"]
                    grid_state.tsurf_rate[:] = timrates(grid_state.tsurf, grid_state.tsurf_next, dtinp)
                    grid_state.tsurf[:] = grid_state.tsurf_next
                    snow_next = met2d["snowewd"]
                    grid_state.snowrat[:] = timrates(grid_state.snow, snow_next, dtinp)
                    grid_state.snow[:] = snow_next
                    grid_state.snowage[:] = met2d["snowage"]

                if idx < len(met3d_ds):
                    met3d = reader.read_met_3d(met3d_ds[idx], date, time_hhmm, dtinp)
                    grid_state.height_next[:] = met3d["z"]
                    _update_field(grid_state.height, grid_state.height_next, grid_state.height_rate)
                    grid_state.press_next[:] = met3d["pressure"]
                    _update_field(grid_state.press, grid_state.press_next, grid_state.press_rate)
                    grid_state.temp_next[:] = met3d["temperature"]
                    _update_field(grid_state.tempk, grid_state.temp_next, grid_state.temp_rate)
                    grid_state.water_next[:] = met3d["humidity"]
                    _update_field(grid_state.water, grid_state.water_next, grid_state.water_rate)
                    grid_state.windu_next[:] = met3d["uwind"]
                    _update_field(grid_state.windu, grid_state.windu_next, grid_state.windu_rate)
                    grid_state.windv_next[:] = met3d["vwind"]
                    _update_field(grid_state.windv, grid_state.windv_next, grid_state.windv_rate)

                if idx < len(cloud_ds):
                    cloud = reader.read_cloud(cloud_ds[idx], date, time_hhmm)
                    grid_state.cldwtr[:] = cloud["cloudwater"]
                    grid_state.ranwtr[:] = cloud["rainwater"]
                    grid_state.snowtr[:] = cloud["snowater"]
                    grid_state.gplwtr[:] = cloud["grplwater"]
                    grid_state.cldod[:] = cloud["cloudod"]
                    grid_state.cigfrc[:] = cloud["kf_cldfrac"]
                    grid_state.cigtim[:] = cloud["kf_tscale"]
                    grid_state.cigwtr[:] = cloud["kf_cldwater"]
                    grid_state.cigpcp[:] = cloud["kf_pcpwater"]
                    grid_state.cigent[:] = cloud["kf_entrain"]
                    grid_state.cigdet[:] = cloud["kf_detrain"]

                if idx < len(kv_ds):
                    kv = reader.read_kv(kv_ds[idx], date, time_hhmm, dtinp)
                    grid_state.rkv_next[:] = kv
                    _update_field(grid_state.rkv, grid_state.rkv_next, grid_state.rkv_rate)

                if idx < len(bnd_ds):
                    bnd_data = reader.read_boundary_conditions(
                        bnd_ds[idx], tuple(s.name for s in config.species), date, time_hhmm
                    )
                    _apply_boundaries(grid_state, bnd_data)

                if idx < len(top_ds):
                    top_data = reader.read_top_conditions(
                        top_ds[idx], tuple(s.name for s in config.species), date, time_hhmm
                    )
                    _apply_top(grid_state, top_data)

            next_met_date, next_met_time = add_minutes(next_met_date, next_met_time, dtinp)

        # TODO: transport/chemistry kernels will run here

        if is_after_or_equal(date, time_hhmm, next_out_date, next_out_time):
            end_date, end_time = add_minutes(date, time_hhmm, dtout)
            if not is_before_or_equal(end_date, end_time, config.end_date, config.end_time):
                end_date, end_time = config.end_date, config.end_time
            start_flag = _flag(date, time_hhmm)
            end_flag = _flag(end_date, end_time)
            for idx, grid_state in enumerate(state.grids):
                ds = out_datasets[idx]
                writer.write_tflags(ds, out_step_index, start_flag, end_flag)
                writer.write_species(
                    ds,
                    grid_state.grid,
                    out_step_index,
                    grid_state.height,
                    [grid_state.conc[:, :, :, i] for i in range(grid_state.conc.shape[3])],
                    [True] * grid_state.conc.shape[3],
                    num_dims=4,
                )
            next_out_date, next_out_time = add_minutes(next_out_date, next_out_time, dtout)
            out_step_index += 1

        date, time_hhmm = add_minutes(date, time_hhmm, step_minutes)

    # TODO: main timestep loop, transport/chemistry, outputs
    raise NotImplementedError("Main model loop not yet ported.")


__all__ = ["DriverSettings", "parse_config", "run_model"]
