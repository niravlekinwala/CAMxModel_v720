"""Map CAMx.in namelist groups into RunConfig and GridSpec."""
from __future__ import annotations

from datetime import date as _date
from typing import Dict, List

from camx_py.state.grid import GridSpec
from camx_py.state.run_config import RunConfig
from camx_py.state.species import SpeciesDef


def _parse_int(val: str, default: int = 0) -> int:
    try:
        return int(float(val.strip()))
    except Exception:
        return default


def _parse_float(val: str, default: float = 0.0) -> float:
    try:
        return float(val.strip())
    except Exception:
        return default


def _parse_bool(val: str) -> bool:
    return val.strip().upper() in {"T", "TRUE", ".T.", ".TRUE.", "Y", "YES", "1"}


def _parse_list(val: str) -> List[str]:
    # Handle comma-separated, strip quotes/spaces
    parts = val.replace("'", "").replace('"', "").split(",")
    return [p.strip() for p in parts if p.strip()]


def _normalize_year(year: int) -> int:
    if year < 100:
        return year + 2000
    return year


def _ymd_to_yyjjj(year: int, month: int, day: int) -> int:
    year = _normalize_year(year)
    doy = _date(year, month, day).timetuple().tm_yday
    return (year % 100) * 1000 + doy


def map_namelist_to_config(groups: Dict[str, Dict[str, str]]) -> RunConfig:
    ctl = groups.get("CAMx_Control", {})

    start_list = _parse_list(str(ctl.get("Start_Date_Hour", "")))
    end_list = _parse_list(str(ctl.get("End_Date_Hour", "")))
    if len(start_list) >= 4:
        start_year = _normalize_year(int(start_list[0]))
        begin_date = _ymd_to_yyjjj(start_year, int(start_list[1]), int(start_list[2]))
        begin_time = float(int(start_list[3]))
    else:
        begin_date = 0
        begin_time = 0.0
        start_year = 0
    if len(end_list) >= 4:
        end_year = _normalize_year(int(end_list[0]))
        end_date = _ymd_to_yyjjj(end_year, int(end_list[1]), int(end_list[2]))
        end_time = float(int(end_list[3]))
    else:
        end_date = 0
        end_time = 0.0
        end_year = 0

    dtout = _parse_float(ctl.get("Output_Frequency", "60.0"))
    dtinp = _parse_float(ctl.get("Met_Input_Frequency", "60.0"))
    ncol = _parse_int(ctl.get("Master_Grid_Columns", "1"))
    nrow = _parse_int(ctl.get("Master_Grid_Rows", "1"))
    nlays = _parse_int(ctl.get("Number_of_Layers", "1"))
    xorig = _parse_float(ctl.get("Master_SW_XCoord", "0.0"))
    yorig = _parse_float(ctl.get("Master_SW_YCoord", "0.0"))
    delx = _parse_float(ctl.get("Master_Cell_XSize", "1.0"))
    dely = _parse_float(ctl.get("Master_Cell_YSize", "1.0"))
    map_proj = ctl.get("Map_Projection", "LATLON").strip().upper()
    proj_map = {
        "LATLON": "latlon",
        "UTM": "utm",
        "LAMBERT": "lambert",
        "ROTMERC": "rotated_polar",
        "POLAR": "polar",
        "MERCATOR": "mercator",
    }
    projection = proj_map.get(map_proj, "latlon")
    polelon = _parse_float(ctl.get("Longitude_Pole", "0.0"))
    polelat = _parse_float(ctl.get("Latitude_Pole", "0.0"))
    tlat1 = _parse_float(ctl.get("True_Latitude1", "0.0"))
    tlat2 = _parse_float(ctl.get("True_Latitude2", "0.0"))

    species_names = _parse_list(str(ctl.get("Output_Species_Names", "")))
    species: List[SpeciesDef] = [SpeciesDef(name=s, units="", long_name=s, description=s) for s in species_names]

    met2d_paths = _parse_list(str(ctl.get("Met2D_Grid", "")))
    met3d_paths = _parse_list(str(ctl.get("Met3D_Grid", "")))
    cloud_paths = _parse_list(str(ctl.get("Cloud_Grid", "")))
    kv_paths = _parse_list(str(ctl.get("Vdiff_Grid", "")))
    ic_paths = _parse_list(str(ctl.get("Initial_Conditions", "")))
    bnd_paths = _parse_list(str(ctl.get("Boundary_Conditions", "")))
    top_paths = _parse_list(str(ctl.get("Top_Concentrations", "")))

    grid = GridSpec(
        ncol=ncol,
        nrow=nrow,
        nlays=nlays,
        xorig_km=xorig,
        yorig_km=yorig,
        delx_km=delx,
        dely_km=dely,
        projection=projection,
        polelon=polelon,
        polelat=polelat,
        tlat1=tlat1,
        tlat2=tlat2,
        mesh_factor=1,
        nested_bounds=None,
        staggered_winds=False,
    )

    return RunConfig(
        begin_date=begin_date,
        begin_time=begin_time,
        end_date=end_date,
        end_time=end_time,
        start_year=start_year,
        end_year=end_year,
        dtout_minutes=dtout,
        dtinp_minutes=dtinp,
        grids=[grid],
        species=species,
        run_message=ctl.get("Run_Message", ""),
        version_string="CAMx v7.20",
        time_zone=_parse_int(ctl.get("Time_Zone", "0")),
        probing_tool=str(ctl.get("Probing_Tool", "")).strip(),
        chem_solver=str(ctl.get("Chemistry_Solver", "")).strip(),
        advection_solver=str(ctl.get("Advection_Solver", "")).strip(),
        drydep_model=str(ctl.get("Drydep_Model", "")).strip(),
        wetdep=_parse_bool(ctl.get("Wet_Deposition", "F")),
        acm2=_parse_bool(ctl.get("ACM2_Diffusion", "F")),
        cig_model=_parse_bool(ctl.get("Subgrid_Convection", "F")),
        surface_model=_parse_bool(ctl.get("Surface_Model", "F")),
        inline_ix_emiss=_parse_bool(ctl.get("Inline_Ix_Emissions", "F")),
        bidi_nh3_drydep=_parse_bool(ctl.get("Bidi_NH3_Drydep", "F")),
        strat_o3_profile=_parse_bool(ctl.get("Strat_Ozone_Profile", "F")),
        super_stepping=_parse_bool(ctl.get("Super_Stepping", "F")),
        gridded_emiss=_parse_bool(ctl.get("Gridded_Emissions", "F")),
        point_emiss=_parse_bool(ctl.get("Point_Emissions", "F")),
        ignore_emiss_dates=_parse_bool(ctl.get("Ignore_Emission_Dates", "F")),
        output_3d=_parse_bool(ctl.get("Output_3D_Grid", "T")),
        pig_sample_grid=_parse_bool(ctl.get("PiG_Sampling_Grid", "F")),
        pig_sample_bckgnd=_parse_bool(ctl.get("Sample_Background", "F")),
        met2d_paths=met2d_paths,
        met3d_paths=met3d_paths,
        cloud_paths=cloud_paths,
        kv_paths=kv_paths,
        ic_paths=ic_paths,
        bnd_paths=bnd_paths,
        top_paths=top_paths,
    )
