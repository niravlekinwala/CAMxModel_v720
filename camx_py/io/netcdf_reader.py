"""NetCDF input reader mirroring CAMx v7.20 IO_NCF layout."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

try:
    import netCDF4 as nc
except ImportError:  # pragma: no cover
    nc = None


@dataclass
class NetCDFInput:
    """Container for open NetCDF datasets."""

    path: str
    dataset: "nc.Dataset"


class NetCDFReader:
    """Helper to read CAMx-format NetCDF inputs (met, IC/BC, emissions)."""

    def __init__(self):
        if nc is None:
            raise RuntimeError("netCDF4 is required to read NetCDF inputs")

    def open(self, path: str) -> NetCDFInput:
        ds = nc.Dataset(path, mode="r")
        return NetCDFInput(path=path, dataset=ds)

    def read_tstep_flags(self, ds: "nc.Dataset", tstep_idx: int) -> Tuple[np.ndarray, np.ndarray]:
        """Return TFLAG and ETFLAG arrays for a timestep index."""
        return (
            np.array(ds.variables["TFLAG"][:, :, tstep_idx]),
            np.array(ds.variables["ETFLAG"][:, :, tstep_idx]),
        )

    def find_tstep_index(
        self,
        ds: "nc.Dataset",
        this_date: int,
        this_time: float,
        ignore_date: bool = False,
        strict: bool = False,
    ) -> int:
        """Find a timestep index that spans the requested date/time.

        Mirrors IO_NCF/ncf_get_tstep.f logic.
        """
        tflag = ds.variables["TFLAG"]
        etflag = ds.variables["ETFLAG"]
        # current date/time combined
        date_time = float(this_date) + float(this_time) / 2400.0

        nsteps = tflag.shape[2]
        for idx in range(nsteps):
            t_in = tflag[:, 0, idx]
            e_in = etflag[:, 0, idx]
            if ignore_date:
                date_time_t = float(this_date) + float(t_in[1]) / 240000.0
                date_time_e = float(this_date) + float(e_in[1]) / 240000.0
                if e_in[1] == 0:
                    date_time_e += 1.0
            else:
                tdate = int(t_in[0])
                edate = int(e_in[0])
                # map 2-digit year to 19xx/20xx like Fortran
                if (tdate // 1000) % 100 <= 80:
                    date_time_t = float(tdate - 2000000) + float(t_in[1]) / 240000.0
                    date_time_e = float(edate - 2000000) + float(e_in[1]) / 240000.0
                else:
                    date_time_t = float(tdate - 1900000) + float(t_in[1]) / 240000.0
                    date_time_e = float(edate - 1900000) + float(e_in[1]) / 240000.0

            if strict:
                if date_time >= date_time_t and date_time < date_time_e:
                    return idx
            else:
                if date_time >= date_time_t and date_time <= date_time_e:
                    return idx

        raise ValueError("No suitable timestep found for requested date/time.")

    def read_species_slice(
        self, ds: "nc.Dataset", var_name: str, tstep_idx: int, four_d: bool = True
    ) -> np.ndarray:
        """Read one species variable at a timestep."""
        if four_d:
            return np.array(ds.variables[var_name][:, :, :, tstep_idx])
        return np.array(ds.variables[var_name][:, :, tstep_idx])

    def read_heights(self, ds: "nc.Dataset", tstep_idx: int) -> np.ndarray:
        """Read z field at a timestep."""
        return np.array(ds.variables["z"][:, :, :, tstep_idx])

    def read_var_2d(self, ds: "nc.Dataset", var_name: str) -> np.ndarray:
        """Read full 2D variable."""
        return np.array(ds.variables[var_name][:, :])

    def read_var_3d_tstep(self, ds: "nc.Dataset", var_name: str, tstep_idx: int) -> np.ndarray:
        """Read 3D variable with time (COL, ROW, TSTEP)."""
        return np.array(ds.variables[var_name][:, :, tstep_idx])

    def read_var_4d_tstep(self, ds: "nc.Dataset", var_name: str, tstep_idx: int) -> np.ndarray:
        """Read 4D variable with time (COL, ROW, LAY, TSTEP)."""
        return np.array(ds.variables[var_name][:, :, :, tstep_idx])

    def read_initial_conditions(
        self,
        ds: "nc.Dataset",
        species_names: Tuple[str, ...],
        date: int,
        time_hhmm: float,
    ) -> Dict[str, np.ndarray]:
        """Read initial conditions for provided species at model start time."""
        tstep = self.find_tstep_index(ds, date, time_hhmm, ignore_date=False, strict=True)
        data: Dict[str, np.ndarray] = {}
        for name in species_names:
            key = name.strip()
            if key in ds.variables:
                data[key] = self.read_var_4d_tstep(ds, key, tstep)
        return data

    def read_boundary_conditions(
        self,
        ds: "nc.Dataset",
        species_names: Tuple[str, ...],
        date: int,
        time_hhmm: float,
    ) -> Dict[str, np.ndarray]:
        """Read lateral boundary conditions for provided species."""
        tstep = self.find_tstep_index(ds, date, time_hhmm, ignore_date=False, strict=True)
        data: Dict[str, np.ndarray] = {}
        for name in species_names:
            key = name.strip()
            if key in ds.variables:
                data[key] = self.read_var_4d_tstep(ds, key, tstep)
        return data

    def read_top_conditions(
        self,
        ds: "nc.Dataset",
        species_names: Tuple[str, ...],
        date: int,
        time_hhmm: float,
    ) -> Dict[str, np.ndarray]:
        """Read top boundary concentrations for provided species (LAY=1)."""
        tstep = self.find_tstep_index(ds, date, time_hhmm, ignore_date=False, strict=True)
        data: Dict[str, np.ndarray] = {}
        for name in species_names:
            key = name.strip()
            if key in ds.variables:
                data[key] = self.read_var_4d_tstep(ds, key, tstep)
        return data

    def read_static_grid(self, ds: "nc.Dataset") -> Dict[str, np.ndarray]:
        """Read static grid variables X, Y, layer, longitude, latitude, topo."""
        return {
            "X": np.array(ds.variables["X"][:]),
            "Y": np.array(ds.variables["Y"][:]),
            "layer": np.array(ds.variables["layer"][:]),
            "longitude": np.array(ds.variables["longitude"][:, :]),
            "latitude": np.array(ds.variables["latitude"][:, :]),
            "topo": np.array(ds.variables["topo"][:, :]),
        }

    def read_met_2d(self, ds: "nc.Dataset", date: int, time_hhmm: float, dt_minutes: float) -> Dict[str, np.ndarray]:
        """Read 2D meteorology fields at the next update time."""
        hdate, htime = self._next_update_time(date, time_hhmm, dt_minutes)
        tstep = self.find_tstep_index(ds, hdate, htime, ignore_date=False, strict=False)
        return {
            "sfctemperature": self.read_var_3d_tstep(ds, "sfctemperature", tstep),
            "snowewd": self._read_optional_3d_tstep(ds, "snowewd", tstep),
            "snowage": self._read_optional_3d_tstep(ds, "snowage", tstep),
        }

    def read_met_3d(self, ds: "nc.Dataset", date: int, time_hhmm: float, dt_minutes: float) -> Dict[str, np.ndarray]:
        """Read 3D meteorology fields at the next update time."""
        hdate, htime = self._next_update_time(date, time_hhmm, dt_minutes)
        tstep = self.find_tstep_index(ds, hdate, htime, ignore_date=False, strict=False)
        return {
            "z": self.read_var_4d_tstep(ds, "z", tstep),
            "pressure": self.read_var_4d_tstep(ds, "pressure", tstep),
            "temperature": self.read_var_4d_tstep(ds, "temperature", tstep),
            "humidity": self.read_var_4d_tstep(ds, "humidity", tstep),
            "uwind": self.read_var_4d_tstep(ds, "uwind", tstep),
            "vwind": self.read_var_4d_tstep(ds, "vwind", tstep),
        }

    def read_cloud(self, ds: "nc.Dataset", date: int, time_hhmm: float) -> Dict[str, np.ndarray]:
        """Read cloud-related fields at the current time."""
        tstep = self.find_tstep_index(ds, date, time_hhmm, ignore_date=False, strict=False)
        return {
            "cloudwater": self._read_optional_4d_tstep(ds, "cloudwater", tstep),
            "rainwater": self._read_optional_4d_tstep(ds, "rainwater", tstep),
            "snowater": self._read_optional_4d_tstep(ds, "snowater", tstep),
            "grplwater": self._read_optional_4d_tstep(ds, "grplwater", tstep),
            "cloudod": self._read_optional_4d_tstep(ds, "cloudod", tstep),
            "kf_cldfrac": self._read_optional_3d_tstep(ds, "kf_cldfrac", tstep),
            "kf_tscale": self._read_optional_3d_tstep(ds, "kf_tscale", tstep),
            "kf_cldwater": self._read_optional_4d_tstep(ds, "kf_cldwater", tstep),
            "kf_pcpwater": self._read_optional_4d_tstep(ds, "kf_pcpwater", tstep),
            "kf_entrain": self._read_optional_4d_tstep(ds, "kf_entrain", tstep),
            "kf_detrain": self._read_optional_4d_tstep(ds, "kf_detrain", tstep),
        }

    def read_kv(self, ds: "nc.Dataset", date: int, time_hhmm: float, dt_minutes: float) -> np.ndarray:
        """Read KV field at next update time."""
        hdate, htime = self._next_update_time(date, time_hhmm, dt_minutes)
        tstep = self.find_tstep_index(ds, hdate, htime, ignore_date=False, strict=False)
        return self.read_var_4d_tstep(ds, "kv", tstep)

    def _read_optional_3d_tstep(self, ds: "nc.Dataset", var_name: str, tstep_idx: int) -> np.ndarray:
        if var_name not in ds.variables:
            shape = (len(ds.dimensions["COL"]), len(ds.dimensions["ROW"]))
            return np.zeros(shape, dtype=np.float32)
        return self.read_var_3d_tstep(ds, var_name, tstep_idx)

    def _read_optional_4d_tstep(self, ds: "nc.Dataset", var_name: str, tstep_idx: int) -> np.ndarray:
        if var_name not in ds.variables:
            shape = (len(ds.dimensions["COL"]), len(ds.dimensions["ROW"]), len(ds.dimensions["LAY"]))
            return np.zeros(shape, dtype=np.float32)
        return self.read_var_4d_tstep(ds, var_name, tstep_idx)

    def _next_update_time(self, date: int, time_hhmm: float, dt_minutes: float) -> Tuple[int, float]:
        """Compute next update time following CAMx dtinp rules."""
        hour = int(time_hhmm // 100)
        minute = int(round(time_hhmm % 100))
        next_minute = minute + int(dt_minutes)
        next_hour = hour + next_minute // 60
        next_minute = next_minute % 60
        next_time = float(next_hour * 100 + next_minute)
        next_date = date
        if next_time >= 2400:
            next_time -= 2400
            next_date = self._add_day(next_date)
        return next_date, next_time

    def _add_day(self, date: int) -> int:
        """Increment YYJJJ date by one day with leap-year handling."""
        year = date // 1000
        jday = date % 1000
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
