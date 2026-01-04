"""Species metadata for CAMx Python port."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SpeciesDef:
    """Species definition aligned with CAMx NetCDF metadata."""

    name: str
    units: str
    long_name: str
    description: str
    coordinates: str = "latitude longitude"

