"""Lightweight parser scaffold for CAMx.in namelist."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict


class NamelistParser:
    """Very small subset parser for CAMx Fortran namelist files."""

    _group_re = re.compile(r"^\s*&(?P<name>[A-Za-z0-9_]+)")
    _assign_re = re.compile(r"^\s*(?P<key>[A-Za-z0-9_]+)\s*=\s*(?P<value>.+?)(,|$)")

    def parse(self, path: str | Path) -> Dict[str, Dict[str, str]]:
        """Parse CAMx.in into a nested dict: group -> key -> raw string value.

        This is intentionally minimal and does not evaluate Fortran syntax.
        """
        content = Path(path).read_text()
        groups: Dict[str, Dict[str, str]] = {}
        current: Dict[str, str] | None = None
        current_name: str | None = None
        for line in content.splitlines():
            if not line.strip() or line.strip().startswith("!"):
                continue
            m = self._group_re.match(line)
            if m:
                current_name = m.group("name")
                current = groups.setdefault(current_name, {})
                continue
            if line.strip().startswith("/"):
                current = None
                current_name = None
                continue
            if current is None:
                continue
            am = self._assign_re.match(line)
            if am:
                key = am.group("key")
                val = am.group("value").strip()
                current[key] = val
        return groups

