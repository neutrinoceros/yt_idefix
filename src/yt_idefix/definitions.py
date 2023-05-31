from typing import Any, Dict

import numpy as np

# Physical constants in c.g.s units defined in $PLUTO_DIR/Src/pluto.h
# The values are copied from PLUTO 4.4-patch-2
pluto_def_constants = {
    "CONST_AH": 1.008,
    "CONST_AHe": 4.004,
    "CONST_AZ": 30.0,
    "CONST_amu": 1.66053886e-24,
    "CONST_au": 1.49597892e13,
    "CONST_c": 2.99792458e10,
    "CONST_e": 4.80320425e-10,
    "CONST_eV": 1.602176463158e-12,
    "CONST_G": 6.6726e-8,
    "CONST_h": 6.62606876e-27,
    "CONST_kB": 1.3806505e-16,
    "CONST_ly": 0.9461e18,
    "CONST_mp": 1.67262171e-24,
    "CONST_mn": 1.67492728e-24,
    "CONST_me": 9.1093826e-28,
    "CONST_mH": 1.6733e-24,
    "CONST_Msun": 2.0e33,
    "CONST_Mearth": 5.9736e27,
    "CONST_NA": 6.0221367e23,
    "CONST_pc": 3.0856775807e18,
    "CONST_PI": 3.14159265358979,
    "CONST_Rearth": 6.378136e8,
    "CONST_Rgas": 8.3144598e7,
    "CONST_Rsun": 6.96e10,
    "CONST_sigma": 5.67051e-5,
    "CONST_sigmaT": 6.6524e-25,
}


class _PlutoBaseUnits:
    """Derive base units from a given combination of units"""

    _base_unit_list = ("mass_unit", "length_unit", "time_unit")

    def __init__(self, unit_combination: Dict[str, Any]):
        # Fow now, unit_combination has been validated before passed in
        # But we still need to check the number of units here for insurance
        if len(unit_combination) != 3:
            raise ValueError(
                "_PlutoBaseUnits requires a combination of 3 units "
                f"({len(unit_combination)} given)"
            )
        self.unit_combination = unit_combination
        self._unit_cache = self.unit_combination
        _data = {}
        for unit in self._base_unit_list:
            self._setup_unit(unit)
            _data[unit] = self._unit_cache[unit]
        self._data = _data

    def __getitem__(self, key):
        return self._data[key]

    # Some unit_setup_functions
    # The condition statements are essential to prevent getting stuck in infinite recursion
    # However they are fragile, one must be cautious to change them
    def _setup_time_unit(self):
        uc = self._unit_cache
        if "time_unit" in uc:
            return
        time = uc["length_unit"] / uc["velocity_unit"]
        uc["time_unit"] = time

    def _setup_length_unit(self):
        uc = self._unit_cache
        if "length_unit" in uc:
            return
        if "time_unit" in uc:
            length = uc["velocity_unit"] * uc["time_unit"]
        else:
            length = (uc["mass_unit"] / uc["density_unit"]) ** (1 / 3)
        uc["length_unit"] = length

    def _setup_mass_unit(self):
        uc = self._unit_cache
        if "mass_unit" in uc:
            return
        if "density_unit" in uc:
            mass = uc["density_unit"] * uc["length_unit"] ** 3
        else:
            mass = (
                (uc["magnetic_unit"] / uc["velocity_unit"]) ** 2
                / 4.0
                / np.pi
                * uc["length_unit"] ** 3
            )
        uc["mass_unit"] = mass

    def _setup_velocity_unit(self):
        uc = self._unit_cache
        if "velocity_unit" in uc:
            return
        if "length_unit" in uc and "time_unit" in uc:
            vel = uc["length_unit"] / uc["time_unit"]
        elif "magnetic_unit" in uc:
            vel = uc["magnetic_unit"] / np.sqrt(4.0 * np.pi * uc["density_unit"])
        else:
            vel = (uc["mass_unit"] / uc["density_unit"]) ** (1 / 3) / uc["time_unit"]
        uc["velocity_unit"] = vel

    def _setup_density_unit(self):
        uc = self._unit_cache
        if "density_unit" in uc:
            return
        if "length_unit" in uc:
            density = uc["mass_unit"] / uc["length_unit"] ** 3
        elif "velocity_unit" in uc:
            density = (uc["magnetic_unit"] / uc["velocity_unit"]) ** 2 / 4.0 / np.pi
        else:
            density = (
                (uc["magnetic_unit"] * uc["time_unit"]) ** 3 / uc["mass_unit"]
            ) ** 2 / (4.0 * np.pi) ** 3
        uc["density_unit"] = density

    def _setup_unit(self, unit):
        """Calculate a unit according to the given unit setup functions"""
        f_name = "_setup_" + unit
        f = getattr(self, f_name)
        try:
            f()
        except KeyError as err:
            # If any unit in the function is missing, calculate that first,
            # then come back and try it again.
            missing_unit = err.args[0]
            self._setup_unit(missing_unit)
            self._setup_unit(unit)
