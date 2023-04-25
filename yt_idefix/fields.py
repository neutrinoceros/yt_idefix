import re

from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases


class BaseFields(FieldInfoContainer):
    # This class can be used later to add functions and attributes common to all fields
    pass


class IdefixVtkFields(BaseFields):
    known_other_fields = (
        ("RHO", ("code_mass / code_length**3", ["density"], None)),  # type: ignore
        ("VX1", ("code_length / code_time", ["velocity_x"], None)),
        ("VX2", ("code_length / code_time", ["velocity_y"], None)),
        ("VX3", ("code_length / code_time", ["velocity_z"], None)),
        ("BX1", ("code_magnetic", [], None)),
        ("BX2", ("code_magnetic", [], None)),
        ("BX3", ("code_magnetic", [], None)),
        ("PRS", ("code_pressure", ["pressure"], None)),
    )

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "idefix-vtk", [f"BX{idir}" for idir in "123"]
        )


class IdefixDmpFields(BaseFields):
    known_other_fields = (
        ("Vc-RHO", ("code_mass / code_length**3", ["density"], None)),  # type: ignore
        ("Vc-VX1", ("code_length / code_time", ["velocity_x"], None)),
        ("Vc-VX2", ("code_length / code_time", ["velocity_y"], None)),
        ("Vc-VX3", ("code_length / code_time", ["velocity_z"], None)),
        ("Vc-BX1", ("code_magnetic", [], None)),
        ("Vc-BX2", ("code_magnetic", [], None)),
        ("Vc-BX3", ("code_magnetic", [], None)),
        ("Vc-PRS", ("code_pressure", ["pressure"], None)),
    )
    # note that velocity '_x', '_y' and '_z' aliases are meant to be
    # overwriten according to geometry in self.setup_fluid_aliases

    known_particle_fields = ()

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "idefix-dmp", [f"Vc-BX{idir}" for idir in "123"]
        )


class PlutoFields(BaseFields):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("RHO", ("code_mass / code_length**3", ["density", "rho"], None)),
        ("VX1", ("code_length/code_time", ["vel_x"], None)),
        ("VX2", ("code_length/code_time", ["vel_y"], None)),
        ("VX3", ("code_length/code_time", ["vel_z"], None)),
        ("PRS", ("code_pressure", ["prs", "pres", "pressure"], None)),
        ("BX1", ("code_magnetic", [], None)),
        ("BX2", ("code_magnetic", [], None)),
        ("BX3", ("code_magnetic", [], None)),
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, self.ds._dataset_type, [f"Bx{idir}" for idir in "123"]
        )

        # Add tracer fields
        _TRC_REGEXP = re.compile(r"^TR(\d+)")
        for _, field in self.field_list:
            if (tr_match := re.fullmatch(_TRC_REGEXP, field)) is not None:
                idx = tr_match.group(1)
                self.add_output_field(
                    (self.ds._dataset_type, f"tr{idx}"),
                    sampling_type="cell",
                    units="",
                )
                self.alias(
                    ("gas", f"tracer_{idx}"),
                    (self.ds._dataset_type, f"tr{idx}"),
                    units="",
                )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
        # {TODO} Starting with version 4.4 PLUTO has particle support and this function needs an implementation in future
