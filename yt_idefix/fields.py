from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases


class BaseVtkFields(FieldInfoContainer):
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


class IdefixVtkFields(BaseVtkFields):
    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "idefix-vtk", [f"BX{idir}" for idir in "123"]
        )


class PlutoVtkFields(BaseVtkFields):
    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(self, "pluto-vtk", [f"BX{idir}" for idir in "123"])


class IdefixDmpFields(FieldInfoContainer):
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
    # overwritten according to geometry in self.setup_fluid_aliases

    # note: I'm not sure about the note above anymore
    # TODO: check that velocity reading works in non-cartesian geometries
    # for now I'm making similar assumptions for particles
    known_particle_fields = (
        ("PX1", ("code_length", ["particle_position_x"], None)),
        ("PX2", ("code_length", ["particle_position_y"], None)),
        ("PX3", ("code_length", ["particle_position_z"], None)),
        ("PVX1", ("code_length / code_time", ["particle_velocity_x"], None)),
        ("PVX2", ("code_length / code_time", ["particle_velocity_y"], None)),
        ("PVX3", ("code_length / code_time", ["particle_velocity_z"], None)),
    )

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "idefix-dmp", [f"Vc-BX{idir}" for idir in "123"]
        )
