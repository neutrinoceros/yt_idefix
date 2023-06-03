from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases


class IdefixVtkFields(FieldInfoContainer):
    known_other_fields = (
        ("RHO", ("code_mass / code_length**3", ["density"], None)),  # type: ignore
        ("VX1", ("code_length / code_time", ["velocity_x"], None)),
        ("VX2", ("code_length / code_time", ["velocity_y"], None)),
        ("VX3", ("code_length / code_time", ["velocity_z"], None)),
        ("BX1", ("code_magnetic", [], None)),
        ("BX2", ("code_magnetic", [], None)),
        ("BX3", ("code_magnetic", [], None)),
        ("PRS", ("code_pressure", ["pressure"], None)),
        (
            "PART_RHO",
            ("code_mass / code_length**3", ["deposited_particle_density"], None),
        ),
        (
            "PART_VX1",
            ("code_length / code_time", ["deposited_particle_velocity_x"], None),
        ),
        (
            "PART_VX2",
            ("code_length / code_time", ["deposited_particle_velocity_y"], None),
        ),
        (
            "PART_VX3",
            ("code_length / code_time", ["deposited_particle_velocity_z"], None),
        ),
    )

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "idefix-vtk", [f"BX{idir}" for idir in "123"]
        )


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
    # overwriten according to geometry in self.setup_fluid_aliases

    known_particle_fields = ()

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, "idefix-dmp", [f"Vc-BX{idir}" for idir in "123"]
        )


class PlutoFields(FieldInfoContainer):
    known_other_fields = (
        # PLUTO 4 standard variable names normalized in upper case,
        # Referring to Tools/IDL/pload.pro in Pluto v4.2-patch2
        # Field names of ion fractions are normalized to upper case
        # Each entry here is of the form
        # ( "name", ("code_unit", ["alias1", "alias2", ...], "display_name")),
        # where alias should be set with yt-standard-name or pluto-macro-name
        ("RHO", ("code_density", ["density"], None)),
        ("VX1", ("code_velocity", ["velocity_x"], None)),
        ("VX2", ("code_velocity", ["velocity_y"], None)),
        ("VX3", ("code_velocity", ["velocity_z"], None)),
        ("PRS", ("code_pressure", ["pressure"], None)),
        ("AX1", ("code_magnetic", [], None)),
        ("AX2", ("code_magnetic", [], None)),
        ("AX3", ("code_magnetic", [], None)),
        ("BX1", ("code_magnetic", [], None)),
        ("BX2", ("code_magnetic", [], None)),
        ("BX3", ("code_magnetic", [], None)),
        ("EX1", ("code_magnetic/c", [], None)),
        ("EX2", ("code_magnetic/c", [], None)),
        ("EX3", ("code_magnetic/c", [], None)),
        ("BX1S", ("code_magnetic", [], None)),
        ("BX2S", ("code_magnetic", [], None)),
        ("BX3S", ("code_magnetic", [], None)),
        ("EX1S", ("code_magnetic/c", [], None)),
        ("EX2S", ("code_magnetic/c", [], None)),
        ("EX3S", ("code_magnetic/c", [], None)),
        ("CRG", ("", ["code_magnetic/code_length/c"], None)),
        ("PSI_GLM", ("", [], None)),
        ("PHI_GLM", ("", [], None)),
        ("RHO_D", ("code_density", [], None)),
        ("VX1_D", ("code_velocity", [], None)),
        ("VX2_D", ("code_velocity", [], None)),
        ("VX3_D", ("code_velocity", [], None)),
        ("ENR", ("code_pressure", [], None)),
        ("FR1", ("code_pressure*code_velocity", [], None)),
        ("FR2", ("code_pressure*code_velocity", [], None)),
        ("FR3", ("code_pressure*code_velocity", [], None)),
        ("ENTROPY", ("", ["entropy"], None)),
        ("TR1", ("", [], None)),
        ("TR2", ("", [], None)),
        ("TR3", ("", [], None)),
        ("TR4", ("", [], None)),
        ("TR5", ("", [], None)),
        ("TR6", ("", [], None)),
        ("TR7", ("", [], None)),
        ("TR8", ("", [], None)),
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        setup_magnetic_field_aliases(
            self, self.ds._dataset_type, [f"BX{idir}" for idir in "123"]
        )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
        # {TODO} Starting with version 4.4 PLUTO has particle support and this function needs an implementation in future
