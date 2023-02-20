from yt.fields.field_info_container import FieldInfoContainer
from yt.fields.magnetic_field import setup_magnetic_field_aliases


class BaseFields(FieldInfoContainer):
    # This class can be used later to add functions and attributes common to all fields
    pass


class BaseVtkFields(BaseFields):
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


class PlutoXdmfFields(BaseFields):
    _pres_units = "code_pressure"
    _erg_units = "code_mass * (code_length/code_time)**2"
    _rho_units = "code_mass / code_length**3"
    _mom_units = "code_mass / code_length**2 / code_time"
    _vel_units = "code_length/code_time"

    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("rho", (_rho_units, ["density", "rho"], None)),
        ("vx1", (_vel_units, ["vel_x"], None)),
        ("vx2", (_vel_units, ["vel_y"], None)),
        ("vx3", (_vel_units, ["vel_z"], None)),
        ("prs", (_pres_units, ["prs", "pres", "pressure"], None)),
        ("Bx1", ("code_magnetic", [], None)),
        ("Bx2", ("code_magnetic", [], None)),
        ("Bx3", ("code_magnetic", [], None)),
    )

    known_particle_fields = ()

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system

        setup_magnetic_field_aliases(
            self, self.ds._dataset_type, [f"Bx{idir}" for idir in "123"]
        )

        # Add tracer fields
        for i in range(1, self.ds._count_tracers() + 1):
            if (self.ds._dataset_type, f"tr{i}") in self.field_list:
                self.add_output_field(
                    (self.ds._dataset_type, f"tr{i}"),
                    sampling_type="cell",
                    units="",
                )
                self.alias(
                    ("gas", f"tracer_{i}"),
                    (self.ds._dataset_type, f"tr{i}"),
                    units="",
                )
        # Add temperature field, these are commonly defined as either "Temp" or "Temperature" by users in PLUTO
        if (self.ds._dataset_type, "Temp") in self.field_list:
            self.add_output_field(
                (self.ds._dataset_type, "Temp"),
                sampling_type="cell",
                units=unit_system["temperature"],
            )
            self.alias(
                ("gas", "temperature"),
                (self.ds._dataset_type, "Temp"),
                units=unit_system["temperature"],
            )
        elif (self.ds._dataset_type, "Temperature") in self.field_list:
            self.add_output_field(
                (self.ds._dataset_type, "Temperature"),
                sampling_type="cell",
                units=unit_system["temperature"],
            )
            self.alias(
                ("gas", "temperature"),
                (self.ds._dataset_type, "Temperature"),
                units=unit_system["temperature"],
            )
        # Optional particle number density field in PLUTO
        if (self.ds._dataset_type, "ndens") in self.field_list:
            self.add_output_field(
                (self.ds._dataset_type, "ndens"),
                sampling_type="cell",
                units=unit_system["length"] ** -3,
            )
            self.alias(
                ("gas", "number_density"),
                (self.ds._dataset_type, "ndens"),
                units=unit_system["length"] ** -3,
            )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
        # {TODO} Starting with version 4.4 PLUTO has particle support and this function needs an implementation in future
