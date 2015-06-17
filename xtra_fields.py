def _Rel_KineticEnergy(field,data):
    bv_x,bv_y,bv_z = data.get_field_parameter("bulk_velocity")
    return 0.5 * (data["Density"] * (
        (data["x-velocity"] - bv_x)**2
        + (data["y-velocity"] - bv_y)**2
        + (data["z-velocity"] - bv_z)**2 ))

add_field("RelKineticEnergy",function = _Rel_KineticEnergy, units=r"\rm{ergs}",
                validators=[ValidateParameter("bulk_velocity")])


