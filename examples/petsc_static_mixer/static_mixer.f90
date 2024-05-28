module user
  use user_intf, only: user_t
  use fluid_user_source_term, only: fluid_user_source_term_t
  use field, only: field_t
  use coefs, only: coef_t
  use json_file_module, only: json_file
  use num_types, only: rp
  use tri, only: tri_t
  use tri_mesh, only: tri_mesh_t
  use brinkman_force, only: brinkman_force_term

  implicit none
  public

contains

  subroutine user_setup(u)
    use initial_conditions, only: scalar_z_split_ic

    type(user_t), intent(inout) :: u

    u%scalar_user_ic => scalar_z_split_ic
    u%fluid_user_f_vector => brinkman_force_term

  end subroutine user_setup

end module user
