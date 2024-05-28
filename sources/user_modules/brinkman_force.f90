module brinkman_force
  use num_types, only: rp
  use fluid_user_source_term, only: fluid_user_source_term_t
  implicit none
  private
  public :: brinkman_force_term

contains

  !> @brief Apply the permeability force term
  !! @details This function applies the permeability force term. This is
  !!          done by computing the permeability at each point in the domain
  !!          and then multiplying it by the velocity at that point. This
  !!          is then added to the force term.
  !!          The force field is read from the field registry named 'brinkman'.
  !!
  !! @param[inout] f The force term
  !! @param[in] t The current time
  subroutine brinkman_force_term(f, t)
    use field, only: field_t
    use field_registry, only: neko_field_registry

    use math, only: col3
    use neko_config, only: NEKO_BCKND_DEVICE
    use device_math, only: device_col3

    implicit none

    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    type(field_t), pointer :: u, v, w, force
    integer :: total_size

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    force => neko_field_registry%get_field('brinkman')

    total_size = f%dm%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(f%u_d, u%x_d, force%x_d, total_size)
       call device_col3(f%v_d, v%x_d, force%x_d, total_size)
       call device_col3(f%w_d, w%x_d, force%x_d, total_size)
    else
       call col3(f%u, u%x, force%x, total_size)
       call col3(f%v, v%x, force%x, total_size)
       call col3(f%w, w%x, force%x, total_size)
    end if

  end subroutine brinkman_force_term

end module brinkman_force
