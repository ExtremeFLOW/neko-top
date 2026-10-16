module floquet_output_utils
   use LightKrylov, only: wp => dp
   use field_math, only: field_cmult
   use field_output, only: field_output_t
   use num_types, only: rp
   use perturb_case, only: perturb_case_t

   implicit none
   private

   public :: sample_floquet_mode_output

contains

   subroutine sample_floquet_mode_output(output, perturb_case, time_value, idx, &
        multiplier, period)
      type(field_output_t), intent(inout) :: output
      type(perturb_case_t), intent(inout) :: perturb_case
      real(rp), intent(in) :: time_value
      integer, intent(in) :: idx
      complex(kind=wp), intent(in) :: multiplier
      real(wp), intent(in) :: period
      real(wp) :: scale_factor

      scale_factor = exp(-real(log(multiplier), kind=wp) * &
           real(time_value, kind=wp) / period)

      call scale_perturb_velocity(perturb_case, scale_factor)
      call output%set_counter(idx - 1)
      call output%sample(time_value)
      call scale_perturb_velocity(perturb_case, 1.0_wp / scale_factor)
   end subroutine sample_floquet_mode_output

   subroutine scale_perturb_velocity(perturb_case, scale_factor)
      type(perturb_case_t), intent(inout) :: perturb_case
      real(wp), intent(in) :: scale_factor

      call field_cmult(perturb_case%fluid_pert%u_pert, &
           real(scale_factor, kind=rp))
      call field_cmult(perturb_case%fluid_pert%v_pert, &
           real(scale_factor, kind=rp))
      call field_cmult(perturb_case%fluid_pert%w_pert, &
           real(scale_factor, kind=rp))
   end subroutine scale_perturb_velocity

end module floquet_output_utils
