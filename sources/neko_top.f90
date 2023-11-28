module neko_top
   use neko
   use develop
   use topology_module
   use json_module
   use initial_conditions

   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

   implicit none

   integer :: max_iter = 4
   real(kind=rp) :: perm = 1.0_rp

contains

   !> Assign user conditions for the neko case
   subroutine set_user_case(neko_case)
      use json_utils, only: json_get, json_get_or_default
      implicit none

      type(case_t), intent(inout) :: neko_case
      character(len=:), allocatable :: json_setting

      ! Set the properties for the fluid
      neko_case%usr%material_properties => set_material_properties

      ! Set the permeability forcing, used to enforce the topology
      neko_case%usr%fluid_user_f_vector => topopt_permeability_force

      ! Set the initial conditions
      neko_case%usr%scalar_user_ic => scalar_z_split_ic

   end subroutine set_user_case

   !> Register user defined functions (see nekos user_intf.f90)
   subroutine neko_top_init(neko_case, topology)
      type(case_t), intent(inout) :: neko_case
      type(topology_t), intent(inout) :: topology

      integer, dimension(3) :: ncells

      call set_user_case(neko_case)

      ! Initialize the neko solver
      call neko_init(neko_case)
      call neko_log%init()

      ! Initialize the topology
      ncells = 5
      call topology%init(neko_case, ncells)

   end subroutine neko_top_init

   subroutine neko_top_solve(neko_case, topology)
      type(case_t), intent(inout) :: neko_case
      type(topology_t), intent(inout) :: topology

      integer :: iter

      do iter = 1, max_iter
         call setup_iteration(neko_case, iter)

         call neko_solve(neko_case)

         ! Update the permeability
         call topology%update(neko_case)
      end do
   end subroutine neko_top_solve

   subroutine neko_top_finalize(neko_case, topology)
      type(case_t), intent(inout) :: neko_case
      type(topology_t), intent(inout) :: topology

      logical :: temperature_enabled

      ! ---------------------------------------------------------------------- !
      ! Compute the outlet area-weighted average temperature
      ! ---------------------------------------------------------------------- !
      ! Read the case file for options
      call json_get_or_default(neko_case%params, 'case.scalar.enabled', &
                               temperature_enabled, .false.)

      if (temperature_enabled) then
         call estimate_temperature(neko_case)
      end if

      call neko_finalize(neko_case)
      call topology%free()
      call neko_log%end()

   end subroutine neko_top_finalize

   subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
      real(kind=rp), intent(in) :: t
      integer, intent(in) :: tstep
      real(kind=rp), intent(inout) :: rho, mu, cp, lambda
      type(json_file), intent(inout) :: params

      real(kind=rp) :: Re, Pe

      call json_get(params, 'case.fluid.Re', Re)
      call json_get(params, 'case.fluid.perm', perm)
      call json_get(params, 'case.scalar.Pe', Pe)

      mu = 1.0_rp/Re
      lambda = 1.0_rp/Pe
      rho = 1.0_rp
      cp = 1.0_rp
   end subroutine set_material_properties

   subroutine setup_iteration(neko_case, iter)
      type(case_t), intent(inout) :: neko_case
      integer, intent(in) :: iter
      character(len=:), allocatable :: dirname
      character(len=80) :: file_name
      type(chkp_output_t) :: initial_checkpoint
      type(sampler_t) :: checkpoint_writer

      if (iter .eq. 1) then
         ! Save the initial configuration
         initial_checkpoint = chkp_output_t(neko_case%fluid%chkp, path='./')
         call checkpoint_writer%init(neko_case%end_time)
         call checkpoint_writer%add(initial_checkpoint, 1.0_rp, 'tsteps')
         call checkpoint_writer%sample(0.0_rp, 0, .true.)
      else
         ! Or load the initial configuration
         call restart(neko_case, 'fluid00000.chkp')
      end if

      call json_get_or_default(neko_case%params, 'case.output_directory', dirname, './')

      write (file_name, '(a,a,i5.5,a)') &
         trim(adjustl(dirname)), '/topopt_', iter, '_.fld'

      neko_case%f_out%output_t%file_%file_type%fname = trim(file_name)
      neko_case%f_out%output_t%file_%file_type%counter = 0
      neko_case%f_out%output_t%file_%file_type%start_counter = 0
      call neko_case%s%sample(0.0_rp, 0, .true.)

   end subroutine setup_iteration

   !> Restart a case @a neko_case from a given checkpoint
   subroutine restart(neko_case, restart_file)
      implicit none
      character(len=*), intent(in) :: restart_file
      type(case_t), intent(inout) :: neko_case
      real(kind=rp):: t
      integer :: i
      type(file_t) :: chkpf
      character(len=LOG_SIZE) :: log_buf

      chkpf = file_t(trim(restart_file))
      call chkpf%read(neko_case%fluid%chkp)
      !Free the previous mesh, dont need it anymore
      call neko_case%fluid%chkp%previous_mesh%free()

      ! Make sure that continuity is maintained (important for interpolation)
      call col2(neko_case%fluid%u%x, neko_case%fluid%c_Xh%mult, neko_case%fluid%u%dof%size())
      call col2(neko_case%fluid%v%x, neko_case%fluid%c_Xh%mult, neko_case%fluid%u%dof%size())
      call col2(neko_case%fluid%w%x, neko_case%fluid%c_Xh%mult, neko_case%fluid%u%dof%size())
      call col2(neko_case%fluid%p%x, neko_case%fluid%c_Xh%mult, neko_case%fluid%u%dof%size())
      select type (fld => neko_case%fluid)
      type is (fluid_pnpn_t)
         do i = 1, fld%ulag%size()
            call col2(fld%ulag%lf(i)%x, fld%c_Xh%mult, fld%u%dof%size())
            call col2(fld%vlag%lf(i)%x, fld%c_Xh%mult, fld%u%dof%size())
            call col2(fld%wlag%lf(i)%x, fld%c_Xh%mult, fld%u%dof%size())
         end do
      end select
      if (allocated(neko_case%scalar)) then
         call col2(neko_case%scalar%s%x, neko_case%scalar%c_Xh%mult, neko_case%scalar%s%dof%size())
      end if

      call neko_case%fluid%chkp%sync_device()
      call neko_case%fluid%gs_Xh%op(neko_case%fluid%u, GS_OP_ADD)
      call neko_case%fluid%gs_Xh%op(neko_case%fluid%v, GS_OP_ADD)
      call neko_case%fluid%gs_Xh%op(neko_case%fluid%w, GS_OP_ADD)
      call neko_case%fluid%gs_Xh%op(neko_case%fluid%p, GS_OP_ADD)
      select type (fld => neko_case%fluid)
      type is (fluid_pnpn_t)
         do i = 1, fld%ulag%size()
            call fld%gs_Xh%op(fld%ulag%lf(i), GS_OP_ADD)
            call fld%gs_Xh%op(fld%vlag%lf(i), GS_OP_ADD)
            call fld%gs_Xh%op(fld%wlag%lf(i), GS_OP_ADD)
         end do
      end select

      if (allocated(neko_case%scalar)) then
         call neko_case%scalar%gs_Xh%op(neko_case%scalar%s, GS_OP_ADD)
      end if
      t = neko_case%fluid%chkp%restart_time()
      call neko_log%section('Restarting from checkpoint')
      write (log_buf, '(A,A)') 'File :   ', &
         trim(restart_file)
      call neko_log%message(log_buf)
      write (log_buf, '(A,E15.7)') 'Time : ', t
      call neko_log%message(log_buf)
      call neko_log%end_section()

      call neko_case%s%set_counter(t)

   end subroutine restart
end module neko_top

