module linop_base
   use LightKrylov, only: abstract_linop_rdp
   use LightKrylov, only: wp => dp
   use case, only: case_t
   use field_math, only: field_copy
   use neko_vector, only: state_vector_t
   use perturb_case, only: perturb_case_t
   use utils, only: neko_error, neko_warning

   implicit none
   private

   type, abstract, extends(abstract_linop_rdp), public :: linop_base_t
      type(case_t), pointer, public :: base_case => null()
      type(perturb_case_t), pointer, public :: lns_case => null()
      type(perturb_case_t), pointer, public :: adjoint_case => null()
      type(state_vector_t), public :: baseflow
      logical, public :: has_baseflow = .false.
      logical, public :: lns_initialized = .false.
      logical, public :: adjoint_initialized = .false.
      logical, public :: if_2d = .false.
      character(len=:), allocatable, public :: base_output_base_fname
      character(len=:), allocatable, public :: lns_output_base_fname
      character(len=:), allocatable, public :: adjoint_output_base_fname
   contains
      private
      procedure, pass(self), public :: init => linop_base_init
      procedure, pass(self), public :: set_baseflow => linop_base_set_baseflow
      procedure, pass(self), public :: free => linop_base_free
   end type linop_base_t

   public :: linop_base_init, linop_base_store_baseflow, &
        linop_base_require_baseflow, linop_base_require_lns, &
        linop_base_require_adjoint, linop_base_free

contains

   subroutine linop_base_init(self, base_case, lns_case, adjoint_case)
      class(linop_base_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case
      type(perturb_case_t), optional, target, intent(inout) :: lns_case
      type(perturb_case_t), optional, target, intent(inout) :: adjoint_case

      self%base_case => base_case
      self%base_output_base_fname = &
           trim(base_case%f_out%file_%get_base_fname())
      self%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)

      self%lns_initialized = .false.
      self%adjoint_initialized = .false.
      nullify(self%lns_case)
      nullify(self%adjoint_case)

      if (present(lns_case)) then
         self%lns_case => lns_case
         self%lns_initialized = .true.
         self%lns_output_base_fname = trim(lns_case%f_out%file_%get_base_fname())
      else
         call neko_warning('Stability linop initialized without an LNS case')
      end if

      if (present(adjoint_case)) then
         self%adjoint_case => adjoint_case
         self%adjoint_initialized = .true.
         self%adjoint_output_base_fname = &
              trim(adjoint_case%f_out%file_%get_base_fname())
      else
         call neko_warning('Stability linop initialized without an adjoint case')
      end if

      call capture_current_baseflow(self)
   end subroutine linop_base_init

   subroutine linop_base_set_baseflow(self, baseflow, time_interval, &
        fixed_step_count)
      class(linop_base_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: baseflow
      real(wp), optional, intent(in) :: time_interval
      integer, optional, intent(in) :: fixed_step_count

      call linop_base_store_baseflow(self, baseflow)
   end subroutine linop_base_set_baseflow

   subroutine linop_base_store_baseflow(self, baseflow)
      class(linop_base_t), intent(inout) :: self
      type(state_vector_t), intent(in) :: baseflow

      call self%baseflow%init_like(baseflow)
      call self%baseflow%axpby(1.0_wp, baseflow, 0.0_wp)
      self%has_baseflow = .true.
   end subroutine linop_base_store_baseflow

   subroutine linop_base_require_baseflow(self)
      class(linop_base_t), intent(in) :: self

      if (.not. self%has_baseflow) then
         call neko_error('Stability linop baseflow has not been set')
      end if
   end subroutine linop_base_require_baseflow

   subroutine linop_base_require_lns(self)
      class(linop_base_t), intent(in) :: self

      if (.not. self%lns_initialized .or. .not. associated(self%lns_case)) then
         call neko_error('Stability linop LNS case has not been initialized')
      end if
   end subroutine linop_base_require_lns

   subroutine linop_base_require_adjoint(self)
      class(linop_base_t), intent(in) :: self

      if (.not. self%adjoint_initialized .or. &
           .not. associated(self%adjoint_case)) then
         call neko_error('Stability linop adjoint case has not been initialized')
      end if
   end subroutine linop_base_require_adjoint

   subroutine capture_current_baseflow(self)
      class(linop_base_t), intent(inout) :: self
      type(state_vector_t) :: mold

      mold%coef => self%base_case%fluid%c_Xh
      mold%if_2d = self%if_2d
      call self%baseflow%init_like(mold)
      call field_copy(self%baseflow%u, self%base_case%fluid%u)
      call field_copy(self%baseflow%v, self%base_case%fluid%v)
      call field_copy(self%baseflow%w, self%base_case%fluid%w)
      self%has_baseflow = .true.
   end subroutine capture_current_baseflow

   subroutine linop_base_free(self)
      class(linop_base_t), intent(inout) :: self

      if (self%has_baseflow) call self%baseflow%free()
      nullify(self%base_case)
      nullify(self%lns_case)
      nullify(self%adjoint_case)
      self%if_2d = .false.
      self%has_baseflow = .false.
      self%lns_initialized = .false.
      self%adjoint_initialized = .false.
      if (allocated(self%base_output_base_fname)) then
         deallocate(self%base_output_base_fname)
      end if
      if (allocated(self%lns_output_base_fname)) then
         deallocate(self%lns_output_base_fname)
      end if
      if (allocated(self%adjoint_output_base_fname)) then
         deallocate(self%adjoint_output_base_fname)
      end if
   end subroutine linop_base_free

end module linop_base
