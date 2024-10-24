! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> A PDE based filter
module PDE_filter
  use num_types, only: rp
  use json_module, only: json_file
  use field_registry, only: neko_field_registry
  use field, only: field_t
  use coefs, only: coef_t
  use ax_product, only: ax_t, ax_helm_factory
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory, &
       krylov_solver_destroy
  use precon, only: pc_t, precon_factory, precon_destroy
  use bc, only: bc_list_add, bc_list_t, bc_list_apply_scalar, bc_list_init, &
       bc_list_free
  use neumann, only: neumann_t
  use profiler, only: profiler_start_region, profiler_end_region
  use gather_scatter, only: gs_t, GS_OP_ADD
  use pnpn_residual, only: pnpn_prs_res_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use field_registry, only: neko_field_registry
  use mapping, only: mapping_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_copy
  use coefs, only: coef_t
  use logger, only: neko_log, LOG_SIZE
  use neko_config, only: NEKO_BCKND_DEVICE
  use dofmap, only: dofmap_t
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use hsmg, only: hsmg_t
  use utils, only: neko_error
  implicit none
  private

  !> A PDE based filter mapping $\rho \mapsto \tilde{\rho}$,
  !! see Lazarov & O. Sigmund 2010,
  !! by solving an equation
  !! of the form $\f -r^2 \nabla^2 \tilde{\rho} + \tilde{\rho} = \rho \f$
  type, public, extends(mapping_t) :: PDE_filter_t
     !> Ax
     class(ax_t), allocatable :: Ax
     !> Solver results monitors ( filter )
     type(ksp_monitor_t) :: ksp_results(1)
     !> Krylov solver for the filter
     class(ksp_t), allocatable :: ksp_filt
     !> Filter Preconditioner
     class(pc_t), allocatable :: pc_filt
     !> They will all be Neumann conditions.
     type(neumann_t) :: filter_bcs
     !> Filter boundary conditions
     type(bc_list_t) :: bclst_filt

     ! Inputs from the user
     !> filter radius
     real(kind=rp) :: r
     !> tolerance for PDE filter
     real(kind=rp) :: abstol_filt
     !> max iterations for PDE filter
     integer :: ksp_max_iter
     !> method for solving PDE
     character(len=:), allocatable :: ksp_solver
     ! > preconditioner type
     character(len=:), allocatable :: precon_type_filt
     integer :: ksp_n, n, i



   contains
     !> Constructor from json.
     procedure, pass(this) :: init => PDE_filter_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          PDE_filter_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => PDE_filter_free
     !> Apply the filter
     procedure, pass(this) :: apply_forward => PDE_filter_apply
     !> Apply the adjoint filter
     ! TODO
     ! TALK TO NIELS, I think this is correct...
     ! but it's not exactly "chain ruling back"
     ! it's filtering the sensitivity

     ! UPDATE:
     ! After an email with Niels, we should be using the chain rule,
     ! not a sensitivity filter
     procedure, pass(this) :: apply_backward => PDE_filter_apply_backward
  end type PDE_filter_t

contains

  !> Constructor from json.
  subroutine PDE_filter_init_from_json(this, json, coef)
    class(PDE_filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout) :: coef
    character(len=20) :: fields(1)
    real(kind=rp) :: tmp_real


    ! TODO
    ! I'll do the json stuff later...
    this%r = 0.1
    this%abstol_filt = 0.0000000001_rp
    this%ksp_max_iter = 200
    this%ksp_solver = "gmres"
    this%precon_type_filt = "ident"

    call this%init_base(json, coef)
    call PDE_filter_init_from_attributes(this, coef)

  end subroutine PDE_filter_init_from_json

  !> Actual constructor.
  subroutine PDE_filter_init_from_attributes(this, coef)
    class(PDE_filter_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    integer :: n
    character(len=NEKO_MSH_MAX_ZLBL_LEN) :: &
         bc_labels_all_neuman(NEKO_MSH_MAX_ZLBLS)



    n = this%coef%dof%size()

    ! initialize the filter BCs
    call this%filter_bcs%init_base(this%coef)

    ! Create list with just Neumann bcs

    ! init the list
    call bc_list_init(this%bclst_filt)

    ! Mark ALL the BCs as Neumann, regardless of what's prescribed
    bc_labels_all_neuman = 'o'
    call this%filter_bcs%mark_zones_from_list(coef%msh%labeled_zones,&
         'o', bc_labels_all_neuman)

    ! set the flux to zero
    call this%filter_bcs%finalize_neumann(0.0_rp)

    ! add them to the filter BCs
    call bc_list_add(this%bclst_filt, this%filter_bcs)


    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax, full_formulation = .false.)

    ! set up krylov solver
    call krylov_solver_factory(this%ksp_filt, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol_filt)

    ! set up preconditioner
    call filter_precon_factory(this%pc_filt, this%ksp_filt, &
         this%coef, this%coef%dof, this%coef%gs_h, this%bclst_filt, &
         this%precon_type_filt)

  end subroutine PDE_filter_init_from_attributes

  !> Destructor.
  subroutine PDE_filter_free(this)
    class(PDE_filter_t), intent(inout) :: this
    if (allocated(this%Ax)) then
       deallocate(this%Ax)
    end if

    call krylov_solver_destroy(this%ksp_filt)

    call precon_destroy(this%pc_filt)

    call this%filter_bcs%free()

    call bc_list_free(this%bclst_filt)

    call this%free_base()

  end subroutine PDE_filter_free

  !> Apply the filter
  !! @param X_out filtered field
  !! @param X_in unfiltered field
  subroutine PDE_filter_apply(this, X_out, X_in)
    class(PDE_filter_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(inout) :: X_out
    integer :: n, i
    type(field_t), pointer :: RHS
    type(field_t), pointer :: ta1, ta2, ta3
    character(len=LOG_SIZE) :: log_buf
    integer :: temp_indices(1)

    n = this%coef%dof%size()
    call neko_scratch_registry%request_field(RHS, temp_indices(1))

    ! set up Helmholtz operators and RHS
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! TODO
       ! I think this is correct but I've never tested it
       call device_cfill(this%coef%h1_d, this%r**2, n)
       call device_cfill(this%coef%h2_d, 1.0_rp, n)
       call device_col3(RHS%x_d, X_in%x_d, this%coef%B_d, n)
    else
       do i = 1, n
          ! h1 is already negative in its definition
          this%coef%h1(i,1,1,1) = this%r**2
          ! ax_helm includes the mass matrix in h2
          this%coef%h2(i,1,1,1) = 1.0_rp
          ! mass matrix should be included here
          RHS%x(i,1,1,1) = X_in%x(i,1,1,1)*this%coef%B(i,1,1,1)
       end do
    end if
    this%coef%ifh2 = .true.

    ! This is a good idea from Niels' email!
    ! copy the unfiltered design as an initial guess for the filtered design
    ! to improved convergence
    call field_copy(X_out, X_in)

    ! gather scatter
    call this%coef%gs_h%op(RHS, GS_OP_ADD)
    ! set BCs
    call bc_list_apply_scalar(this%bclst_filt, RHS%x, n)

    ! Solve Helmholtz equation
    call profiler_start_region('filter solve')
    this%ksp_results(1) = &
         this%ksp_filt%solve(this%Ax, X_out, RHS%x, n, this%coef, &
         this%bclst_filt, this%coef%gs_h)

    call profiler_end_region

    ! update preconditioner (needed?)
    call this%pc_filt%update()

    ! write it all out
    call neko_log%message('Filter')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') this%ksp_results%iter, &
         this%ksp_results%res_start, this%ksp_results%res_final
    call neko_log%message(log_buf)

    call neko_scratch_registry%relinquish_field(temp_indices)



  end subroutine PDE_filter_apply

  !> Apply the adjoint filter
  !! @param X_in unfiltered field
  !! @param dF_dX_in is the sensitivity with respect to the unfiltered design
  !! @param dF_dX_out is the sensitivity with respect to the filtered design
  ! TODO
  ! this really confuses me!
  ! it's not really a chain rule back, it's just a filtering of the sensitivity
  ! ?
  !
  ! Update:
  ! After an email exchange with Niels:
  ! We DON'T want to be filtering the sensitivity, this IS just the chain rule.
  ! So to the best of my knowledge, we're just applying the same filter
  ! on the sensitivity field.
  !
  ! Niels did mention the order of the RHS assembly should be reversed however.
  ! I'm not exactly sure how this applies to us, but it should be brought up
  ! in the next group meeting.
  subroutine PDE_filter_apply_backward(this, dF_dX_in, dF_dX_out, X_in)
    class(PDE_filter_t), intent(inout) :: this
    type(field_t), intent(in) :: X_in
    type(field_t), intent(in) :: df_dX_out
    type(field_t), intent(inout) :: df_dX_in
    integer :: n, i
    type(field_t), pointer :: RHS
    integer :: temp_indices(1)
    type(field_t), pointer :: ta1, ta2, ta3
    character(len=LOG_SIZE) :: log_buf

    n = this%coef%dof%size()

    call neko_scratch_registry%request_field(RHS, temp_indices(1))

    ! set up Helmholtz operators and RHS
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! TODO
       ! I think this is correct but I've never tested it
       call device_cfill(this%coef%h1_d, this%r**2, n)
       call device_cfill(this%coef%h2_d, 1.0_rp, n)
       call device_col3(RHS%x_d, dF_dX_out%x_d, this%coef%B_d, n)
    else
       do i = 1, n
          ! h1 is already negative in its definition
          this%coef%h1(i,1,1,1) = this%r**2
          ! ax_helm includes the mass matrix in h2
          this%coef%h2(i,1,1,1) = 1.0_rp
          ! mass matrix should be included here
          RHS%x(i,1,1,1) = dF_dX_out%x(i,1,1,1) * this%coef%B(i,1,1,1)
       end do
    end if
    this%coef%ifh2 = .true.

    ! Same trick as before
    call field_copy(dF_dX_in, dF_dX_out)

    ! set BCs
    call bc_list_apply_scalar(this%bclst_filt, RHS%x, n)

    ! gather scatter
    call this%coef%gs_h%op(RHS, GS_OP_ADD)

    ! Solve Helmholtz equation
    call profiler_start_region('filter solve')
    this%ksp_results(1) = &
         this%ksp_filt%solve(this%Ax, dF_dX_in, RHS%x, n, this%coef, &
         this%bclst_filt, this%coef%gs_h)

    call profiler_end_region

    ! update preconditioner (needed?)
    call this%pc_filt%update()

    ! write it all out
    call neko_log%message('Filter')

    write(log_buf, '(A,A,A)') 'Iterations:   ',&
         'Start residual:     ', 'Final residual:'
    call neko_log%message(log_buf)
    write(log_buf, '(I11,3x, E15.7,5x, E15.7)') this%ksp_results%iter, &
         this%ksp_results%res_start, this%ksp_results%res_final
    call neko_log%message(log_buf)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine PDE_filter_apply_backward

  subroutine filter_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)

    implicit none
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(inout) :: coef
    type(dofmap_t), target, intent(inout) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type (pcp => pc)
      type is (jacobi_t)
       call pcp%init(coef, dof, gs)
      type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
      type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
      type is (hsmg_t)
       if (len_trim(pctype) .gt. 4) then
          if (index(pctype, '+') .eq. 5) then
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, &
                  bclst, trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    end select

    call ksp%set_pc(pc)

  end subroutine filter_precon_factory


end module PDE_filter
