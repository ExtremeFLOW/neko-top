!> @file adjoint_curl_curl.f90
!! @copyright
!! Copyright (c) 2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Adjoint of the curl-curl term in the Pn/Pn pressure residual.
!! @details The Karniadakis-Israeli-Orszag splitting of Neko's
!! \f$ P_N-P_N \f$ scheme puts the rotational form of the viscous term,
!! evaluated on the extrapolated velocity \f$ \mathbf{u}_e \f$, into the
!! pressure residual (`pnpn_prs_res`). Before gather-scatter, that
!! contribution is
!! \f[
!!   r_{cc}(\mathbf{u}_e) = -\frac{\mu}{\rho} G^T M \mathbf{W}, \qquad
!!   \mathbf{W} = M D_c M D_c \mathbf{u}_e,
!! \f]
!! where \f$ D_c \f$ is the element-local, pointwise curl,
!! \f$ M \mathbf{y} = \bar{B}^{-1} \mathrm{gs}(B \mathbf{y}) \f$ is the
!! mass-weighted average that Neko's `curl` applies to it (so that
!! `curl` is \f$ M D_c \f$), with \f$ B \f$ the element-local diagonal mass
!! matrix and \f$ \bar{B} = \mathrm{gs}(B) \f$ the assembled one (`coef%Binv`
!! holds \f$ \bar{B}^{-1} \f$), and \f$ G \f$ is the weak gradient `opgrad`,
!! \f$ G_i p = B D_i p \f$, whose transpose is `cdtp`. Since
!! \f$ \bar{B} = \mathrm{gs}(B) \f$, \f$ M \f$ maps every field to a
!! continuous one and leaves continuous fields unchanged, so it is a
!! projection, \f$ M^2 = M \f$, and
!! \f[
!!   r_{cc}(\mathbf{u}_e) = -\frac{\mu}{\rho} G^T M D_c M D_c \mathbf{u}_e .
!! \f]
!! The adjoint of this term with respect to \f$ \mathbf{u}_e \f$, acting on
!! the adjoint pressure \f$ p^\dagger \f$, therefore carries two factors
!! \f$ M^T \f$ and is the volume load
!! \f[
!!   \mathbf{L} = \frac{\mu}{\rho} D_c^T M^T D_c^T M^T G p^\dagger, \qquad
!!   M^T \mathbf{y} = B \, \mathrm{gs}(\bar{B}^{-1} \mathbf{y}),
!! \f]
!! i.e. for every continuous \f$ \mathbf{u}_e \f$ and \f$ p^\dagger \f$,
!! \f$ \mathbf{u}_e \cdot \mathbf{L} = -p^\dagger \cdot
!! r_{cc}(\mathbf{u}_e) \f$, with both dot products summed over all local
!! (non-unique) nodes and neither \f$ \mathbf{L} \f$ nor \f$ r_{cc} \f$
!! gather-scattered.
!!
!! The transpose is exact in three dimensions only. For a two-dimensional
!! mesh Neko's `curl` drops the \f$ z \f$ derivatives and does not
!! mass-average the \f$ z \f$ component, whereas \f$ \mathbf{L} \f$ applies
!! the full three-dimensional \f$ D_c^T \f$ and \f$ M^T \f$ to all three
!! components. The adjoint case therefore rejects two-dimensional meshes,
!! and `adjoint_curl_curl_check` rejects the boundary conditions for which
!! the load is not the transpose of the primal term.
!!
!! In the continuous setting \f$ \nabla \times \nabla p^\dagger = 0 \f$ and
!! the term reduces to a boundary integral. The discrete operators do not
!! satisfy that identity, so a boundary-only form is not the transpose of
!! what the primal computes; the load acts on the whole volume.
module adjoint_curl_curl
  use num_types, only: rp
  use field, only: field_t
  use coefs, only: coef_t
  use gather_scatter, only: gs_t
  use gs_ops, only: GS_OP_ADD
  use operators, only: opgrad, cdtp
  use field_math, only: field_sub2, field_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: col2
  use device_math, only: device_col2
  use facet_normal, only: facet_normal_t
  use comm, only: NEKO_COMM
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR
  use utils, only: neko_error
  implicit none
  private

  public :: adjoint_curl_curl_load, adjoint_curl_curl_check

contains

  !> Stop with an error if the boundary conditions make the adjoint
  !! curl-curl load differ from the transpose of the primal's curl-curl
  !! pressure term.
  !! @details The load does not transpose the symmetry-surface part of the
  !! primal term (`bc_sym_surface`), nor the cyclic rotations applied inside
  !! the primal's `curl`. The symmetry test is reduced over all ranks, so
  !! that every rank stops, also those that hold no symmetry facet. The
  !! restriction to three-dimensional meshes (see the module description) is
  !! not tested here: it has to be enforced before the adjoint fluid is
  !! initialised, which assumes a three-dimensional mesh.
  !! @param prim_sym The symmetry surface of the primal pressure residual.
  !! @param adj_sym The symmetry surface of the adjoint pressure residual.
  !! @param cyclic Whether the primal has cyclic periodic boundaries.
  subroutine adjoint_curl_curl_check(prim_sym, adj_sym, cyclic)
    type(facet_normal_t), intent(in) :: prim_sym, adj_sym
    logical, intent(in) :: cyclic
    logical :: has_sym

    has_sym = prim_sym%marked_facet%size() .gt. 0 .or. &
         adj_sym%marked_facet%size() .gt. 0
    call MPI_Allreduce(MPI_IN_PLACE, has_sym, 1, MPI_LOGICAL, MPI_LOR, &
         NEKO_COMM)

    if (has_sym) then
       call neko_error("Symmetry boundaries " // &
            "(case.fluid.boundary_conditions or " // &
            "case.adjoint_fluid.boundary_conditions) are not supported " // &
            "by the adjoint: the transpose of the curl-curl pressure " // &
            "term on symmetry surfaces is not implemented.")
    end if

    if (cyclic) then
       call neko_error("Cyclic boundaries (case.fluid.cyclic) are not " // &
            "supported by the adjoint: the transpose of the cyclic " // &
            "rotations in the curl-curl pressure term is not implemented.")
    end if

  end subroutine adjoint_curl_curl_check

  !> Compute the adjoint curl-curl load
  !! \f$ \mathbf{L} = \frac{\mu}{\rho} D_c^T M^T D_c^T M^T G p^\dagger \f$.
  !! @details See the module description for the notation. The load is
  !! element-local (not gather-scattered), since it is added to the adjoint
  !! forcing, which the velocity residual gather-scatters. It is added to the
  !! forcing before `makeabf`, which scales the forcing by \f$ \rho \f$ and
  !! extrapolates it in time, so the velocity right-hand side receives
  !! \f$ \mu D_c^T M^T D_c^T M^T G p^\dagger \f$, weighted over the lagged
  !! adjoint pressures as the primal's extrapolation weights
  !! \f$ \mathbf{u}_e \f$. The material properties are taken as constant, as
  !! in the primal residual. The cyclic rotations and the symmetry-surface
  !! part of the primal term are not transposed, so the load is not valid for
  !! cyclic or symmetry boundaries (see `adjoint_curl_curl_check`), nor for a
  !! two-dimensional mesh. All field arguments must be distinct.
  !! @param load_x x component of the load.
  !! @param load_y y component of the load.
  !! @param load_z z component of the load.
  !! @param p_adj The adjoint pressure (continuous).
  !! @param coef The SEM coefficients.
  !! @param gs The gather-scatter operator of the velocity space.
  !! @param mu The constant dynamic viscosity.
  !! @param rho The constant density.
  !! @param work1 Work field.
  !! @param work2 Work field.
  !! @param work3 Work field.
  !! @param work4 Work field.
  subroutine adjoint_curl_curl_load(load_x, load_y, load_z, p_adj, coef, gs, &
       mu, rho, work1, work2, work3, work4)
    type(field_t), intent(inout) :: load_x, load_y, load_z
    type(field_t), intent(in) :: p_adj
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: mu, rho
    type(field_t), intent(inout) :: work1, work2, work3, work4

    ! g = Binv gs(G p), so that B g = M^T G p; held in load
    call opgrad(load_x%x, load_y%x, load_z%x, p_adj%x, coef)
    call adjoint_curl_curl_gs_binv(load_x, load_y, load_z, coef, gs)

    ! c = Binv gs(D_c^T B g), so that B c = M^T D_c^T M^T G p; held in
    ! work1..3
    call adjoint_curl_transpose(work1, work2, work3, load_x, load_y, load_z, &
         work4, coef)
    call adjoint_curl_curl_gs_binv(work1, work2, work3, coef, gs)

    ! L = (mu / rho) D_c^T B c
    call adjoint_curl_transpose(load_x, load_y, load_z, work1, work2, work3, &
         work4, coef)
    call field_cmult(load_x, mu / rho)
    call field_cmult(load_y, mu / rho)
    call field_cmult(load_z, mu / rho)

  end subroutine adjoint_curl_curl_load

  !> Apply the transpose of the pointwise curl to mass-weighted data,
  !! \f$ \mathbf{w} = D_c^T B \mathbf{y} \f$, element by element (no
  !! gather-scatter): \f$ w_x = (B D_z)^T y_y - (B D_y)^T y_z \f$ and
  !! cyclically.
  !! @param w1 x component of the result.
  !! @param w2 y component of the result.
  !! @param w3 z component of the result.
  !! @param y1 x component of the input.
  !! @param y2 y component of the input.
  !! @param y3 z component of the input.
  !! @param work Work field.
  !! @param coef The SEM coefficients.
  subroutine adjoint_curl_transpose(w1, w2, w3, y1, y2, y3, work, coef)
    type(field_t), intent(inout) :: w1, w2, w3
    type(field_t), intent(inout) :: y1, y2, y3
    type(field_t), intent(inout) :: work
    type(coef_t), intent(in) :: coef

    call cdtp(w1%x, y2%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call cdtp(work%x, y3%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call field_sub2(w1, work)

    call cdtp(w2%x, y3%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(work%x, y1%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call field_sub2(w2, work)

    call cdtp(w3%x, y1%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(work%x, y2%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call field_sub2(w3, work)

  end subroutine adjoint_curl_transpose

  !> Apply \f$ \mathbf{y} \leftarrow \bar{B}^{-1} \mathrm{gs}(\mathbf{y}) \f$
  !! (`coef%Binv` times the gather-scatter sum) to the three components of a
  !! vector field.
  !! @param y1 x component.
  !! @param y2 y component.
  !! @param y3 z component.
  !! @param coef The SEM coefficients.
  !! @param gs The gather-scatter operator.
  subroutine adjoint_curl_curl_gs_binv(y1, y2, y3, coef, gs)
    type(field_t), intent(inout) :: y1, y2, y3
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n

    n = y1%size()

    call gs%op(y1, GS_OP_ADD)
    call gs%op(y2, GS_OP_ADD)
    call gs%op(y3, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(y1%x_d, coef%Binv_d, n)
       call device_col2(y2%x_d, coef%Binv_d, n)
       call device_col2(y3%x_d, coef%Binv_d, n)
    else
       call col2(y1%x, coef%Binv, n)
       call col2(y2%x, coef%Binv, n)
       call col2(y3%x, coef%Binv, n)
    end if

  end subroutine adjoint_curl_curl_gs_binv

end module adjoint_curl_curl
