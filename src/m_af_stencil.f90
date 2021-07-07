#include "cpp_macros.h"
!> This module contains functionality for dealing with numerical stencils
module m_af_stencil
  use m_af_types

  implicit none
  private

    !> Number of predefined stencil shapes
  integer, parameter, private :: num_shapes = 1

  !> 3/5/7 point stencil in 1D/2D/3D
  integer, parameter, public :: af_stencil_357 = 1

  !> Number of coefficients in the stencils
  integer, parameter, public :: af_stencil_sizes(num_shapes) = [2*NDIM+1]

  abstract interface
     !> Subroutine for setting a stencil on a box
     subroutine af_subr_stencil(box, stencil)
       import
       type(box_t), intent(in)        :: box     !< Box to set stencil for
       type(stencil_t), intent(inout) :: stencil !< Stencil to set
     end subroutine af_subr_stencil
  end interface

  public :: af_stencil_index
  public :: af_stencil_prepare_store
  public :: af_stencil_store_box
  public :: af_stencil_store
  public :: af_stencil_apply
  public :: af_stencil_apply_box
  public :: af_stencil_gsrb_box
  public :: af_stencil_get_box
  public :: af_stencil_try_constant

contains

  !> Get index of a stencil, or af_stencil_none is not present
  pure integer function af_stencil_index(box, key)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: key
    integer                 :: i

    af_stencil_index = af_stencil_none
    do i = 1, box%n_stencils
       if (box%stencils(i)%key == key) then
          af_stencil_index = i
          exit
       end if
    end do
  end function af_stencil_index

  !> Store a constant stencil
  subroutine af_stencil_store(tree, key, set_stencil)
    type(af_t), intent(inout)  :: tree
    integer, intent(in)        :: key !< Stencil key
    !> Method to set a stencil
    procedure(af_subr_stencil) :: set_stencil
    integer                    :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call af_stencil_store_box(tree%boxes(id), key, set_stencil)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_stencil_store

  !> Prepare to store a stencil
  subroutine af_stencil_prepare_store(box, key, ix)
    type(box_t), intent(inout)   :: box
    integer, intent(in)          :: key !< Stencil key
    integer, intent(out)         :: ix  !< Index to store stencil
    type(stencil_t), allocatable :: tmp(:)

    ix = af_stencil_index(box, key)

    if (ix == af_stencil_none) then
       ! Allocate storage if necessary
       if (box%n_stencils == 0) then
          allocate(box%stencils(1))
       else if (box%n_stencils == size(box%stencils)) then
          allocate(tmp(2 * box%n_stencils))
          tmp(1:box%n_stencils) = box%stencils
          call move_alloc(tmp, box%stencils)
       end if

       box%n_stencils       = box%n_stencils + 1
       ix                   = box%n_stencils
       box%stencils(ix)%key = key
    end if
  end subroutine af_stencil_prepare_store

  !> Store a stencil for a box
  subroutine af_stencil_store_box(box, key, set_stencil)
    type(box_t), intent(inout)   :: box
    integer, intent(in)          :: key !< Stencil key
    !> Method to set a stencil
    procedure(af_subr_stencil)   :: set_stencil
    integer                      :: ix

    call af_stencil_prepare_store(box, key, ix)
    call set_stencil(box, box%stencils(ix))
  end subroutine af_stencil_store_box

  subroutine af_stencil_apply(tree, key, iv, i_out)
    type(af_t), intent(inout) :: tree  !< Operate on this box
    integer, intent(in)       :: key   !< Stencil key
    integer, intent(in)       :: iv    !< Input variable
    integer, intent(in)       :: i_out !< Output variable

    integer :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call af_stencil_apply_box(tree%boxes(id), key, iv, i_out)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_stencil_apply

  subroutine af_stencil_apply_box(box, key, iv, i_out)
    type(box_t), intent(inout) :: box   !< Operate on this box
    integer, intent(in)        :: key   !< Stencil key
    integer, intent(in)        :: iv    !< Input variable
    integer, intent(in)        :: i_out !< Output variable
    integer                    :: i

    i = af_stencil_index(box, key)
    if (i == af_stencil_none) error stop "Stencil not stored"

    select case (box%stencils(i)%shape)
    case (af_stencil_357)
       call stencil_apply_357(box, box%stencils(i), iv, i_out)
    case default
       error stop "Unknown stencil"
    end select
  end subroutine af_stencil_apply_box

  !> Apply a 3/5/7-point stencil in 1D/2D/3D
  subroutine stencil_apply_357(box, stencil, iv, i_out)
    type(box_t), intent(inout)  :: box     !< Operate on this box
    type(stencil_t), intent(in) :: stencil !< Stencil
    integer, intent(in)         :: iv      !< Input variable
    integer, intent(in)         :: i_out   !< Output variable
    real(dp)                    :: c(2*NDIM+1)
    real(dp)                    :: tmp(DTIMES(box%n_cell))
    integer                     :: IJK
#if NDIM == 2
    real(dp)                    :: rfac(2, box%n_cell), c_cyl(2*NDIM+1)
#endif

    associate (cc => box%cc, nc => box%n_cell)
      if (stencil%constant) c = stencil%c

#if NDIM == 1
      do KJI_DO(1, nc)
         if (.not. stencil%constant) c = stencil%v(:, IJK)
         tmp(i) = sum(c * &
              [cc(i, iv), cc(i-1, iv), cc(i+1, iv)])
      end do; CLOSE_DO
#elif NDIM == 2
      if (stencil%cylindrical_gradient) then
         ! Correct for cylindrical coordinates, assuming the elements correspond
         ! to a gradient operation
         call af_cyl_flux_factors(box, rfac)
         do KJI_DO(1, nc)
            if (.not. stencil%constant) c = stencil%v(:, IJK)
            c_cyl(2:3) = rfac(1:2, i) * c(2:3)
            c_cyl(1) = c(1) - (c_cyl(2) - c(2)) - (c_cyl(3) - c(3))
            c_cyl(4:) = c(4:)
            tmp(i, j) = sum(c_cyl * &
                 [cc(i, j, iv), cc(i-1, j, iv), cc(i+1, j, iv), &
                 cc(i, j-1, iv), cc(i, j+1, iv)])
         end do; CLOSE_DO
      else
         do KJI_DO(1, nc)
            if (.not. stencil%constant) c = stencil%v(:, IJK)
            tmp(i, j) = sum(c * &
                 [cc(i, j, iv), cc(i-1, j, iv), cc(i+1, j, iv), &
                 cc(i, j-1, iv), cc(i, j+1, iv)])
         end do; CLOSE_DO
      end if
#elif NDIM == 3
      do KJI_DO(1, nc)
         if (.not. stencil%constant) c = stencil%v(:, IJK)
         tmp(i, j, k) = sum(c * &
              [cc(i, j, k, iv), cc(i-1, j, k, iv), cc(i+1, j, k, iv), &
              cc(i, j-1, k, iv), cc(i, j+1, k, iv), &
              cc(i, j, k-1, iv), cc(i, j, k+1, iv)])
      end do; CLOSE_DO
#endif

      ! if (box%n_bc > 0) then
      !    call stencil_correct_bc_357(box, stencil, iv, tmp)
      ! end if
      if (allocated(stencil%bc_correction)) then
         tmp = tmp - stencil%bc_correction
      end if

      cc(DTIMES(1:nc), i_out) = tmp
    end associate

  end subroutine stencil_apply_357

  subroutine stencil_correct_bc_357(box, stencil, iv, tmp)
    type(box_t), intent(inout)  :: box     !< Operate on this box
    type(stencil_t), intent(in) :: stencil !< Stencil
    integer, intent(in)         :: iv      !< Input variable
    real(dp), intent(inout)     :: tmp(DTIMES(box%n_cell))
    integer                     :: bc_ix, bc_type, nb, lo(NDIM), hi(NDIM)
    integer                     :: nb_dim, olo(NDIM), ohi(NDIM)
    real(dp)                    :: bc_val(box%n_cell**(NDIM-1))
    real(dp)                    :: stencil_coeffs(box%n_cell**(NDIM-1))
    real(dp)                    :: values_inside(box%n_cell**(NDIM-1))
    real(dp)                    :: values_outside(box%n_cell**(NDIM-1))

    do bc_ix = 1, box%n_bc
       nb = box%bc_index_to_nb(bc_ix)

       bc_val = box%bc_val(:, iv, bc_ix)
       bc_type = box%bc_type(iv, bc_ix)

       ! Determine index range next to boundary
       call af_get_index_bc_inside(nb, box%n_cell, 1, lo, hi)
       call af_get_index_bc_outside(nb, box%n_cell, 1, olo, ohi)

       ! Get stencil coefficients near boundary
       if (stencil%constant) then
          stencil_coeffs = stencil%c(nb+1)
       else
          stencil_coeffs = pack(stencil%v(nb+1, &
               DSLICE(lo, hi)), .true.)
       end if

       values_inside = pack(box%cc(DSLICE(lo, hi), iv), .true.)
       values_outside = pack(box%cc(DSLICE(olo, ohi), iv), .true.)

       select case (bc_type)
       case (af_bc_dirichlet)
          ! Dirichlet value at cell face, so compute gradient over h/2
          ! E.g. 1 -2 1 becomes 0 -3 1 for a 1D Laplacian
          tmp(DSLICE(lo, hi)) = tmp(DSLICE(lo, hi)) &
               + reshape(2 * stencil_coeffs * bc_val - &
               stencil_coeffs * (values_inside + values_outside), &
               hi - lo + 1)
       case (af_bc_neumann)
          nb_dim = af_neighb_dim(nb)
          ! E.g. 1 -2 1 becomes 0 -1 1 for a 1D Laplacian
          tmp(DSLICE(lo, hi)) = tmp(DSLICE(lo, hi)) &
               - reshape(stencil_coeffs * (values_outside - values_inside) &
               - stencil_coeffs * af_neighb_high_pm(nb) * box%dr(nb_dim) * &
               bc_val, hi - lo + 1)
       case default
          error stop "unsupported boundary condition"
       end select
    end do
  end subroutine stencil_correct_bc_357

  !> Perform Gauss-Seidel red-black on a stencil
  subroutine af_stencil_gsrb_box(box, key, redblack, iv, i_rhs)
    type(box_t), intent(inout) :: box      !< Operate on this box
    integer, intent(in)        :: key      !< Stencil key
    integer, intent(in)        :: redblack !< Even or odd integer
    integer, intent(in)        :: iv       !< Solve for variable
    integer, intent(in)        :: i_rhs    !< Right-hand side
    integer                    :: i

    i = af_stencil_index(box, key)
    if (i == af_stencil_none) error stop "Stencil not stored"

    select case (box%stencils(i)%shape)
    case (af_stencil_357)
       call stencil_gsrb_357(box, box%stencils(i), redblack, iv, i_rhs)
    case default
       error stop "Unknown stencil"
    end select
  end subroutine af_stencil_gsrb_box

  !> Apply a constant 3/5/7-point stencil in 1D/2D/3D
  subroutine stencil_gsrb_357(box, stencil, redblack, iv, i_rhs)
    type(box_t), intent(inout)  :: box      !< Operate on this box
    type(stencil_t), intent(in) :: stencil  !< Stencil
    integer, intent(in)         :: redblack !< Even or odd integer
    integer, intent(in)         :: iv       !< Solve for variable
    integer, intent(in)         :: i_rhs    !< Right-hand side

    real(dp) :: rhs(DTIMES(box%n_cell))
    real(dp) :: c(2*NDIM+1)
    real(dp) :: tmp(DTIMES(box%n_cell))
    integer  :: IJK, i0, nc
#if NDIM == 2
    real(dp) :: rfac(2, box%n_cell), c_cyl(2*NDIM+1)
#endif

    nc = box%n_cell

    if (stencil%constant) then
       c = stencil%c
    else
       c = 0                    ! avoid warning
    end if

    associate (cc => box%cc, nc => box%n_cell)
      tmp = cc(DTIMES(1:nc), iv)
      rhs = cc(DTIMES(1:nc), i_rhs)

      if (allocated(stencil%bc_correction)) then
         rhs = rhs + stencil%bc_correction
      end if

#if NDIM == 1
      i0 = 2 - iand(redblack, 1)
      do i = i0, nc, 2
         if (.not. stencil%constant) c = stencil%v(:, IJK)
         tmp(IJK) = (rhs(IJK) - &
              sum(c(2:) * [cc(i-1, iv), cc(i+1, iv)])) / c(1)
      end do
#elif NDIM == 2
      if (stencil%cylindrical_gradient) then
         ! Correct for cylindrical coordinates, assuming the elements correspond
         ! to a gradient operation
         call af_cyl_flux_factors(box, rfac)
         do j = 1, nc
            i0 = 2 - iand(ieor(redblack, j), 1)
            do i = i0, nc, 2
               if (.not. stencil%constant) c = stencil%v(:, IJK)
               c_cyl(2:3) = rfac(1:2, i) * c(2:3)
               ! Correct central coefficient
               c_cyl(1) = c(1) - (c_cyl(2) - c(2)) - (c_cyl(3) - c(3))
               c_cyl(4:) = c(4:)

               tmp(IJK) = (rhs(IJK) - &
                    sum(c_cyl(2:) * [cc(i-1, j, iv), cc(i+1, j, iv), &
                    cc(i, j-1, iv), cc(i, j+1, iv)])) / c_cyl(1)
            end do
         end do
      else
         do j = 1, nc
            i0 = 2 - iand(ieor(redblack, j), 1)
            do i = i0, nc, 2
               if (.not. stencil%constant) c = stencil%v(:, IJK)
               tmp(IJK) = (rhs(IJK) - &
                    sum(c(2:) * [cc(i-1, j, iv), cc(i+1, j, iv), &
                    cc(i, j-1, iv), cc(i, j+1, iv)])) / c(1)
            end do
         end do
      end if
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            i0 = 2 - iand(ieor(redblack, k+j), 1)
            do i = i0, nc, 2
               if (.not. stencil%constant) c = stencil%v(:, IJK)
               tmp(IJK) = (rhs(IJK) - sum(c(2:) * &
                    [cc(i-1, j, k, iv), cc(i+1, j, k, iv), &
                    cc(i, j-1, k, iv), cc(i, j+1, k, iv), &
                    cc(i, j, k-1, iv), cc(i, j, k+1, iv)])) / c(1)
            end do
         end do
      end do
#endif

      cc(DTIMES(1:nc), iv) = tmp
    end associate
  end subroutine stencil_gsrb_357

  !> Convert a variable stencil to constant one if possible
  subroutine af_stencil_try_constant(box, ix, abs_tol)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: ix  !< Stencil index
    real(dp), intent(in)       :: abs_tol !< Absolute tolerance
    integer                    :: IJK, nc, n_coeff

    if (.not. allocated(box%stencils(ix)%v)) &
         error stop "No variable stencil present"

    nc = box%n_cell

    ! Check if all coefficients are the same, otherwise return
    do KJI_DO(1, nc)
       if (any(abs(box%stencils(ix)%v(:, IJK) - &
            box%stencils(ix)%v(:, DTIMES(1))) > abs_tol)) return
    end do; CLOSE_DO

    ! print *, "CONSTANT"
    box%stencils(ix)%constant = .true.
    n_coeff = size(box%stencils(ix)%v, 1)
    allocate(box%stencils(ix)%c(n_coeff))
    box%stencils(ix)%c = box%stencils(ix)%v(:, DTIMES(1))
    deallocate(box%stencils(ix)%v)
  end subroutine af_stencil_try_constant

  !> Get the stencil for a box
  subroutine af_stencil_get_box(box, key, v)
    type(box_t), intent(inout) :: box !< Operate on this box
    integer, intent(in)        :: key !< Stencil key
    !> Stencil coefficients
    real(dp), intent(inout)    :: v(:, DTIMES(:))
    integer                    :: i

    i = af_stencil_index(box, key)
    if (i == af_stencil_none) error stop "Stencil not stored"

    select case (box%stencils(i)%shape)
    case (af_stencil_357)
       call stencil_get_357(box, box%stencils(i), v)
    case default
       error stop "Unknown stencil"
    end select
  end subroutine af_stencil_get_box

  !> Get the stencil for all cells in a box
  subroutine stencil_get_357(box, stencil, v)
    type(box_t), intent(inout)  :: box     !< Box
    type(stencil_t), intent(in) :: stencil !< Stencil
    !> Stencil coefficients
    real(dp), intent(inout)     :: v(:, DTIMES(:))
#if NDIM == 2
    real(dp)                    :: c_cyl(size(v, 1)), rfac(2, box%n_cell)
#endif
    integer                     :: IJK, nc

    nc = box%n_cell

    if (size(v, 1) /= af_stencil_sizes(stencil%shape)) &
         error stop "Argument v has wrong size"

    if (stencil%constant) then
       do KJI_DO(1, nc)
          v(:, IJK) = stencil%c
       end do; CLOSE_DO
    else
       v = stencil%v
    end if

#if NDIM == 2
    if (stencil%cylindrical_gradient) then
       ! Correct for cylindrical coordinates, assuming the elements correspond
       ! to a gradient operation
       call af_cyl_flux_factors(box, rfac)
       do KJI_DO(1, nc)
          c_cyl(2:3) = rfac(1:2, i) * v(2:3, IJK)
          c_cyl(1) = v(1, IJK) - (c_cyl(2) - v(2, IJK)) &
               - (c_cyl(3) - v(3, IJK))
          c_cyl(4:) = v(4:, IJK)
          v(:, IJK) = c_cyl
       end do; CLOSE_DO
    end if
#endif
  end subroutine stencil_get_357

end module m_af_stencil
