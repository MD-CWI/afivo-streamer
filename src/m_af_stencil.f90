#include "cpp_macros.h"
!> This module contains functionality for dealing with numerical stencils
module m_af_stencil
  use m_af_types

  implicit none
  private

  integer, parameter, public :: stencil_constant = 1 !< Constant stencil
  integer, parameter, public :: stencil_variable = 2 !< Variable stencil
  integer, parameter, public :: stencil_sparse   = 3 !< Sparse stencil

  !> Number of predefined stencil shapes
  integer, parameter, private :: num_shapes = 5

  !> 3/5/7 point stencil in 1D/2D/3D
  integer, parameter, public :: af_stencil_357 = 1

  !> Prolongation stencil using nearest 2, 3, 4 neighbors in 1D-3D
  integer, parameter, public :: af_stencil_p234 = 2

  !> Prolongation stencil using nearest 2, 4, 8 neighbors in 1D-3D
  integer, parameter, public :: af_stencil_p248 = 3

  !> Stencil for direct neighbors
  integer, parameter, public :: af_stencil_246 = 4

  !> Stencil for masking
  integer, parameter, public :: af_stencil_mask = 5

  !> Number of coefficients in the stencils
  integer, parameter, public :: af_stencil_sizes(num_shapes) = &
       [2*NDIM+1, NDIM+1, 2**NDIM, 2*NDIM, 1]

  abstract interface
     !> Subroutine for setting a stencil on a box
     subroutine af_subr_stencil(box, stencil)
       import
       type(box_t), intent(in)        :: box     !< Box to set stencil for
       type(stencil_t), intent(inout) :: stencil !< Stencil to set
     end subroutine af_subr_stencil
  end interface

  public :: af_stencil_print_info
  public :: af_stencil_index
  public :: af_stencil_prepare_store
  public :: af_stencil_store_box
  public :: af_stencil_check_box
  public :: af_stencil_store
  public :: af_stencil_apply
  public :: af_stencil_apply_box
  public :: af_stencil_gsrb_box
  public :: af_stencil_prolong_box
  public :: af_stencil_get_box
  public :: af_stencil_try_constant

contains

  !> Print statistics about the stencils in the tree
  subroutine af_stencil_print_info(tree)
    type(af_t), intent(in) :: tree
    integer                :: id, n_stencils_stored
    integer                :: n_boxes_with_stencils
    integer                :: n_constant_stored, n_variable_stored
    integer                :: n, n_stencils, n_sparse_stored

    n_stencils_stored     = 0
    n_constant_stored     = 0
    n_variable_stored     = 0
    n_sparse_stored       = 0
    n_boxes_with_stencils = 0

    do id = 1, tree%highest_id
       if (.not. tree%boxes(id)%in_use) cycle

       n_stencils = tree%boxes(id)%n_stencils
       if (n_stencils > 0) then
          n_boxes_with_stencils = n_boxes_with_stencils + 1
          n_stencils_stored = n_stencils_stored + n_stencils

          do n = 1, n_stencils
             select case (tree%boxes(id)%stencils(n)%stype)
             case (stencil_constant)
                n_constant_stored = n_constant_stored + 1
             case (stencil_variable)
                n_variable_stored = n_variable_stored + 1
             case (stencil_sparse)
                n_sparse_stored = n_sparse_stored + 1
             end select
          end do
       end if
    end do

    write(*, '(A)')     ''
    write(*, '(A)')     ' ** Information about stencils **'
    write(*, '(A, I0)') ' #boxes with stencils: ', n_boxes_with_stencils
    write(*, '(A, I0)') ' #stencils stored:     ', n_stencils_stored
    write(*, '(A, I0)') ' #constant stencils:   ', n_constant_stored
    write(*, '(A, I0)') ' #variable stencils:   ', n_variable_stored
    write(*, '(A, I0)') ' #sparse stencils:     ', n_sparse_stored
  end subroutine af_stencil_print_info

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
    else
       error stop "Stencil with this key has already been stored"
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

  !> Check if a stencil was correctly stored for a box
  subroutine af_stencil_check_box(box, key, ix)
    type(box_t), intent(in) :: box !< Box
    integer, intent(in)     :: key !< Stencil key
    integer, intent(in)     :: ix  !< Stencil index

    if (key == af_stencil_none) &
         error stop "key == af_stencil_none"

    if (af_stencil_index(box, key) /= ix) &
         error stop "Stencil with key not found at index"

    associate (stencil => box%stencils(ix))
      if (stencil%shape < 1 .or. stencil%shape > num_shapes) &
           error stop "Unknown stencil shape"
      select case (stencil%stype)
      case (stencil_constant)
         if (.not. allocated(stencil%c)) &
              error stop "stencil%c not allocated"
      case (stencil_variable)
         if (.not. allocated(stencil%v)) &
              error stop "stencil%v not allocated"
      case (stencil_sparse)
         if (.not. allocated(stencil%sparse_ix)) &
              error stop "stencil%sparse_ix not allocated"
         if (.not. allocated(stencil%sparse_v)) &
              error stop "stencil%sparse_v not allocated"
      case default
         error stop "Unknow stencil%stype"
      end select
    end associate
  end subroutine af_stencil_check_box

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
    integer                     :: IJK
#if NDIM == 2
    real(dp) :: rfac(2, box%n_cell), c_cyl(2*NDIM+1)
    real(dp) :: cc_cyl(2*NDIM+1, box%n_cell)
#endif

    if (iv == i_out) error stop "Cannot have iv == i_out"
    if (stencil%stype == stencil_sparse) error stop "sparse not implemented"

    associate (cc => box%cc, nc => box%n_cell)

#if NDIM == 1
      if (stencil%stype == stencil_constant) then
         c = stencil%c
         do KJI_DO(1, nc)
            cc(i, i_out) = &
                 c(1) * cc(i, iv) + c(2) * cc(i-1, iv) + c(3) * cc(i+1, iv)
         end do; CLOSE_DO
      else
         do KJI_DO(1, nc)
            c = stencil%v(:, IJK)
            cc(i, i_out) = &
                 c(1) * cc(i, iv) + c(2) * cc(i-1, iv) + c(3) * cc(i+1, iv)
         end do; CLOSE_DO
      end if
#elif NDIM == 2
      if (stencil%cylindrical_gradient) then
         ! Correct for cylindrical coordinates, assuming the elements correspond
         ! to a gradient operation
         call af_cyl_flux_factors(box, rfac)

         if (stencil%stype == stencil_constant) then
            c = stencil%c

            ! Pre-compute coefficients for each i-index
            do i = 1, nc
               cc_cyl(2:3, i) = rfac(1:2, i) * c(2:3)
               cc_cyl(1, i)   = c(1) - (cc_cyl(2, i) - c(2)) &
                    - (cc_cyl(3, i) - c(3))
               cc_cyl(4:, i)  = c(4:)
            end do

            do KJI_DO(1, nc)
               cc(i, j, i_out) = &
                    cc_cyl(1, i) * cc(i, j, iv) + &
                    cc_cyl(2, i) * cc(i-1, j, iv) + &
                    cc_cyl(3, i) * cc(i+1, j, iv) + &
                    cc_cyl(4, i) * cc(i, j-1, iv) + &
                    cc_cyl(5, i) * cc(i, j+1, iv)
            end do; CLOSE_DO
         else
            ! Variable stencil
            do KJI_DO(1, nc)
               c = stencil%v(:, IJK)
               c_cyl(2:3) = rfac(1:2, i) * c(2:3)
               c_cyl(1) = c(1) - (c_cyl(2) - c(2)) - (c_cyl(3) - c(3))
               c_cyl(4:) = c(4:)
               cc(i, j, i_out) = &
                    c_cyl(1) * cc(i, j, iv) + &
                    c_cyl(2) * cc(i-1, j, iv) + &
                    c_cyl(3) * cc(i+1, j, iv) + &
                    c_cyl(4) * cc(i, j-1, iv) + &
                    c_cyl(5) * cc(i, j+1, iv)
            end do; CLOSE_DO
         end if
      else
         ! No cylindrical gradient correction
         if (stencil%stype == stencil_constant) then
            c = stencil%c
            do KJI_DO(1, nc)
               cc(i, j, i_out) = &
                    c(1) * cc(i, j, iv) + &
                    c(2) * cc(i-1, j, iv) + &
                    c(3) * cc(i+1, j, iv) + &
                    c(4) * cc(i, j-1, iv) + &
                    c(5) * cc(i, j+1, iv)
            end do; CLOSE_DO
         else
            do KJI_DO(1, nc)
               c = stencil%v(:, IJK)
               cc(i, j, i_out) = &
                    c(1) * cc(i, j, iv) + &
                    c(2) * cc(i-1, j, iv) + &
                    c(3) * cc(i+1, j, iv) + &
                    c(4) * cc(i, j-1, iv) + &
                    c(5) * cc(i, j+1, iv)
            end do; CLOSE_DO
         end if
      end if
#elif NDIM == 3
      if (stencil%stype == stencil_constant) then
         c = stencil%c
         do KJI_DO(1, nc)
            cc(i, j, k, i_out) = &
                 c(1) * cc(i, j, k, iv) + &
                 c(2) * cc(i-1, j, k, iv) + &
                 c(3) * cc(i+1, j, k, iv) + &
                 c(4) * cc(i, j-1, k, iv) + &
                 c(5) * cc(i, j+1, k, iv) + &
                 c(6) * cc(i, j, k-1, iv) + &
                 c(7) * cc(i, j, k+1, iv)
         end do; CLOSE_DO
      else
         do KJI_DO(1, nc)
            c = stencil%v(:, IJK)
            cc(i, j, k, i_out) = &
                 c(1) * cc(i, j, k, iv) + &
                 c(2) * cc(i-1, j, k, iv) + &
                 c(3) * cc(i+1, j, k, iv) + &
                 c(4) * cc(i, j-1, k, iv) + &
                 c(5) * cc(i, j+1, k, iv) + &
                 c(6) * cc(i, j, k-1, iv) + &
                 c(7) * cc(i, j, k+1, iv)
         end do; CLOSE_DO
      end if
#endif

      if (allocated(stencil%bc_correction)) then
         cc(DTIMES(1:nc), i_out) = cc(DTIMES(1:nc), i_out) - &
              stencil%bc_correction
      end if
    end associate

  end subroutine stencil_apply_357

  ! subroutine stencil_correct_bc_357(box, stencil, iv, tmp)
  !   type(box_t), intent(inout)  :: box     !< Operate on this box
  !   type(stencil_t), intent(in) :: stencil !< Stencil
  !   integer, intent(in)         :: iv      !< Input variable
  !   real(dp), intent(inout)     :: tmp(DTIMES(box%n_cell))
  !   integer                     :: bc_ix, bc_type, nb, lo(NDIM), hi(NDIM)
  !   integer                     :: nb_dim, olo(NDIM), ohi(NDIM)
  !   real(dp)                    :: bc_val(box%n_cell**(NDIM-1))
  !   real(dp)                    :: stencil_coeffs(box%n_cell**(NDIM-1))
  !   real(dp)                    :: values_inside(box%n_cell**(NDIM-1))
  !   real(dp)                    :: values_outside(box%n_cell**(NDIM-1))

  !   do bc_ix = 1, box%n_bc
  !      nb = box%bc_index_to_nb(bc_ix)

  !      bc_val = box%bc_val(:, iv, bc_ix)
  !      bc_type = box%bc_type(iv, bc_ix)

  !      ! Determine index range next to boundary
  !      call af_get_index_bc_inside(nb, box%n_cell, 1, lo, hi)
  !      call af_get_index_bc_outside(nb, box%n_cell, 1, olo, ohi)

  !      ! Get stencil coefficients near boundary
  !      if (stencil%constant) then
  !         stencil_coeffs = stencil%c(nb+1)
  !      else
  !         stencil_coeffs = pack(stencil%v(nb+1, &
  !              DSLICE(lo, hi)), .true.)
  !      end if

  !      values_inside = pack(box%cc(DSLICE(lo, hi), iv), .true.)
  !      values_outside = pack(box%cc(DSLICE(olo, ohi), iv), .true.)

  !      select case (bc_type)
  !      case (af_bc_dirichlet)
  !         ! Dirichlet value at cell face, so compute gradient over h/2
  !         ! E.g. 1 -2 1 becomes 0 -3 1 for a 1D Laplacian
  !         tmp(DSLICE(lo, hi)) = tmp(DSLICE(lo, hi)) &
  !              + reshape(2 * stencil_coeffs * bc_val - &
  !              stencil_coeffs * (values_inside + values_outside), &
  !              hi - lo + 1)
  !      case (af_bc_neumann)
  !         nb_dim = af_neighb_dim(nb)
  !         ! E.g. 1 -2 1 becomes 0 -1 1 for a 1D Laplacian
  !         tmp(DSLICE(lo, hi)) = tmp(DSLICE(lo, hi)) &
  !              - reshape(stencil_coeffs * (values_outside - values_inside) &
  !              - stencil_coeffs * af_neighb_high_pm(nb) * box%dr(nb_dim) * &
  !              bc_val, hi - lo + 1)
  !      case default
  !         error stop "unsupported boundary condition"
  !      end select
  !   end do
  ! end subroutine stencil_correct_bc_357

  subroutine af_stencil_prolong_box(box_p, box_c, key, iv, iv_to, add)
    type(box_t), intent(in)       :: box_p    !< Parent box
    type(box_t), intent(inout)    :: box_c !< Child box
    integer, intent(in)           :: key   !< Stencil key
    integer, intent(in)           :: iv    !< Input variable
    integer, intent(in)           :: iv_to !< Destination variable
    logical, intent(in), optional :: add   !< Add to old values
    logical                       :: add_to
    integer                       :: i

    add_to = .false.
    if (present(add)) add_to = add

    i = af_stencil_index(box_c, key)
    if (i == af_stencil_none) error stop "Stencil not stored"

    select case (box_c%stencils(i)%shape)
    case (af_stencil_p234)
       call stencil_prolong_234(box_p, box_c, box_c%stencils(i), &
            iv, iv_to, add_to)
    case (af_stencil_p248)
       call stencil_prolong_248(box_p, box_c, box_c%stencils(i), &
            iv, iv_to, add_to)
    case default
       print *, "box_c%stencils(i)%shape: ", box_c%stencils(i)%shape
       error stop "Unknown stencil"
    end select
  end subroutine af_stencil_prolong_box

  !> Linear prolongation using nearest NDIM+1 neighbors
  subroutine stencil_prolong_234(box_p, box_c, stencil, iv, iv_to, add)
    type(box_t), intent(in)     :: box_p      !< Parent box
    type(box_t), intent(inout)  :: box_c   !< Child box
    type(stencil_t), intent(in) :: stencil !< Stencil
    integer, intent(in)         :: iv      !< Input variable
    integer, intent(in)         :: iv_to   !< Destination variable
    logical, intent(in)         :: add     !< Add to old values

    integer  :: nc, ix_offset(NDIM), IJK
    integer  :: IJK_(c1), IJK_(c2)
    real(dp) :: c(NDIM+1)

    nc        = box_c%n_cell
    ix_offset = af_get_child_offset(box_c)
    if (stencil%stype == stencil_sparse) error stop "sparse not implemented"
    if (.not. add) box_c%cc(DTIMES(1:nc), iv_to) = 0

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
    if (stencil%stype == stencil_constant) then
       c = stencil%c
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
               c(1) * box_p%cc(i_c1, iv) + &
               c(2) * box_p%cc(i_c2, iv)
       end do
    else
       do i = 1, nc
          c = stencil%v(:, IJK)
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
               c(1) * box_p%cc(i_c1, iv) + &
               c(2) * box_p%cc(i_c2, iv)
       end do
    end if
#elif NDIM == 2
    if (stencil%stype == stencil_constant) then
       c = stencil%c
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                  c(1) * box_p%cc(i_c1, j_c1, iv) + &
                  c(2) * box_p%cc(i_c2, j_c1, iv) + &
                  c(3) * box_p%cc(i_c1, j_c2, iv)
          end do
       end do
    else
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             c = stencil%v(:, IJK)
             box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                  c(1) * box_p%cc(i_c1, j_c1, iv) + &
                  c(2) * box_p%cc(i_c2, j_c1, iv) + &
                  c(3) * box_p%cc(i_c1, j_c2, iv)
          end do
       end do
    end if
#elif NDIM == 3
    if (stencil%stype == stencil_constant) then
       c = stencil%c
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             do i = 1, nc
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
                box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                     c(1) * box_p%cc(i_c1, j_c1, k_c1, iv) + &
                     c(2) * box_p%cc(i_c2, j_c1, k_c1, iv) + &
                     c(3) * box_p%cc(i_c1, j_c2, k_c1, iv) + &
                     c(4) * box_p%cc(i_c1, j_c1, k_c2, iv)
             end do
          end do
       end do
    else
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             do i = 1, nc
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
                c = stencil%v(:, IJK)
                box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                     c(1) * box_p%cc(i_c1, j_c1, k_c1, iv) + &
                     c(2) * box_p%cc(i_c2, j_c1, k_c1, iv) + &
                     c(3) * box_p%cc(i_c1, j_c2, k_c1, iv) + &
                     c(4) * box_p%cc(i_c1, j_c1, k_c2, iv)
             end do
          end do
       end do
    end if
#endif
  end subroutine stencil_prolong_234

  !> (Bi-/tri-)linear prolongation using nearest 2**NDIM neighbors
  subroutine stencil_prolong_248(box_p, box_c, stencil, iv, iv_to, add)
    type(box_t), intent(in)     :: box_p      !< Parent box
    type(box_t), intent(inout)  :: box_c   !< Child box
    type(stencil_t), intent(in) :: stencil !< Stencil
    integer, intent(in)         :: iv      !< Input variable
    integer, intent(in)         :: iv_to   !< Destination variable
    logical, intent(in)         :: add     !< Add to old values

    integer  :: nc, ix_offset(NDIM), IJK
    integer  :: IJK_(c1), IJK_(c2)
    real(dp) :: c(2**NDIM)

    nc        = box_c%n_cell
    ix_offset = af_get_child_offset(box_c)
    if (stencil%stype == stencil_sparse) error stop "sparse not implemented"
    if (.not. add) box_c%cc(DTIMES(1:nc), iv_to) = 0

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
    if (stencil%stype == stencil_constant) then
       c = stencil%c
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
               c(1) * box_p%cc(i_c1, iv) + &
               c(2) * box_p%cc(i_c2, iv)
       end do
    else
       do i = 1, nc
          c = stencil%v(:, IJK)
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
               c(1) * box_p%cc(i_c1, iv) + &
               c(2) * box_p%cc(i_c2, iv)
       end do
    end if
#elif NDIM == 2
    if (stencil%stype == stencil_constant) then
       c = stencil%c
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                  c(1) * box_p%cc(i_c1, j_c1, iv) + &
                  c(2) * box_p%cc(i_c2, j_c1, iv) + &
                  c(3) * box_p%cc(i_c1, j_c2, iv) + &
                  c(4) * box_p%cc(i_c2, j_c2, iv)
          end do
       end do
    else
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             c = stencil%v(:, IJK)
             box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                  c(1) * box_p%cc(i_c1, j_c1, iv) + &
                  c(2) * box_p%cc(i_c2, j_c1, iv) + &
                  c(3) * box_p%cc(i_c1, j_c2, iv) + &
                  c(4) * box_p%cc(i_c2, j_c2, iv)
          end do
       end do
    end if
#elif NDIM == 3
    if (stencil%stype == stencil_constant) then
       c = stencil%c
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             do i = 1, nc
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
                box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                     c(1) * box_p%cc(i_c1, j_c1, k_c1, iv) + &
                     c(2) * box_p%cc(i_c2, j_c1, k_c1, iv) + &
                     c(3) * box_p%cc(i_c1, j_c2, k_c1, iv) + &
                     c(4) * box_p%cc(i_c2, j_c2, k_c1, iv) + &
                     c(5) * box_p%cc(i_c1, j_c1, k_c2, iv) + &
                     c(6) * box_p%cc(i_c2, j_c1, k_c2, iv) + &
                     c(7) * box_p%cc(i_c1, j_c2, k_c2, iv) + &
                     c(8) * box_p%cc(i_c2, j_c2, k_c2, iv)
             end do
          end do
       end do
    else
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             do i = 1, nc
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
                c = stencil%v(:, IJK)
                box_c%cc(IJK, iv_to) = box_c%cc(IJK, iv_to) + &
                     c(1) * box_p%cc(i_c1, j_c1, k_c1, iv) + &
                     c(2) * box_p%cc(i_c2, j_c1, k_c1, iv) + &
                     c(3) * box_p%cc(i_c1, j_c2, k_c1, iv) + &
                     c(4) * box_p%cc(i_c2, j_c2, k_c1, iv) + &
                     c(5) * box_p%cc(i_c1, j_c1, k_c2, iv) + &
                     c(6) * box_p%cc(i_c2, j_c1, k_c2, iv) + &
                     c(7) * box_p%cc(i_c1, j_c2, k_c2, iv) + &
                     c(8) * box_p%cc(i_c2, j_c2, k_c2, iv)
             end do
          end do
       end do
    end if
#endif
  end subroutine stencil_prolong_248

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

    real(dp) :: c(2*NDIM+1), inv_c1
    integer  :: IJK, i0, nc
#if NDIM == 2
    real(dp) :: rfac(2, box%n_cell), c_cyl(2*NDIM+1)
    real(dp) :: cc_cyl(2*NDIM+1, box%n_cell), inv_cc1(box%n_cell)
#endif

    if (stencil%stype == stencil_sparse) error stop "sparse not implemented"
    nc = box%n_cell

    associate (cc => box%cc, nc => box%n_cell)
      if (allocated(stencil%bc_correction)) then
         cc(DTIMES(1:nc), i_rhs) = cc(DTIMES(1:nc), i_rhs) + &
              stencil%bc_correction
      end if

#if NDIM == 1
      i0 = 2 - iand(redblack, 1)
      if (stencil%stype == stencil_constant) then
         c = stencil%c
         inv_c1 = 1 / c(1)

         do i = i0, nc, 2
            cc(IJK, iv) = (cc(IJK, i_rhs) &
                 - c(2) * cc(i-1, iv) &
                 - c(3) * cc(i+1, iv)) * inv_c1
         end do
      else
         do i = i0, nc, 2
            c = stencil%v(:, IJK)
            cc(IJK, iv) = (cc(IJK, i_rhs) &
                 - c(2) * cc(i-1, iv) &
                 - c(3) * cc(i+1, iv)) / c(1)
         end do
      end if
#elif NDIM == 2
      if (stencil%cylindrical_gradient) then
         ! Correct for cylindrical coordinates, assuming the elements correspond
         ! to a gradient operation
         call af_cyl_flux_factors(box, rfac)

         if (stencil%stype == stencil_constant) then
            c = stencil%c

            ! Pre-compute coefficients for each i-index
            do i = 1, nc
               cc_cyl(2:3, i) = rfac(1:2, i) * c(2:3)
               cc_cyl(1, i)   = c(1) - (cc_cyl(2, i) - c(2)) &
                    - (cc_cyl(3, i) - c(3))
               cc_cyl(4:, i)  = c(4:)
               inv_cc1(i)     = 1 / cc_cyl(1, i)
            end do

            do j = 1, nc
               i0 = 2 - iand(ieor(redblack, j), 1)
               do i = i0, nc, 2
                  cc(IJK, iv) = (cc(IJK, i_rhs) &
                       - cc_cyl(2, i) * cc(i-1, j, iv) &
                       - cc_cyl(3, i) * cc(i+1, j, iv) &
                       - cc_cyl(4, i) * cc(i, j-1, iv) &
                       - cc_cyl(5, i) * cc(i, j+1, iv)) * inv_cc1(i)
               end do
            end do
         else
            ! Variable stencil
            do j = 1, nc
               i0 = 2 - iand(ieor(redblack, j), 1)
               do i = i0, nc, 2
                  c = stencil%v(:, IJK)
                  c_cyl(2:3) = rfac(1:2, i) * c(2:3)
                  c_cyl(1) = c(1) - (c_cyl(2) - c(2)) - (c_cyl(3) - c(3))
                  c_cyl(4:) = c(4:)

                  cc(IJK, iv) = (cc(IJK, i_rhs) &
                       - c_cyl(2) * cc(i-1, j, iv) &
                       - c_cyl(3) * cc(i+1, j, iv) &
                       - c_cyl(4) * cc(i, j-1, iv) &
                       - c_cyl(5) * cc(i, j+1, iv)) / c_cyl(1)
               end do
            end do
         end if
      else                      ! No cylindrical gradient correction
         if (stencil%stype == stencil_constant) then
            c = stencil%c
            inv_c1 = 1 / c(1)

            do j = 1, nc
               i0 = 2 - iand(ieor(redblack, j), 1)
               do i = i0, nc, 2
                  cc(IJK, iv) = (cc(IJK, i_rhs) &
                       - c(2) * cc(i-1, j, iv) &
                       - c(3) * cc(i+1, j, iv) &
                       - c(4) * cc(i, j-1, iv) &
                       - c(5) * cc(i, j+1, iv)) * inv_c1
               end do
            end do
         else
            do j = 1, nc
               i0 = 2 - iand(ieor(redblack, j), 1)
               do i = i0, nc, 2
                  c = stencil%v(:, IJK)
                  cc(IJK, iv) = (cc(IJK, i_rhs) &
                       - c(2) * cc(i-1, j, iv) &
                       - c(3) * cc(i+1, j, iv) &
                       - c(4) * cc(i, j-1, iv) &
                       - c(5) * cc(i, j+1, iv)) / c(1)
               end do
            end do
         end if
      end if
#elif NDIM == 3
      if (stencil%stype == stencil_constant) then
         c = stencil%c
         inv_c1 = 1 / c(1)

         do k = 1, nc
            do j = 1, nc
               i0 = 2 - iand(ieor(redblack, k+j), 1)
               do i = i0, nc, 2
                  cc(IJK, iv) = (cc(IJK, i_rhs) &
                       - c(2) * cc(i-1, j, k, iv) &
                       - c(3) * cc(i+1, j, k, iv) &
                       - c(4) * cc(i, j-1, k, iv) &
                       - c(5) * cc(i, j+1, k, iv) &
                       - c(6) * cc(i, j, k-1, iv) &
                       - c(7) * cc(i, j, k+1, iv)) * inv_c1
               end do
            end do
         end do
      else
         do k = 1, nc
            do j = 1, nc
               i0 = 2 - iand(ieor(redblack, k+j), 1)
               do i = i0, nc, 2
                  c = stencil%v(:, IJK)
                  cc(IJK, iv) = (cc(IJK, i_rhs) &
                       - c(2) * cc(i-1, j, k, iv) &
                       - c(3) * cc(i+1, j, k, iv) &
                       - c(4) * cc(i, j-1, k, iv) &
                       - c(5) * cc(i, j+1, k, iv) &
                       - c(6) * cc(i, j, k-1, iv) &
                       - c(7) * cc(i, j, k+1, iv)) / c(1)
               end do
            end do
         end do
      end if
#endif

      if (allocated(stencil%bc_correction)) then
         cc(DTIMES(1:nc), i_rhs) = cc(DTIMES(1:nc), i_rhs) - &
              stencil%bc_correction
      end if
    end associate
  end subroutine stencil_gsrb_357

  !> Convert a variable stencil to constant one if possible
  subroutine af_stencil_try_constant(box, ix, abs_tol, success)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: ix       !< Stencil index
    real(dp), intent(in)       :: abs_tol  !< Absolute tolerance
    logical                    :: success !< Whether the stencil was converted
    integer                    :: IJK, nc, n_coeff

    if (box%stencils(ix)%stype == stencil_sparse) &
         error stop "sparse not implemented"
    if (.not. allocated(box%stencils(ix)%v)) &
         error stop "No variable stencil present"

    nc = box%n_cell
    success = .false.

    ! Check if all coefficients are the same, otherwise return
    do KJI_DO(1, nc)
       if (any(abs(box%stencils(ix)%v(:, IJK) - &
            box%stencils(ix)%v(:, DTIMES(1))) > abs_tol)) return
    end do; CLOSE_DO

    box%stencils(ix)%stype = stencil_constant
    n_coeff = size(box%stencils(ix)%v, 1)
    allocate(box%stencils(ix)%c(n_coeff))
    box%stencils(ix)%c = box%stencils(ix)%v(:, DTIMES(1))
    deallocate(box%stencils(ix)%v)
    success = .true.
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

    if (stencil%stype == stencil_constant) then
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
