#include "cpp_macros.h"
!> This module contains routines for including dielectrics
module m_dielectric
  use m_af_types

  implicit none
  private

  type surface_t
     logical               :: in_use    !< Whether the surface is active
     integer               :: id_in     !< Id of box inside dielectric
     integer               :: id_out    !< Id of box outside dielectric
     !> Which neighbor of the outside box is the dielectric
     integer               :: direction
     integer               :: ix_parent !< Index of parent surface
     real(dp)              :: eps       !< Permittivity of dielectric
     integer               :: offset_parent(NDIM-1) !< Index of parent surface
     !> Surface densities
#if NDIM == 2
     real(dp), allocatable :: sd(:, :)
#elif NDIM == 3
     real(dp), allocatable :: sd(:, :, :)
#endif
  end type surface_t

  !> Value indicating there is no surface
  integer, parameter :: no_surface = -1

  type dielectric_t
     !> Number of variables to store on the surface
     integer :: n_variables = 0
     !> Size of boxes (and thus also surfaces)
     integer :: n_cell = -1
     !> Index of epsilon in tree
     integer :: i_eps = -1
     !> Index of electric potential in tree
     integer :: i_phi = -1
     !> Whether the electric field is stored at cell faces
     logical :: field_at_faces = .true.
     !> Index of face-centered electric field variable
     integer :: ifc_E = -1
     !> Indices of cell-centered electric field variables
     integer :: icc_E(NDIM) = -1
     !> Highest surface id
     integer :: max_ix = 0
     !> Maximum number of surfaces
     integer :: surface_limit = -1

     !> List of all the surfaces
     type(surface_t), allocatable :: surfaces(:)
     !> Number of removed surfaces
     integer                      :: n_removed = 0
     !> Indices of removed surfaces
     integer, allocatable         :: removed_surfaces(:)
     !> Mapping of boxes to surfaces
     integer, allocatable         :: box_id_to_surface_ix(:)
  end type dielectric_t

  interface
     function eps_func(x) result(eps)
       import
       real(dp), intent(in) :: x(NDIM)
       real(dp)             :: eps
     end function eps_func

     function value_func(x) result(my_val)
       import
       real(dp), intent(in) :: x(NDIM)
       real(dp)             :: my_val
     end function value_func
  end interface

  public :: dielectric_t
  public :: dielectric_initialize
  public :: dielectric_set_values
  public :: dielectric_update_after_refinement
  public :: dielectric_surface_charge_to_rhs
  public :: dielectric_correct_field_fc

contains

  subroutine dielectric_initialize(tree, dielectric, n_variables)
    use m_af_ghostcell, only: af_gc_tree
    type(af_t), intent(inout)         :: tree
    type(dielectric_t), intent(inout) :: dielectric
    integer, intent(in)               :: n_variables
    integer                           :: nc, id, i, ix, i_eps
    integer                           :: nb, nb_id
    real(dp)                          :: curr_eps, nb_eps

    if (.not. tree%ready) error stop "Tree not ready"

    if (tree%highest_lvl > 1) then
       ! Check that there are no partially refined grid levels
       if (size(tree%lvls(tree%highest_lvl-1)%leaves) > 0) &
            error stop "Tree not uniformly refined"
    end if

    if (dielectric%i_eps < 0) error stop "dielectric%i_eps not set"

    i_eps = dielectric%i_eps
    nc    = tree%n_cell

    dielectric%max_ix = 0
    dielectric%n_removed       = 0
    dielectric%n_variables     = n_variables
    dielectric%n_cell          = tree%n_cell
    dielectric%surface_limit   = tree%box_limit / 10 ! Conservative estimate

    allocate(dielectric%surfaces(dielectric%surface_limit))
    allocate(dielectric%removed_surfaces(dielectric%surface_limit))
    allocate(dielectric%box_id_to_surface_ix(tree%box_limit))
    dielectric%box_id_to_surface_ix = no_surface

    ! Construct a list of surfaces
    do i = 1, size(tree%lvls(tree%highest_lvl)%ids)
       id = tree%lvls(tree%highest_lvl)%ids(i)

       ! Get eps on the interior of the current box
       curr_eps = tree%boxes(id)%cc(DTIMES(1), i_eps)

       do nb = 1, af_num_neighbors
          nb_id = tree%boxes(id)%neighbors(nb)

          ! Skip physical boundaries
          if (nb_id <= af_no_box) cycle

          ! Get eps on the interior of the neighbour box
          nb_eps = tree%boxes(nb_id)%cc(DTIMES(1), i_eps)

          if (nb_eps > curr_eps) then
             ! There is a surface
             ix = get_new_surface_id(dielectric)

             dielectric%surfaces(ix)%in_use    = .true.
             dielectric%surfaces(ix)%ix_parent = no_surface
             dielectric%surfaces(ix)%offset_parent(:) = 0
             dielectric%surfaces(ix)%id_in     = nb_id
             dielectric%surfaces(ix)%id_out    = id
             dielectric%surfaces(ix)%direction = nb
             dielectric%surfaces(ix)%eps       = nb_eps

             ! TODO: only a single surface per box for now
             dielectric%box_id_to_surface_ix(id) = ix
          end if
       end do
    end do

  end subroutine dielectric_initialize

  function get_new_surface_id(dielectric) result(ix)
    type(dielectric_t), intent(inout) :: dielectric
    logical                           :: use_removed
    integer                           :: ix, nc

    !$omp critical
    use_removed = (dielectric%n_removed > 0)
    if (use_removed) then
       ix = dielectric%removed_surfaces(dielectric%n_removed)
       dielectric%n_removed = dielectric%n_removed - 1
    else
       dielectric%max_ix = dielectric%max_ix + 1
       ix = dielectric%max_ix
    end if
    !$omp end critical

    if (.not. use_removed) then
       nc = dielectric%n_cell
#if NDIM == 2
       allocate(dielectric%surfaces(ix)%sd(nc, dielectric%n_variables))
#elif NDIM == 3
       allocate(dielectric%surfaces(ix)%sd(nc, nc, dielectric%n_variables))
#endif
       dielectric%surfaces(ix)%sd = 0.0_dp
    end if
  end function get_new_surface_id

  subroutine dielectric_set_values(tree, dielectric, iv, user_func)
    use m_af_ghostcell, only: af_gc_get_boundary_coords
    type(af_t), intent(in)            :: tree
    type(dielectric_t), intent(inout) :: dielectric
    integer, intent(in)               :: iv !< Surface variable
    procedure(value_func)             :: user_func !< User supplied function
    integer                           :: i, ix, id_out
    real(dp) :: coords(NDIM, tree%n_cell**(NDIM-1))

    ! Loop over the surfaces and call the user function to set values
    do ix = 1, dielectric%max_ix
       if (dielectric%surfaces(ix)%in_use) then
          id_out = dielectric%surfaces(ix)%id_out
          call af_gc_get_boundary_coords(tree%boxes(id_out), &
               dielectric%surfaces(ix)%direction, coords)
#if NDIM == 2
          do i = 1, tree%n_cell
             dielectric%surfaces(ix)%sd(i, iv) = user_func(coords(:, i))
          end do
#elif NDIM == 3
          error stop
#endif
       end if
    end do
  end subroutine dielectric_set_values

  subroutine dielectric_update_after_refinement(tree, dielectric, ref_info)
    type(af_t), intent(in)            :: tree
    type(dielectric_t), intent(inout) :: dielectric
    type(ref_info_t), intent(in)      :: ref_info
    integer                           :: lvl, i, id, p_id, ix, p_ix
    integer                           :: nc, id_in, id_out

    nc = tree%n_cell

    ! Handle removed surfaces
    do i = 1, size(ref_info%rm)
       id = ref_info%rm(i)
       ix = dielectric%box_id_to_surface_ix(id)
       if (ix > no_surface) then
          call restrict_surface_to_parent(tree, dielectric, ix)

          !$omp critical
          dielectric%n_removed = dielectric%n_removed + 1
          dielectric%removed_surfaces(dielectric%n_removed) = ix
          !$omp end critical
          id_in = dielectric%surfaces(ix)%id_in
          id_out = dielectric%surfaces(ix)%id_out
          dielectric%box_id_to_surface_ix([id_in, id_out]) = no_surface
          dielectric%surfaces(ix)%in_use = .false.
       end if
    end do

    ! Add new surfaces
    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id   = ref_info%lvls(lvl)%add(i)

          ! Check if parent has a surface
          p_id = tree%boxes(id)%parent
          p_ix = dielectric%box_id_to_surface_ix(p_id)

          if (p_ix > no_surface) then
             ! We only need to prolong once per parent surface, so check if it
             ! is the first child outside the dielectric
             if (p_id == dielectric%surfaces(p_ix)%id_out .and. &
                  af_ix_to_ichild(tree%boxes(id)%ix) == 1) then
                call prolong_surface_from_parent(tree, dielectric, p_ix, p_id)
             end if
          end if
       end do
    end do

  end subroutine dielectric_update_after_refinement

  subroutine prolong_surface_from_parent(tree, dielectric, p_ix, p_id)
    type(af_t), intent(in)            :: tree
    type(dielectric_t), intent(inout) :: dielectric
    integer, intent(in)               :: p_ix !< Index of parent surface
    integer, intent(in)               :: p_id !< Index of parent box (on outside)
    integer                           :: i, n, ix, ix_offset(NDIM)
    integer                           :: nc, dix(NDIM-1), direction
    integer                           :: i_child, id_child, id_in, id_out

    nc = tree%n_cell

    do n = 1, size(af_child_adj_nb, 1)
       ix = get_new_surface_id(dielectric)
       direction = dielectric%surfaces(p_ix)%direction

       i_child  = af_child_adj_nb(n, direction)
       id_child = tree%boxes(p_id)%children(i_child)
       id_out   = tree%boxes(p_id)%children(i_child)
       id_in    = tree%boxes(id_child)%neighbors(direction)

       dielectric%surfaces(ix)%in_use    = .true.
       dielectric%surfaces(ix)%ix_parent = p_ix
       dielectric%surfaces(ix)%id_in     = id_in
       dielectric%surfaces(ix)%id_out    = id_out
       dielectric%surfaces(ix)%direction = direction
       dielectric%surfaces(ix)%eps       = dielectric%surfaces(p_ix)%eps
       dielectric%box_id_to_surface_ix([id_out, id_in]) = ix

       ix_offset = af_get_child_offset(tree%boxes(id_child))

       ! Compute index offset of child surface relative to parent. This is nc/2
       ! times the child offset, but only in the dimensions parallel to the
       ! surface.
       dix = nc/2 * pack(af_child_dix(:, i_child), &
            [(i, i=1,NDIM)] /= af_neighb_dim(direction))

       associate (sd => dielectric%surfaces(ix)%sd, &
            sd_p => dielectric%surfaces(p_ix)%sd)
#if NDIM == 2
         ! Copy the values from the parent
         sd(1:nc:2, :) = sd_p(dix(1)+1:dix(1)+nc/2, :)
         sd(2:nc:2, :) = sd(1:nc:2, :)
#elif NDIM == 3
         error stop
#endif
       end associate
    end do

    dielectric%surfaces(p_ix)%in_use = .false.

  end subroutine prolong_surface_from_parent

  subroutine restrict_surface_to_parent(tree, dielectric, ix)
    type(af_t), intent(in)            :: tree
    type(dielectric_t), intent(inout) :: dielectric
    integer, intent(in)               :: ix
    integer                           :: p_ix, nc, dix(NDIM-1)

    p_ix = dielectric%surfaces(ix)%ix_parent
    dix  = dielectric%surfaces(ix)%offset_parent
    nc   = tree%n_cell

    associate (sd => dielectric%surfaces(ix)%sd, &
         sd_p => dielectric%surfaces(p_ix)%sd)
#if NDIM == 2
      ! Average the value on the children
      sd_p(dix(1)+1:dix(1)+nc/2, 1) = 0.5_dp * (sd(1:nc:2, 1) + sd(2:nc:2, 1))
#elif NDIM == 3
      error stop
#endif
    end associate

    dielectric%surfaces(p_ix)%in_use = .true.

  end subroutine restrict_surface_to_parent

  ! subroutine dielectric_fix_refine(tree, dielectric, ref_flags)
  !   type(af_t), intent(in) :: tree
  !   integer, intent(inout) :: ref_flags(:)
  !   integer                :: lvl, i, id, nb_id, ix

  !   do lvl = 1, tree%highest_lvl
  !      do i = 1, size(tree%lvls(lvl)%ids)
  !         id = tree%lvls(lvl)%ids(i)

  !         ix = box_id_to_surface_ix(id)
  !         if (ix > 0) then
  !            if (id == surfaces(ix)%id_gas) then
  !               nb_id = surfaces(ix)%id_diel
  !               ref_flags([id, nb_id]) = &
  !                    maxval([ref_flags([id, nb_id]), af_keep_ref])
  !            end if
  !         end if
  !      end do
  !   end do
  ! end subroutine dielectric_fix_refine

  subroutine dielectric_surface_charge_to_rhs(tree, dielectric, &
       i_sigma, i_rhs, fac)
    type(af_t), intent(inout)      :: tree
    type(dielectric_t), intent(in) :: dielectric
    integer, intent(in)            :: i_sigma !< Surface charage variable
    integer, intent(in)            :: i_rhs   !< Rhs variable (in the tree)
    real(dp), intent(in)           :: fac     !< Multiplication factor
    integer                        :: n

    do n = 1, dielectric%max_ix
       if (.not. dielectric%surfaces(n)%in_use) cycle

       call surface_charge_to_rhs(tree%boxes, dielectric%surfaces(n), &
            i_sigma, i_rhs, fac)
    end do
  end subroutine dielectric_surface_charge_to_rhs

  subroutine surface_charge_to_rhs(boxes, surface, i_sigma, i_rhs, fac)
    type(box_t), intent(inout)  :: boxes(:)
    type(surface_t), intent(in) :: surface
    integer, intent(in)         :: i_sigma !< Surface charage variable
    integer, intent(in)         :: i_rhs   !< Rhs variable (in the tree)
    real(dp), intent(in)        :: fac     !< Multiplication factor
    integer                     :: nb, nc, id_in, id_out
    integer                     :: glo(NDIM), ghi(NDIM)
    integer                     :: dlo(NDIM), dhi(NDIM)
    real(dp)                    :: frac_gas, dr, fac_to_volume

    nb     = surface%direction
    id_in  = surface%id_in
    id_out = surface%id_out
    nc     = boxes(id_out)%n_cell
    dr     = boxes(id_out)%dr(af_neighb_dim(nb))

    ! Factor to convert surface charge density to volume density
    fac_to_volume = 1 / dr

    ! Get index range for cells adjacent to the boundary of the dielectric box
    call af_get_index_bc_inside(af_neighb_rev(nb), nc, 1, dlo, dhi)

    ! Get similar index range for gas-phase box
    call af_get_index_bc_inside(nb, nc, 1, glo, ghi)

    ! How much of the rhs to put on the gas and dielectric side
    frac_gas  = 1.0_dp / (1.0_dp + surface%eps)

#if NDIM == 2
    boxes(id_out)%cc(glo(1):ghi(1), glo(2):ghi(2), i_rhs) = &
         boxes(id_out)%cc(glo(1):ghi(1), glo(2):ghi(2), i_rhs) &
         + frac_gas * fac * fac_to_volume * &
         reshape(surface%sd(:, i_sigma), [ghi(1)-glo(1)+1, ghi(2)-glo(2)+1])

    boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_rhs) = &
         boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_rhs) &
         + (1-frac_gas) * fac * fac_to_volume * &
         reshape(surface%sd(:, i_sigma), [dhi(1)-dlo(1)+1, dhi(2)-dlo(2)+1])
#elif NDIM == 3
    error stop
#endif

  end subroutine surface_charge_to_rhs

  !   !> Recompute the electric field near the surface. This routine should be used
  !   !> for electric fields defined at cell centers.
  !   subroutine dielectric_correct_field_cc(tree, i_phi, i_Ex, ...)
  !     type(af_t), intent(inout)  :: tree
  !     integer                     :: nc
  !     real(dp)  :: fac_fld, fac_charge, inv_dr(NDIM)
  !     integer                 :: lvl, i, id

  !     if (.not. tree%ready) stop "dielectric_correct_field: set_base has not been called"

  !     do lvl = 1, tree%highest_lvl
  !        do i = 1, size(tree%lvls(lvl)%leaves)
  !           id = tree%lvls(lvl)%leaves(i)
  !           if (surf%box_id_to_direction(id) > 0) then ! If the box is adjacent to boundary
  !              associate(box => tree%boxes(id), cc => tree%boxes(id)%cc)
  !                fac_fld = 2 / (1 + dielectric_eps) ! Assuming eps_gas = 1
  !                fac_charge = 1 / (1 + dielectric_eps) ! Assuming eps_gas = 1
  !                inv_dr = 1 / box%dr
  !                nc     = box%n_cell
  ! #if NDIM == 2
  !                select case (surf%box_id_to_direction(id))
  !                   ! Depending on the case, correct Ex or Ey, and update E accordingly
  !                case (af_neighb_lowx)
  !                   cc(1, 1:nc, i_Ex)  = fac_fld * cc(0, 1:nc, i_eps) * inv_dr(1) * (cc(0, 1:nc, i_phi) - cc(1, 1:nc, i_phi))&
  !                        - fac_charge * surf%box_id_to_charge(id, :)
  !                   cc(1, 1:nc, i_E) = sqrt(cc(1, 1:nc, i_Ex)**2 + cc(1, 1:nc, i_Ey)**2)
  !                case (af_neighb_highx)
  !                   cc(nc, 1:nc, i_Ex)  = fac_fld * cc(nc+1, 1:nc, i_eps) * inv_dr(1) * (cc(nc, 1:nc, i_phi) - cc(nc+1, 1:nc, i_phi))&
  !                        + fac_charge * surf%box_id_to_charge(id, :)
  !                   cc(nc, 1:nc, i_E) = sqrt(cc(nc, 1:nc, i_Ex)**2 + cc(nc, 1:nc, i_Ey)**2)
  !                case (af_neighb_lowy)
  !                   cc(1:nc, 1, i_Ey)  = fac_fld * cc(1:nc, 0, i_eps) * inv_dr(2) * (cc(1:nc, 0, i_phi) - cc(1:nc, 1, i_phi))&
  !                        - fac_charge * surf%box_id_to_charge(id, :)
  !                   cc(1:nc, 1, i_E) = sqrt(cc(1:nc, 1, i_Ex)**2 + cc(1:nc, 1, i_Ey)**2)
  !                case (af_neighb_highy)
  !                   cc(1:nc, nc, i_Ey) = fac_fld * cc(1:nc, nc+1, i_eps) * inv_dr(2) * (cc(1:nc, nc, i_phi) - cc(1:nc, nc+1, i_phi))&
  !                        + fac_charge * surf%box_id_to_charge(id, :)
  !                   cc(1:nc, nc, i_E) = sqrt(cc(1:nc, nc, i_Ex)**2 + cc(1:nc, nc, i_Ey)**2)
  !                end select
  ! #elif NDIM == 3
  !                error stop
  ! #endif
  !              end associate
  !           end if
  !        end do
  !     end do
  !   end subroutine dielectric_correct_field_cc

  !> Correct the electric field at face centers, assuming it has already been
  !> computed for a dielectric permittivity of one everywhere
  subroutine dielectric_correct_field_fc(tree, dielectric, i_sigma, i_fld, eps0)
    type(af_t), intent(inout)      :: tree
    type(dielectric_t), intent(in) :: dielectric
    integer, intent(in)            :: i_sigma !< Surface charge variable
    integer, intent(in)            :: i_fld   !< Face-centered field variable
    real(dp), intent(in)           :: eps0    !< Dielectric permittivity (vacuum)
    integer                        :: id_in, id_out, nc, nb, ix
    real(dp)                       :: eps, fac_fld(2), fac_charge

    nc = tree%n_cell

    do ix = 1, dielectric%max_ix
       if (dielectric%surfaces(ix)%in_use) then
          nb = dielectric%surfaces(ix)%direction
          eps = dielectric%surfaces(ix)%eps
          id_out = dielectric%surfaces(ix)%id_out
          id_in = dielectric%surfaces(ix)%id_in

          fac_fld = [2 * eps, 2.0_dp] / (1 + eps)
          fac_charge = 1 / (eps0 * (1 + eps))

          associate (fc_out => tree%boxes(id_out)%fc, &
               fc_in => tree%boxes(id_in)%fc, &
               sd => dielectric%surfaces(ix)%sd)
            select case (nb)
#if NDIM == 2
            case (af_neighb_lowx)
               fc_out(1, 1:nc, 1, i_fld) = fac_fld(1) * &
                    fc_out(1, 1:nc, 1, i_fld) - fac_charge * sd(:, i_sigma)
               fc_in(nc+1, 1:nc, 1, i_fld) = fac_fld(2) * &
                    fc_in(nc+1, 1:nc, 1, i_fld) + fac_charge * sd(:, i_sigma)
            case (af_neighb_highx)
               fc_out(nc+1, 1:nc, 1, i_fld) = fac_fld(1) * &
                    fc_out(nc+1, 1:nc, 1, i_fld) + fac_charge * sd(:, i_sigma)
               fc_in(1, 1:nc, 1, i_fld) = fac_fld(2) * &
                    fc_in(1, 1:nc, 1, i_fld) - fac_charge * sd(:, i_sigma)
            case (af_neighb_lowy)
               fc_out(1:nc, 1, 2, i_fld) = fac_fld(1) * &
                    fc_out(1:nc, 1, 2, i_fld) - fac_charge * sd(:, i_sigma)
               fc_in(1:nc, nc+1, 2, i_fld) = fac_fld(2) * &
                    fc_in(1:nc, nc+1, 2, i_fld) + fac_charge * sd(:, i_sigma)
            case (af_neighb_highy)
               fc_out(1:nc, nc+1, 2, i_fld) = fac_fld(1) * &
                    fc_out(1:nc, nc+1, 2, i_fld) + fac_charge * sd(:, i_sigma)
               fc_in(1:nc, 1, 2, i_fld) = fac_fld(2) * &
                    fc_in(1:nc, 1, 2, i_fld) - fac_charge * sd(:, i_sigma)
#elif NDIM == 3
            case default
               error stop
#endif
            end select
          end associate
       end if
    end do

  end subroutine dielectric_correct_field_fc

end module m_dielectric
