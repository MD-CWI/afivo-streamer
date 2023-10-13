#include "cpp_macros.h"
!> This module contains routines for including flat surfaces between changes in
!> epsilon (some material property). This can for example be used to include
!> flat dielectrics in electrostatic computations.
module m_af_surface
  use m_af_types

  implicit none
  private

  !> Type for a single surface
  type surface_t
     logical               :: in_use    !< Whether the surface is active
     integer               :: id_in     !< Id of box inside dielectric
     integer               :: id_out    !< Id of box outside dielectric
     !> Which neighbor of the outside box is the dielectric
     integer               :: direction
     integer               :: ix_parent !< Index of parent surface
     real(dp)              :: eps       !< Permittivity of dielectric
     integer               :: offset_parent(NDIM-1) !< Index of parent surface
     real(dp)              :: dr(NDIM-1) !< Grid spacing on surface
     !> Surface densities
     real(dp), allocatable :: sd(DTIMES(:))
  end type surface_t

  !> Value indicating there is no surface
  integer, public, parameter :: surface_none = -1

  !> Type for storing all the surfaces on a mesh
  type surfaces_t
     !> Whether the dielectric is initialized
     logical :: initialized = .false.
     !> Number of variables to store on the surface
     integer :: n_variables = 0
     !> Size of boxes (and thus also surfaces)
     integer :: n_cell = -1
     !> Index of epsilon in tree
     integer :: i_eps = -1
     !> Index of electric potential in tree
     integer :: i_phi = -1
     !> Highest surface id
     integer :: max_ix = 0
     !> Maximum number of surfaces
     integer :: surface_limit = -1
     !> Whether there is a 2D cylindrical geometry
     logical :: cylindrical = .false.

     !> List of all the surfaces
     type(surface_t), allocatable :: surfaces(:)
     !> Number of removed surfaces
     integer                      :: n_removed = 0
     !> Indices of removed surfaces
     integer, allocatable         :: removed_surfaces(:)
     !> Mapping of boxes next to a dielectric to surfaces
     integer, allocatable         :: box_id_out_to_surface_ix(:)
     !> Mapping of boxes inside a dielectric to surfaces
     integer, allocatable         :: box_id_in_to_surface_ix(:)
  end type surfaces_t

  interface
     function value_func(x) result(my_val)
       import
       real(dp), intent(in) :: x(NDIM)
       real(dp)             :: my_val
     end function value_func
  end interface

  public :: surface_t
  public :: surfaces_t
  public :: surface_initialize
  public :: surface_set_values
  public :: surface_set_weighted_sum
  public :: surface_copy_variable
  public :: surface_get_integral
  public :: surface_update_after_refinement
  public :: surface_inside_layer_to_surface
  public :: surface_surface_charge_to_rhs
  public :: surface_correct_field_fc
  public :: surface_get_refinement_links
  public :: surface_get_surface_cell
  public :: surface_write_output
#if NDIM == 2
  public :: surface_correct_field_cc
#endif

contains

  !> Initialize a set of surfaces based on the value of epsilon
  subroutine surface_initialize(tree, i_eps, diel, n_variables)
    use m_af_ghostcell, only: af_gc_tree
    type(af_t), intent(inout)         :: tree        !< Initialized grid
    integer, intent(in)               :: i_eps       !< Which variable stores epsilon
    type(surfaces_t), intent(inout) :: diel        !< The dielectric surface
    integer, intent(in)               :: n_variables !< Number of surface variables
    integer                           :: nc, id, i, j, ix
    integer                           :: nb, nb_id
    real(dp)                          :: curr_eps, nb_eps

    if (.not. tree%ready) error stop "Tree not ready"

    if (tree%highest_lvl > 1) then
       ! Check that there are no partially refined grid levels
       if (size(tree%lvls(tree%highest_lvl-1)%leaves) > 0) &
            error stop "Tree not uniformly refined"
    end if

    diel%initialized   = .true.
    diel%cylindrical   = (tree%coord_t == af_cyl)
    diel%i_eps         = i_eps
    nc                 = tree%n_cell
    diel%max_ix        = 0
    diel%n_removed     = 0
    diel%n_variables   = n_variables
    diel%n_cell        = tree%n_cell
    ! Assume that at most 1/10th of boxes has a surface
    diel%surface_limit = tree%box_limit / 10

    allocate(diel%surfaces(diel%surface_limit))
    allocate(diel%removed_surfaces(diel%surface_limit))
    allocate(diel%box_id_out_to_surface_ix(tree%box_limit))
    allocate(diel%box_id_in_to_surface_ix(tree%box_limit))
    diel%box_id_out_to_surface_ix = surface_none
    diel%box_id_in_to_surface_ix = surface_none

    ! Construct a list of surfaces
    do i = 1, size(tree%lvls(tree%highest_lvl)%ids)
       id = tree%lvls(tree%highest_lvl)%ids(i)

       ! Check whether epsilon is set consistently on this box
       if (maxval(tree%boxes(id)%cc(DTIMES(1:nc), i_eps)) > &
            minval(tree%boxes(id)%cc(DTIMES(1:nc), i_eps))) &
            error stop "epsilon not uniform on a box"

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
             ix = get_new_surface_ix(diel)

             diel%surfaces(ix)%in_use    = .true.
             diel%surfaces(ix)%ix_parent = surface_none
             diel%surfaces(ix)%offset_parent(:) = 0
             diel%surfaces(ix)%id_in     = nb_id
             diel%surfaces(ix)%id_out    = id
             diel%surfaces(ix)%direction = nb
             diel%surfaces(ix)%eps       = nb_eps

             ! Extract the grid spacing parallel to the surface
             diel%surfaces(ix)%dr = &
                  pack(tree%boxes(id)%dr, [(j, j=1,NDIM)] /= af_neighb_dim(nb))

             if (diel%box_id_out_to_surface_ix(id) /= surface_none) &
                  error stop "box with multiple adjacent surfaces"
             if (diel%box_id_in_to_surface_ix(nb_id) /= surface_none) &
                  error stop "box with multiple adjacent surfaces"
             diel%box_id_out_to_surface_ix(id) = ix
             diel%box_id_in_to_surface_ix(nb_id) = ix
          end if
       end do
    end do

  end subroutine surface_initialize

  !> Get index for new surface
  function get_new_surface_ix(diel) result(ix)
    type(surfaces_t), intent(inout) :: diel
    logical                           :: use_removed
    integer                           :: ix, nc

    !$omp critical
    use_removed = (diel%n_removed > 0)
    if (use_removed) then
       ix = diel%removed_surfaces(diel%n_removed)
       diel%n_removed = diel%n_removed - 1
    else
       diel%max_ix = diel%max_ix + 1
       ix = diel%max_ix
    end if
    !$omp end critical

    if (.not. use_removed) then
       nc = diel%n_cell
#if NDIM == 1
       allocate(diel%surfaces(ix)%sd(diel%n_variables))
#elif NDIM == 2
       allocate(diel%surfaces(ix)%sd(nc, diel%n_variables))
#elif NDIM == 3
       allocate(diel%surfaces(ix)%sd(nc, nc, diel%n_variables))
#endif
       diel%surfaces(ix)%sd = 0.0_dp
    end if
  end function get_new_surface_ix

  !> Set values on a dielectric with a user-defined function
  subroutine surface_set_values(tree, diel, iv, user_func)
    type(af_t), intent(in)            :: tree
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)               :: iv        !< Surface variable
    procedure(value_func)             :: user_func !< User supplied function
    integer                           :: ix, id_out
#if NDIM > 1
    integer                           :: i
#endif
#if NDIM == 3
    integer                           :: j, n
#endif
    real(dp)                          :: coords(NDIM, tree%n_cell**(NDIM-1))

    if (.not. diel%initialized) error stop "dielectric not initialized"

    ! Loop over the surfaces and call the user function to set values
    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
          id_out = diel%surfaces(ix)%id_out
          call af_get_face_coords(tree%boxes(id_out), &
               diel%surfaces(ix)%direction, coords)
#if NDIM == 1
          diel%surfaces(ix)%sd(iv) = user_func(coords(1, 1))
#elif NDIM == 2
          do i = 1, tree%n_cell
             diel%surfaces(ix)%sd(i, iv) = user_func(coords(:, i))
          end do
#elif NDIM == 3
          n = 0
          do j = 1, tree%n_cell
             do i = 1, tree%n_cell
                n = n + 1
                diel%surfaces(ix)%sd(i, j, iv) = user_func(coords(:, n))
             end do
          end do
#endif
       end if
    end do
  end subroutine surface_set_values

  !> Set surface variable to a weighted sum of other variables
  subroutine surface_set_weighted_sum(diel, i_out, i_in, w_in)
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)               :: i_out   !< Output surface density
    integer, intent(in)               :: i_in(:) !< List of input surface densities
    real(dp), intent(in)              :: w_in(:) !< Weights of input densities
    integer                           :: ix
#if NDIM == 3
    integer                           :: ii
#endif

    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
#if NDIM == 2
          diel%surfaces(ix)%sd(:, i_out) = &
               matmul(diel%surfaces(ix)%sd(:, i_in), w_in)
#elif NDIM == 3
          ! TODO improve by removing do-loop??
          diel%surfaces(ix)%sd(:, :, i_out) = 0.0_dp
          do ii=1, size(i_in)
            diel%surfaces(ix)%sd(:, :, i_out) =  diel%surfaces(ix)%sd(:, :, i_out) &
              + w_in(ii) * diel%surfaces(ix)%sd(:, :, i_in(ii))
          end do
#endif
       end if
    end do
  end subroutine surface_set_weighted_sum

  !> Copy surface variable to another index
  subroutine surface_copy_variable(diel, i_in, i_out)
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)             :: i_in  !< Input surface density
    integer, intent(in)             :: i_out !< Output surface density
    integer                         :: ix

    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
#if NDIM == 2
          diel%surfaces(ix)%sd(:, i_out) = diel%surfaces(ix)%sd(:, i_in)
#elif NDIM == 3
          diel%surfaces(ix)%sd(:, :, i_out) = diel%surfaces(ix)%sd(:, :, i_in)
#endif
       end if
    end do
  end subroutine surface_copy_variable

  !> Compute integral of surface variable
  subroutine surface_get_integral(tree, diel, i_surf, surf_int)
    type(af_t), intent(in)       :: tree
    type(surfaces_t), intent(in) :: diel
    integer, intent(in)          :: i_surf   !< Surface variables
    real(dp), intent(out)        :: surf_int !< Surface integral
    integer                      :: ix
#if NDIM == 2
    real(dp), parameter          :: two_pi = 2 * acos(-1.0_dp)
    real(dp)                     :: coords(NDIM, diel%n_cell**(NDIM-1))
#endif

    surf_int = 0
    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
#if NDIM == 2
          if (diel%cylindrical) then
             ! Multiply surface elements with 2 * pi * r
             call af_get_face_coords(tree%boxes(diel%surfaces(ix)%id_out), &
                  diel%surfaces(ix)%direction, coords)
             surf_int = surf_int + two_pi * diel%surfaces(ix)%dr(1) * &
                  sum(diel%surfaces(ix)%sd(:, i_surf) * coords(1, :))
          else
             surf_int = surf_int + product(diel%surfaces(ix)%dr) * &
                  sum(diel%surfaces(ix)%sd(:, i_surf))
          end if
#elif NDIM == 3
          surf_int = surf_int + product(diel%surfaces(ix)%dr) * &
               sum(diel%surfaces(ix)%sd(:, :, i_surf))
#endif
       end if
    end do
  end subroutine surface_get_integral

  !> Update the dielectric surface after the mesh has been refined
  subroutine surface_update_after_refinement(tree, diel, ref_info)
    type(af_t), intent(in)            :: tree
    type(surfaces_t), intent(inout) :: diel
    type(ref_info_t), intent(in)      :: ref_info
    integer                           :: lvl, i, id, p_id, ix, p_ix, nc

    if (.not. diel%initialized) error stop "dielectric not initialized"
    nc = tree%n_cell

    ! Handle removed surfaces
    do i = 1, size(ref_info%rm)
       id = ref_info%rm(i)
       ix = diel%box_id_out_to_surface_ix(id)
       if (ix > surface_none) then
          call restrict_surface_to_parent(tree, diel, ix)
       end if
    end do

    ! Add new surfaces
    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id   = ref_info%lvls(lvl)%add(i)

          ! Check if parent has a surface
          p_id = tree%boxes(id)%parent
          p_ix = diel%box_id_out_to_surface_ix(p_id)

          if (p_ix > surface_none) then
             ! We only need to prolong once per parent surface
             if (af_ix_to_ichild(tree%boxes(id)%ix) == 1) then
                call prolong_surface_from_parent(tree, diel, p_ix, p_id)
             end if
          end if
       end do
    end do

  end subroutine surface_update_after_refinement

  !> Prolong a parent surface to newly created child surfaces
  subroutine prolong_surface_from_parent(tree, diel, p_ix, p_id)
    type(af_t), intent(in)            :: tree
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)               :: p_ix !< Index of parent surface
    integer, intent(in)               :: p_id !< Index of parent box (on outside)
    integer                           :: i, n, ix, ix_offset(NDIM)
    integer                           :: nc, dix(NDIM-1), direction
    integer                           :: i_child, id_in, id_out

    nc = tree%n_cell

    do n = 1, size(af_child_adj_nb, 1)
       ix = get_new_surface_ix(diel)
       direction = diel%surfaces(p_ix)%direction

       ! Select a child adjacent to the neighbor containing the dielectric
       i_child  = af_child_adj_nb(n, direction)
       id_out   = tree%boxes(p_id)%children(i_child)
       ! The neighbor inside the dielectric should exist
       id_in    = tree%boxes(id_out)%neighbors(direction)
       if (id_out <= af_no_box .or. id_in <= af_no_box) &
            error stop "Error in prolongation: children do not exist"

       diel%surfaces(ix)%in_use    = .true.
       diel%surfaces(ix)%ix_parent = p_ix
       diel%surfaces(ix)%id_in     = id_in
       diel%surfaces(ix)%id_out    = id_out
       diel%surfaces(ix)%direction = direction
       diel%surfaces(ix)%eps       = diel%surfaces(p_ix)%eps
       diel%box_id_out_to_surface_ix(id_out) = ix
       diel%box_id_in_to_surface_ix(id_in) = ix

       ix_offset = af_get_child_offset(tree%boxes(id_out))

       ! Compute index offset of child surface relative to parent. This is nc/2
       ! times the child offset, but only in the dimensions parallel to the
       ! surface.
       dix = nc/2 * pack(af_child_dix(:, i_child), &
            [(i, i=1,NDIM)] /= af_neighb_dim(direction))
       diel%surfaces(ix)%offset_parent(:) = dix

       ! Extract the grid spacing parallel to the surface
       diel%surfaces(ix)%dr = pack(tree%boxes(id_out)%dr, &
            [(i, i=1,NDIM)] /= af_neighb_dim(direction))

       associate (sd => diel%surfaces(ix)%sd, &
            sd_p => diel%surfaces(p_ix)%sd)
#if NDIM == 2
         ! Copy the values from the parent
         sd(1:nc:2, :) = sd_p(dix(1)+1:dix(1)+nc/2, :)
         sd(2:nc:2, :) = sd(1:nc:2, :)
#elif NDIM == 3
        ! Copy the values from the parent
!FLAG
         sd(1:nc:2, 1:nc:2, :) = sd_p(dix(1)+1:dix(1)+nc/2, dix(2)+1:dix(2)+nc/2, :)
         sd(2:nc:2, 1:nc:2, :) = sd(1:nc:2, 1:nc:2, :)
         sd(:, 2:nc:2, :)      = sd(:, 1:nc:2, :)
#endif
       end associate
    end do

    diel%surfaces(p_ix)%in_use = .false.

  end subroutine prolong_surface_from_parent

  !> Restrict a child surface to its parent
  subroutine restrict_surface_to_parent(tree, diel, ix)
    type(af_t), intent(in)            :: tree
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)               :: ix
    integer                           :: p_ix, nc, dix(NDIM-1), id_out, id_in

    p_ix = diel%surfaces(ix)%ix_parent
    if (p_ix == surface_none) error stop "Too much derefinement on surface"
    dix  = diel%surfaces(ix)%offset_parent
    nc   = tree%n_cell

    associate (sd => diel%surfaces(ix)%sd, &
         sd_p => diel%surfaces(p_ix)%sd)
#if NDIM == 2
      ! Average the value on the children
      sd_p(dix(1)+1:dix(1)+nc/2, :) = 0.5_dp * (sd(1:nc:2, :) + sd(2:nc:2, :))
#elif NDIM == 3
!FLAG
      sd_p(dix(1)+1:dix(1)+nc/2, dix(2)+1:dix(2)+nc/2, :) = 0.25_dp * &
      (sd(1:nc:2, 1:nc:2, :) + sd(2:nc:2, 1:nc:2, :) + sd(1:nc:2, 2:nc:2, :) + sd(2:nc:2, 1:nc:2, :) + sd(2:nc:2, 2:nc:2, :))
#endif
    end associate

    !$omp critical
    ! Remove child surface
    diel%n_removed = diel%n_removed + 1
    diel%removed_surfaces(diel%n_removed) = ix
    diel%surfaces(p_ix)%in_use = .true.
    !$omp end critical
    id_out = diel%surfaces(ix)%id_out
    id_in = diel%surfaces(ix)%id_in
    diel%box_id_out_to_surface_ix(id_out) = surface_none
    diel%box_id_in_to_surface_ix(id_in) = surface_none
    diel%surfaces(ix)%in_use = .false.

  end subroutine restrict_surface_to_parent

  !> Get an array of pairs of boxes (their indices) across a surface. This can
  !> be used in af_adjust_refinement to prevent refinement jumps over the
  !> surface.
  subroutine surface_get_refinement_links(diel, refinement_links)
    type(surfaces_t), intent(in)      :: diel
    !> Array of linked boxes, on output of size (2, n_links)
    integer, allocatable, intent(inout) :: refinement_links(:, :)
    integer                             :: max_ix, n, ix

    if (allocated(refinement_links)) deallocate(refinement_links)
    max_ix = diel%max_ix
    n      = count(diel%surfaces(1:max_ix)%in_use)
    allocate(refinement_links(2, n))

    n = 0
    do ix = 1, max_ix
       if (diel%surfaces(ix)%in_use) then
          n = n + 1
          refinement_links(:, n) = [diel%surfaces(ix)%id_out, &
               diel%surfaces(ix)%id_in]
       end if
    end do
  end subroutine surface_get_refinement_links

  !> Map surface charge to a cell-centered right-hand side
  subroutine surface_surface_charge_to_rhs(tree, diel, i_sigma, i_rhs, fac)
    type(af_t), intent(inout)      :: tree
    type(surfaces_t), intent(in) :: diel
    integer, intent(in)            :: i_sigma !< Surface charage variable
    integer, intent(in)            :: i_rhs   !< Rhs variable (in the tree)
    real(dp), intent(in)           :: fac     !< Multiplication factor
    integer                        :: n

    if (.not. diel%initialized) error stop "dielectric not initialized"

    do n = 1, diel%max_ix
       if (.not. diel%surfaces(n)%in_use) cycle

       call surface_charge_to_rhs(tree%boxes, diel%surfaces(n), &
            i_sigma, i_rhs, fac)
    end do
  end subroutine surface_surface_charge_to_rhs

  !> Routine that implements the mapping of surface charge to a cell-centered
  !> right-hand side
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
!FLAG
    boxes(id_out)%cc(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), i_rhs) = &
      boxes(id_out)%cc(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), i_rhs) &
      + frac_gas * fac * fac_to_volume * &
      reshape(surface%sd(:, :, i_sigma), [ghi(1)-glo(1)+1, ghi(2)-glo(2)+1, ghi(3)-glo(3)+1])

    boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), i_rhs) = &
      boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), i_rhs) &
      + (1-frac_gas) * fac * fac_to_volume * &
      reshape(surface%sd(:, :, i_sigma), [dhi(1)-dlo(1)+1, dhi(2)-dlo(2)+1, dhi(3)-dlo(3)+1])
#endif

  end subroutine surface_charge_to_rhs

  !> Map first layer inside the dielectric to surface density
  subroutine surface_inside_layer_to_surface(tree, diel, i_cc, i_sigma, &
       fac, clear_cc, clear_surf)
    type(af_t), intent(inout)         :: tree
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)               :: i_cc       !< Cell-centered variable (in the tree)
    integer, intent(in)               :: i_sigma    !< Surface charage variable
    real(dp), intent(in)              :: fac        !< Multiplication factor
    logical, intent(in)               :: clear_cc   !< Set density to zero afterwards
    logical, intent(in)               :: clear_surf !< Set surface density to zero initially
    integer                           :: n

    if (.not. diel%initialized) error stop "dielectric not initialized"

    do n = 1, diel%max_ix
       if (.not. diel%surfaces(n)%in_use) cycle

       call inside_layer_to_surface(tree%boxes, diel%surfaces(n), &
            i_cc, i_sigma, fac, clear_cc, clear_surf)
    end do
  end subroutine surface_inside_layer_to_surface

  !> Map first layer inside the dielectric to surface density
  subroutine inside_layer_to_surface(boxes, surface, i_cc, i_sigma, fac, &
       clear_cc, clear_surf)
    type(box_t), intent(inout)     :: boxes(:)
    type(surface_t), intent(inout) :: surface
    integer, intent(in)            :: i_cc       !< Cell-centered variable (in the tree)
    integer, intent(in)            :: i_sigma    !< Surface charage variable
    real(dp), intent(in)           :: fac        !< Multiplication factor
    logical, intent(in)            :: clear_cc   !< Set density to zero afterwards
    logical, intent(in)            :: clear_surf !< Set surface density to zero initially
    integer                        :: nb, nc, id_in
    integer                        :: dlo(NDIM), dhi(NDIM)
    real(dp)                       :: dr

    nb     = surface%direction
    id_in  = surface%id_in
    nc     = boxes(id_in)%n_cell
    dr     = boxes(id_in)%dr(af_neighb_dim(nb))

    ! Get index range for cells adjacent to the boundary of the dielectric box
    call af_get_index_bc_inside(af_neighb_rev(nb), nc, 1, dlo, dhi)

#if NDIM == 2
    if (clear_surf) surface%sd(:, i_sigma) = 0.0_dp

    surface%sd(:, i_sigma) = surface%sd(:, i_sigma) + fac * dr * &
         pack(boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_cc), .true.)
    if (clear_cc) boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_cc) = 0.0_dp
#elif NDIM == 3
    if (clear_surf) surface%sd(:, :, i_sigma) = 0.0_dp
!FLAG
    surface%sd(:,:, i_sigma) = surface%sd(:, :, i_sigma) + fac * dr * &
         reshape(boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), i_cc), [nc, nc])
    if (clear_cc) boxes(id_in)%cc(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), i_cc) = 0.0_dp
#endif

  end subroutine inside_layer_to_surface

  !> Compute the electric field at face centers near surfaces
  subroutine surface_correct_field_fc(tree, diel, i_sigma, i_fld, i_phi, fac)
    type(af_t), intent(inout)      :: tree
    type(surfaces_t), intent(in) :: diel
    integer, intent(in)            :: i_sigma !< Surface charge variable
    integer, intent(in)            :: i_fld   !< Face-centered field variable
    integer, intent(in)            :: i_phi   !< Cell-centered potential variable
    real(dp), intent(in)           :: fac     !< Elementary charge over eps0
    integer                        :: id_in, id_out, nc, nb, ix
    real(dp)                       :: eps, fac_fld(2), fac_charge, inv_dr(NDIM)
    integer                        :: dlo(NDIM), dhi(NDIM)
    integer                        :: glo(NDIM), ghi(NDIM)
    integer                        :: dlo_shft(NDIM), dhi_shft(NDIM)
    integer                        :: glo_shft(NDIM), ghi_shft(NDIM)
    integer                        :: shft(NDIM, 2*NDIM)

    if (.not. diel%initialized) error stop "dielectric not initialized"

    nc = tree%n_cell

    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
          nb = diel%surfaces(ix)%direction
          eps = diel%surfaces(ix)%eps
          id_out = diel%surfaces(ix)%id_out
          id_in = diel%surfaces(ix)%id_in
          inv_dr = 1/tree%boxes(id_out)%dr

          fac_fld = [2 * eps, 2.0_dp] / (1 + eps)
          fac_charge = fac / (1 + eps)

          ! Get index range for cells adjacent to the boundary of the dielectric box
          call af_get_index_bface_inside(af_neighb_rev(nb), nc, 1, dlo, dhi)

          ! Get similar index range for gas-phase box
          call af_get_index_bface_inside(nb, nc, 1, glo, ghi)

          ! Compute the offset for finite difference
          shft = - abs(af_neighb_dix(:, :))

          dlo_shft(:) = dlo(:) + shft(:, nb)
          dhi_shft(:) = dhi(:) + shft(:, nb)
          glo_shft(:) = glo(:) + shft(:, nb)
          ghi_shft(:) = ghi(:) + shft(:, nb)

          associate (fc_out => tree%boxes(id_out)%fc, &
               fc_in => tree%boxes(id_in)%fc, &
               cc_out => tree%boxes(id_out)%cc, &
               cc_in => tree%boxes(id_in)%cc, &
               sd => diel%surfaces(ix)%sd)
            select case (nb)
#if NDIM == 2
            case (af_neighb_lowx)
               fc_out(1, 1:nc, 1, i_fld) = fac_fld(1) * inv_dr(1) * &
                    (cc_out(0, 1:nc, i_phi) - cc_out(1, 1:nc, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
               fc_in(nc+1, 1:nc, 1, i_fld)  = fac_fld(2) * inv_dr(1) * &
                    (cc_in(nc, 1:nc, i_phi) - cc_in(nc+1, 1:nc, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
            case (af_neighb_highx)
               fc_out(nc+1, 1:nc, 1, i_fld)  = fac_fld(1) * inv_dr(1) * &
                    (cc_out(nc, 1:nc, i_phi) - cc_out(nc+1, 1:nc, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
               fc_in(1, 1:nc, 1, i_fld) = fac_fld(2) * inv_dr(1) * &
                    (cc_in(0, 1:nc, i_phi) - cc_in(1, 1:nc, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
            case (af_neighb_lowy)
               fc_out(1:nc, 1, 2, i_fld)  = fac_fld(1) * inv_dr(2) * &
                    (cc_out(1:nc, 0, i_phi) - cc_out(1:nc, 1, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
               fc_in(1:nc, nc+1, 2, i_fld) = fac_fld(2) * inv_dr(2) * &
                    (cc_in(1:nc, nc, i_phi) - cc_in(1:nc, nc+1, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
            case (af_neighb_highy)
               fc_out(1:nc, nc+1, 2, i_fld) = fac_fld(1) * inv_dr(2) * &
                    (cc_out(1:nc, nc, i_phi) - cc_out(1:nc, nc+1, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
               fc_in(1:nc, 1, 2, i_fld)  = fac_fld(2) * inv_dr(2) * &
                    (cc_in(1:nc, 0, i_phi) - cc_in(1:nc, 1, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
#elif NDIM == 3
            case default
                fc_out(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), af_neighb_dim(nb), i_fld) = &
                  fac_fld(1) * inv_dr(af_neighb_dim(nb)) * &
                  (cc_out(glo_shft(1):ghi_shft(1), glo_shft(2):ghi_shft(2), glo_shft(3):ghi_shft(3), i_phi) - &
                   cc_out(glo(1):ghi(1), glo(2):ghi(2), glo(3):ghi(3), i_phi)) &
                   + fac_charge * reshape(sd(:, :, i_sigma), [ghi(1)-glo(1)+1, ghi(2)-glo(2)+1, ghi(3)-glo(3)+1])

                fc_in(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), af_neighb_dim(nb), i_fld) = &
                  fac_fld(2) * inv_dr(af_neighb_dim(nb)) * &
                  (cc_in(dlo_shft(1):dhi_shft(1), dlo_shft(2):dhi_shft(2), dlo_shft(3):dhi_shft(3), i_phi) - &
                  cc_in(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3), i_phi)) &
                   - fac_charge * reshape(sd(:, :, i_sigma), [dhi(1)-dlo(1)+1, dhi(2)-dlo(2)+1, dhi(3)-dlo(3)+1])
#endif
            end select
          end associate
       end if
    end do

  end subroutine surface_correct_field_fc

#if NDIM == 2
  !> Compute the electric field at cell centers near surfaces
  subroutine surface_correct_field_cc(tree, diel, i_sigma, i_fld, i_phi, fac)
    type(af_t), intent(inout)      :: tree
    type(surfaces_t), intent(in) :: diel
    integer, intent(in)            :: i_sigma     !< Surface charge variable
    integer, intent(in)            :: i_fld(NDIM) !< Cell-centered field variables
    integer, intent(in)            :: i_phi       !< Cell-centered potential variable
    real(dp), intent(in)           :: fac         !< Elementary charge over eps0
    integer                        :: id_in, id_out, nc, nb, ix
    real(dp)                       :: eps, fac_fld(2), fac_charge, inv_dr(NDIM)
    real(dp)                       :: f_out(tree%n_cell), f_in(tree%n_cell)

    if (.not. diel%initialized) error stop "dielectric not initialized"

    nc = tree%n_cell

    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
          nb = diel%surfaces(ix)%direction
          eps = diel%surfaces(ix)%eps
          id_out = diel%surfaces(ix)%id_out
          id_in = diel%surfaces(ix)%id_in
          inv_dr = 1/tree%boxes(id_out)%dr

          fac_fld = [2 * eps, 2.0_dp] / (1 + eps)
          fac_charge = fac / (1 + eps)

          associate (cc_out => tree%boxes(id_out)%cc, &
               cc_in => tree%boxes(id_in)%cc, &
               sd => diel%surfaces(ix)%sd)

            ! Compute field at two cell faces and average them
            select case (nb)
            case (af_neighb_lowx)
               f_out = fac_fld(1) * inv_dr(1) * &
                    (cc_out(0, 1:nc, i_phi) - cc_out(1, 1:nc, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
               f_in = fac_fld(2) * inv_dr(1) * &
                    (cc_in(nc, 1:nc, i_phi) - cc_in(nc+1, 1:nc, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
               cc_out(1, 1:nc, i_fld(1)) = 0.5_dp * (f_out + inv_dr(1) * &
                    (cc_out(1, 1:nc, i_phi) - cc_out(2, 1:nc, i_phi)))
               cc_out(0, 1:nc, i_fld(1)) = f_out
               cc_in(nc, 1:nc, i_fld(1))  = 0.5_dp * (f_in + inv_dr(1) * &
                    (cc_in(nc-1, 1:nc, i_phi) - cc_in(nc, 1:nc, i_phi)))
               cc_in(nc+1, 1:nc, i_fld(1)) = f_in
            case (af_neighb_highx)
               f_out = fac_fld(1) * inv_dr(1) * &
                    (cc_out(nc, 1:nc, i_phi) - cc_out(nc+1, 1:nc, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
               f_in = fac_fld(2) * inv_dr(1) * &
                    (cc_in(0, 1:nc, i_phi) - cc_in(1, 1:nc, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
               cc_out(nc, 1:nc, i_fld(1))  = 0.5_dp * (f_out + inv_dr(1) * &
                    (cc_out(nc-1, 1:nc, i_phi) - cc_out(nc, 1:nc, i_phi)))
               cc_out(nc+1, 1:nc, i_fld(1)) = f_out
               cc_in(1, 1:nc, i_fld(1)) = 0.5_dp * (f_in + inv_dr(1) * &
                    (cc_in(1, 1:nc, i_phi) - cc_in(2, 1:nc, i_phi)))
               cc_in(0, 1:nc, i_fld(1)) = f_in
            case (af_neighb_lowy)
               f_out = fac_fld(1) * inv_dr(2) * &
                    (cc_out(1:nc, 0, i_phi) - cc_out(1:nc, 1, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
               f_in = fac_fld(2) * inv_dr(2) * &
                    (cc_in(1:nc, nc, i_phi) - cc_in(1:nc, nc+1, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
               cc_out(1:nc, 1, i_fld(2))  = 0.5_dp * (f_out + inv_dr(2) * &
                    (cc_out(1:nc, 1, i_phi) - cc_out(1:nc, 2, i_phi)))
               cc_out(1:nc, 0, i_fld(2)) = f_out
               cc_in(1:nc, nc, i_fld(2)) = 0.5_dp * (f_in + inv_dr(2) * &
                    (cc_in(1:nc, nc-1, i_phi) - cc_in(1:nc, nc, i_phi)))
               cc_in(1:nc, nc+1, i_fld(2)) = f_in
            case (af_neighb_highy)
               f_out = fac_fld(1) * inv_dr(2) * &
                    (cc_out(1:nc, nc, i_phi) - cc_out(1:nc, nc+1, i_phi)) &
                    - fac_charge * sd(:, i_sigma)
               f_in = fac_fld(2) * inv_dr(2) * &
                    (cc_in(1:nc, 0, i_phi) - cc_in(1:nc, 1, i_phi)) &
                    + fac_charge * sd(:, i_sigma)
               cc_out(1:nc, nc, i_fld(2)) = 0.5_dp * (f_out + inv_dr(2) * &
                    (cc_out(1:nc, nc-1, i_phi) - cc_out(1:nc, nc, i_phi)))
               cc_out(1:nc, nc+1, i_fld(2)) = f_out
               cc_in(1:nc, 1, i_fld(2))  = 0.5_dp * (f_in + inv_dr(2) * &
                    (cc_in(1:nc, 1, i_phi) - cc_in(1:nc, 2, i_phi)))
               cc_in(1:nc, 0, i_fld(2)) = f_in
            end select
          end associate
       end if
    end do
  end subroutine surface_correct_field_cc
#endif

  subroutine surface_get_surface_cell(tree, diel, x, ix_surf, ix_cell)
    use m_af_utils
    type(af_t), intent(in)         :: tree
    type(surfaces_t), intent(in) :: diel
    real(dp), intent(in)           :: x(NDIM)         !< Coordinate inside dielectric
    integer, intent(out)           :: ix_surf         !< Index of surface
    integer, intent(out)           :: ix_cell(NDIM-1) !< Index of cell on surface
    type(af_loc_t)                 :: loc
    integer                        :: i, id, dim, direction

    ! Find location in box
    loc = af_get_loc(tree, x)

    ! Check whether id is valid
    id = loc%id
    if (id == -1) error stop "Coordinate out of domain"

    ! Check for a surface from both sides
    ix_surf = diel%box_id_out_to_surface_ix(id)
    if (ix_surf == surface_none) ix_surf = diel%box_id_in_to_surface_ix(id)
    if (ix_surf == surface_none) error stop "No surface found"

    ! Extract the cell index but only parallel to the surface
    direction = diel%surfaces(ix_surf)%direction
    dim       = af_neighb_dim(direction)
    ix_cell   = pack(loc%ix, [(i, i=1,NDIM)] /= dim)
  end subroutine surface_get_surface_cell

    !> Write surface quantities to a separate output file
  subroutine surface_write_output(tree, diel, i_vars, var_names, output_name, &
       output_cnt)
    use m_npy
    type(af_t), intent(inout)    :: tree
    type(surfaces_t), intent(in) :: diel
    integer, intent(in)          :: i_vars(:) !< Indices of the variables
    character(len=*), intent(in) :: var_names(:) !< Names of the variables
    character(len=*), intent(in) :: output_name !< Base name for output
    integer, intent(in)          :: output_cnt !< Index of output
    integer                      :: n, i, ix, nc, n_vars
    integer                      :: lo(NDIM-1), hi(NDIM-1)
    character(len=200)           :: tmpname, filename

    real(dp)              :: coords(NDIM, tree%n_cell**(NDIM-1))
    real(dp), allocatable :: r(:, :), dr(:, :)
    real(dp), allocatable :: surface_vars(:, :)
    integer, allocatable  :: surf_dim(:)

    nc = diel%n_cell
    n = count(diel%surfaces(1:diel%max_ix)%in_use)
    n_vars = size(i_vars)
    allocate(r(NDIM, n*nc), surface_vars(n*nc, n_vars))
    allocate(dr(NDIM-1, n), surf_dim(n))

    write(filename, "(A,I6.6,A)") trim(output_name) // "_", &
         output_cnt, "_surface.npz"
    write(tmpname, "(A,I6.6,A)") trim(output_name) // "_", &
         output_cnt, "_tmp.npy"

    i = 0
    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
          i = i + 1
          lo = (i-1) * nc + 1
          hi = i * nc

          associate(box => tree%boxes(diel%surfaces(ix)%id_out), &
               surf => diel%surfaces(ix))
            call af_get_face_coords(box, surf%direction, coords)
#if NDIM == 2
            r(:, lo(1):hi(1)) = coords
            dr(:, i) = surf%dr
            surf_dim(i) = af_neighb_dim(surf%direction)
            surface_vars(lo(1):hi(1), :) = surf%sd(:, i_vars)
#elif NDIM == 3
            error stop "not implemented"
#endif
          end associate
       end if
    end do

    call save_npy(tmpname, [n])
    call add_to_zip(filename, tmpname, .false., "n_surfaces")
    call save_npy(tmpname, [nc])
    call add_to_zip(filename, tmpname, .false., "n_cell")
    call save_npy(tmpname, r)
    call add_to_zip(filename, tmpname, .false., "r")
    call save_npy(tmpname, dr)
    call add_to_zip(filename, tmpname, .false., "dr")
    call save_npy(tmpname, surf_dim)
    call add_to_zip(filename, tmpname, .false., "normal_dim")

    do n = 1, n_vars
       call save_npy(tmpname, surface_vars(:, n))
       call add_to_zip(filename, tmpname, .false., trim(var_names(n)))
    end do
    print *, "surface_write_output: written " // trim(filename)

  end subroutine surface_write_output

end module m_af_surface
