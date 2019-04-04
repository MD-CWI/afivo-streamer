#include "../afivo/src/cpp_macros.h"
!> Module with settings and routines to handle dielectrics
module m_dielectric
  use m_af_all

  implicit none
  private

  type surface_data_t
     integer :: id_gas
     integer :: id_diel
     integer :: direction
#if NDIM == 2
     real(dp), allocatable :: charge(:)
#elif NDIM == 3
     real(dp), allocatable :: charge(:, :)
#endif
  end type surface_data_t

  integer                           :: num_surfaces = 0
  type(surface_data_t), allocatable :: surface_list(:)
  integer                           :: n_removed_surfaces = 0
  integer, allocatable              :: removed_surfaces(:)
  integer, allocatable              :: box_id_to_surface_id(:)

  public :: dielectric_initialize
  public :: dielectric_get_surfaces
  public :: dielectric_update_after_refinement
  public :: dielectric_fix_refine

contains

  subroutine dielectric_initialize(tree, cfg)
    use m_config
    type(af_t), intent(in) :: tree
    type(CFG_t), intent(in) :: cfg

    allocate(box_id_to_surface_id(tree%box_limit))

    ! Maximum number of surfaces
    allocate(surface_list(tree%box_limit/2))
    allocate(removed_surfaces(tree%box_limit/2))

    box_id_to_surface_id(:) = -1
  end subroutine dielectric_initialize

  integer function get_new_surface()
    if (n_removed_surfaces > 0) then
       get_new_surface = removed_surfaces(n_removed_surfaces)
       n_removed_surfaces = n_removed_surfaces - 1
    else
       num_surfaces = num_surfaces + 1
       get_new_surface = num_surfaces
    end if
  end function get_new_surface

  subroutine dielectric_get_surfaces(tree)
    type(af_t), intent(in) :: tree
    integer                :: lvl, i, id, nb, nc, ix
    integer                :: id_diel

    nc = tree%n_cell

    ! Locate all boxes at the boundary of the dielectric
    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          call find_dielectric_neighbor(tree, id, id_diel, nb)

          if (id_diel > 0) then
             ix = get_new_surface()
             surface_list(ix)%id_gas = id
             surface_list(ix)%id_diel = id_diel
             surface_list(ix)%direction = nb

             if (.not. allocated(surface_list(ix)%charge)) then
#if NDIM == 2
                allocate(surface_list(ix)%charge(nc))
#elif NDIM == 3
                allocate(surface_list(ix)%charge(nc, nc))
#endif
             end if

          end if
       end do
    end do

  end subroutine dielectric_get_surfaces

  subroutine find_dielectric_neighbor(tree, id, id_diel, nb)
    use m_streamer
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: id
    integer, intent(out)   :: id_diel
    integer, intent(out)   :: nb
    integer                :: nb_id
    real(dp)               :: nb_eps

    id_diel = -1

    if (eps_changes_in_box(tree%boxes(id)) .and. &
         tree%boxes(id)%cc(DTIMES(1), i_eps) <= 1.0_dp) then

       do nb = 1, af_num_neighbors
          nb_id = tree%boxes(id)%neighbors(nb)
          if (nb_id <= af_no_box) cycle

          nb_eps = tree%boxes(nb_id)%cc(DTIMES(1), i_eps)
          if (nb_eps > 1.0_dp) then
             id_diel = nb_id
             exit
          end if
       end do
    end if

  end subroutine find_dielectric_neighbor

  logical function eps_changes_in_box(box)
    use m_streamer
    type(box_t), intent(in) :: box
    real(dp)                :: a, b
#if NDIM == 2
    a = minval(box%cc(:, :, i_eps))
    b = maxval(box%cc(:, :, i_eps))
#elif NDIM == 3
    a = minval(box%cc(:, :, :, i_eps))
    b = maxval(box%cc(:, :, :, i_eps))
#endif

    eps_changes_in_box = (b > a)
  end function eps_changes_in_box

  subroutine dielectric_update_after_refinement(tree, ref_info)
    type(af_t), intent(in)       :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, id_diel, nb, ix
    integer                      :: nc

    nc = tree%n_cell

    ! Handle removed surfaces
    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%rm)
          id = ref_info%lvls(lvl)%rm(i)

          ix = box_id_to_surface_id(id)
          if (ix > 0) then
             n_removed_surfaces = n_removed_surfaces + 1
             removed_surfaces(n_removed_surfaces) = ix
             box_id_to_surface_id(surface_list%id_gas) = -1
             box_id_to_surface_id(surface_list%id_diel) = -1

             call restrict_surface_to_parent(tree, surface_list(ix))
          end if
       end do
    end do

    ! Add new surfaces
    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)

          call find_dielectric_neighbor(tree, id, id_diel, nb)
          if (id_diel > 0) then
             ix = get_new_surface()
             surface_list(ix)%id_gas = id
             surface_list(ix)%id_diel = id_diel
             surface_list(ix)%direction = nb

             box_id_to_surface_id(id) = ix
             box_id_to_surface_id(id_diel) = ix

             if (.not. allocated(surface_list(ix)%charge)) then
#if NDIM == 2
                allocate(surface_list(ix)%charge(nc))
#elif NDIM == 3
                allocate(surface_list(ix)%charge(nc, nc))
#endif
             end if

             call prolong_surface_from_parent(tree, surface_list(ix))
          end if
       end do
    end do

  end subroutine dielectric_update_after_refinement

  subroutine prolong_surface_from_parent(tree, surface)
    type(af_t), intent(in) :: tree
    type(surface_data_t), intent(inout) :: surface

    ! error stop "todo"
  end subroutine prolong_surface_from_parent

  subroutine restrict_surface_to_parent(tree, surface)
    type(af_t), intent(in)           :: tree
    type(surface_data_t), intent(in) :: surface

    ! error stop "todo"
  end subroutine restrict_surface_to_parent

  subroutine dielectric_fix_refine(tree, ref_flags)
    use m_streamer
    type(af_t), intent(in) :: tree
    integer, intent(inout) :: ref_flags(:)
    integer                :: lvl, i, id, nb_id, ix

    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          ix = box_id_to_surface_id(id)
          if (ix > 0) then
             if (id == surface_list(ix)%id_gas) then
                nb_id = surface_list(ix)%id_diel
                ref_flags([id, nb_id]) = maxval(ref_flags([id, nb_id]))
             end if
          end if
       end do
    end do
  end subroutine dielectric_fix_refine

end module m_dielectric
