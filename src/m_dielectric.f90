#include "../afivo/src/cpp_macros.h"
!> Module with settings and routines to handle dielectrics
module m_dielectric
  use m_af_all

  implicit none
  private

  type surface_data_t
     logical :: in_use = .false.
     integer :: id_gas
     integer :: id_diel
     integer :: direction
     real(dp) :: eps
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
  public :: dielectric_rearrange_charge
  public :: dielectric_adjust_field
  public :: dielectric_copy_fluxes

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
    real(dp)               :: eps

    nc = tree%n_cell

    ! Locate all boxes at the boundary of the dielectric
    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          call find_dielectric_neighbor(tree, id, id_diel, nb, eps)

          if (id_diel > 0) then
             ix = get_new_surface()
             surface_list(ix)%in_use = .true.
             surface_list(ix)%id_gas = id
             surface_list(ix)%id_diel = id_diel
             surface_list(ix)%direction = nb
             surface_list(ix)%eps = eps

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

  subroutine find_dielectric_neighbor(tree, id, id_diel, nb, nb_eps)
    use m_streamer
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: id
    integer, intent(out)   :: id_diel
    integer, intent(out)   :: nb
    real(dp), intent(out)  :: nb_eps
    integer                :: nb_id
    real(dp), parameter    :: eps_threshold = 1e-8_dp

    id_diel = -1
    nb      = -1
    nb_eps  = -1.0_dp

    if (eps_changes_in_box(tree%boxes(id)) .and. &
         tree%boxes(id)%cc(DTIMES(1), i_eps) < 1 + eps_threshold) then

       do nb = 1, af_num_neighbors
          nb_id = tree%boxes(id)%neighbors(nb)
          if (nb_id <= af_no_box) cycle

          nb_eps = tree%boxes(nb_id)%cc(DTIMES(1), i_eps)
          if (nb_eps > 1 + eps_threshold) then
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
    real(dp)                     :: eps

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
             surface_list(ix)%in_use = .false.

             call restrict_surface_to_parent(tree, surface_list(ix))
          end if
       end do
    end do

    ! Add new surfaces
    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)

          call find_dielectric_neighbor(tree, id, id_diel, nb, eps)

          if (id_diel > 0) then
             ix = get_new_surface()
             surface_list(ix)%in_use = .true.
             surface_list(ix)%id_gas = id
             surface_list(ix)%id_diel = id_diel
             surface_list(ix)%direction = nb
             surface_list(ix)%eps = eps

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

  subroutine dielectric_rearrange_charge(tree)
    type(af_t), intent(inout) :: tree
    integer                   :: n, id_gas

    do n = 1, num_surfaces
       if (.not. surface_list(n)%in_use) cycle
       id_gas = surface_list(n)%id_gas

       if (af_has_children(tree%boxes(id_gas))) cycle

       call rearrange_charge_on_surface(tree%boxes, &
            surface_list(n), tree%n_cell)
    end do
  end subroutine dielectric_rearrange_charge

  subroutine rearrange_charge_on_surface(boxes, surface, nc)
    use m_streamer
    use m_units_constants
    type(box_t), intent(inout)          :: boxes(:)
    type(surface_data_t), intent(inout) :: surface
    integer, intent(in)                 :: nc
    integer                             :: nb, id_gas, id_diel
    integer                             :: glo(NDIM), ghi(NDIM)
    integer                             :: dlo(NDIM), dhi(NDIM)
    real(dp)                            :: frac_gas, dr

    nb      = surface%direction
    id_gas  = surface%id_gas
    id_diel = surface%id_diel
    dr      = boxes(id_gas)%dr(af_neighb_dim(nb))

    ! Get index range for dielectric box
    call af_get_index_bc_inside(af_neighb_rev(nb), nc, dlo, dhi)

    ! Get index range for gas-phase box
    call af_get_index_bc_inside(nb, nc, glo, ghi)

    ! How much of the rhs to put on the gas and dielectric side
    frac_gas  = 1.0_dp / (1.0_dp + surface%eps)

    ! The rhs is defined as -rho/eps_0, hence the minus sign. Multiply with the
    ! grid spacing to transform a charge density to a surface charge.
    surface%charge = -reshape(boxes(id_diel)%cc(dlo(1):dhi(1), &
         dlo(2):dhi(2), i_rhs), shape(surface%charge)) * dr

#if NDIM == 2
    boxes(id_gas)%cc(glo(1):ghi(1), glo(2):ghi(2), i_rhs) = &
         boxes(id_gas)%cc(glo(1):ghi(1), glo(2):ghi(2), i_rhs) &
         - frac_gas / dr * &
         reshape(surface%charge, [ghi(1)-glo(1)+1, ghi(2)-glo(2)+1])

    boxes(id_diel)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_rhs) = &
         -(1-frac_gas) / dr * &
         reshape(surface%charge, [dhi(1)-dlo(1)+1, dhi(2)-dlo(2)+1])
#elif NDIM == 3
    error stop
#endif

  end subroutine rearrange_charge_on_surface

  subroutine dielectric_adjust_field(boxes, id)
    use m_streamer
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)        :: id
    integer                    :: nc, nb, ix
    real(dp)                   :: eps, fac_fld, fac_charge

    nc = boxes(id)%n_cell

    ix = box_id_to_surface_id(id)

    if (ix <= 0) return

    eps = surface_list(ix)%eps

    if (id == surface_list(ix)%id_gas) then
       nb = surface_list(ix)%direction
       fac_fld = 2 * eps / (1 + eps)
       fac_charge = 1 / (1 + eps)
    else
       nb = af_neighb_rev(surface_list(ix)%direction)
       fac_fld = 2 / (1 + eps)
       fac_charge = 1 / (1 + eps)
    end if

#if NDIM == 2
    ! Compute fields at the boundaries of the box, where eps can change
    select case (nb)
    case (af_neighb_lowx)
       boxes(id)%fc(1, 1:nc, 1, electric_fld) = fac_fld * &
            boxes(id)%fc(1, 1:nc, 1, electric_fld) + &
            fac_charge * surface_list(ix)%charge
    case (af_neighb_highx)
       boxes(id)%fc(nc+1, 1:nc, 1, electric_fld) = fac_fld * &
            boxes(id)%fc(nc+1, 1:nc, 1, electric_fld) - &
            fac_charge * surface_list(ix)%charge
    case (af_neighb_lowy)
       boxes(id)%fc(1:nc, 1, 2, electric_fld) = fac_fld * &
            boxes(id)%fc(1:nc, 1, 2, electric_fld) + &
            fac_charge * surface_list(ix)%charge
    case (af_neighb_highy)
       boxes(id)%fc(1:nc, nc+1, 2, electric_fld) = fac_fld * &
            boxes(id)%fc(1:nc, nc+1, 2, electric_fld) - &
            fac_charge * surface_list(ix)%charge
    end select
#elif NDIM == 3
    error stop
    ! box%fc(1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
    !      (box%cc(0, 1:nc, 1:nc, i_phi) - box%cc(1, 1:nc, 1:nc, i_phi)) * &
    !      box%cc(0, 1:nc, 1:nc, i_eps) / &
    !      (box%cc(1, 1:nc, 1:nc, i_eps) + box%cc(0, 1:nc, 1:nc, i_eps))
    ! box%fc(nc+1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
    !      (box%cc(nc, 1:nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, 1:nc, i_phi)) * &
    !      box%cc(nc+1, 1:nc, 1:nc, i_eps) / &
    !      (box%cc(nc+1, 1:nc, 1:nc, i_eps) + box%cc(nc, 1:nc, 1:nc, i_eps))
    ! box%fc(1:nc, 1, 1:nc, 2, electric_fld) = 2 * inv_dr(2) * &
    !      (box%cc(1:nc, 0, 1:nc, i_phi) - box%cc(1:nc, 1, 1:nc, i_phi)) * &
    !      box%cc(1:nc, 0, 1:nc, i_eps) / &
    !      (box%cc(1:nc, 1, 1:nc, i_eps) + box%cc(1:nc, 0, 1:nc, i_eps))
    ! box%fc(1:nc, nc+1, 1:nc, 2, electric_fld) = 2 * inv_dr(2) * &
    !      (box%cc(1:nc, nc, 1:nc, i_phi) - box%cc(1:nc, nc+1, 1:nc, i_phi)) * &
    !      box%cc(1:nc, nc+1, 1:nc, i_eps) / &
    !      (box%cc(1:nc, nc+1, 1:nc, i_eps) + box%cc(1:nc, nc, 1:nc, i_eps))
    ! box%fc(1:nc, 1:nc, 1, 3, electric_fld) = 2 * inv_dr(3) * &
    !      (box%cc(1:nc, 1:nc, 0, i_phi) - box%cc(1:nc, 1:nc, 1, i_phi)) * &
    !      box%cc(1:nc, 1:nc, 0, i_eps) / &
    !      (box%cc(1:nc, 1:nc, 1, i_eps) + box%cc(1:nc, 1:nc, 0, i_eps))
    ! box%fc(1:nc, 1:nc, nc+1, 3, electric_fld) = 2 * inv_dr(3) * &
    !      (box%cc(1:nc, 1:nc, nc, i_phi) - box%cc(1:nc, 1:nc, nc+1, i_phi)) * &
    !      box%cc(1:nc, 1:nc, nc+1, i_eps) / &
    !      (box%cc(1:nc, 1:nc, nc+1, i_eps) + box%cc(1:nc, 1:nc, nc, i_eps))
#endif

  end subroutine dielectric_adjust_field

  subroutine dielectric_copy_fluxes(tree, iflux)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: iflux
    integer                   :: n, id_diel, nb

    do n = 1, num_surfaces
       if (.not. surface_list(n)%in_use) cycle
       id_diel = surface_list(n)%id_diel

       if (af_has_children(tree%boxes(id_diel))) cycle

       nb = af_neighb_rev(surface_list(n)%direction)
       call copy_flux_from_neighbor(tree%boxes, id_diel, nb, iflux)
    end do
  end subroutine dielectric_copy_fluxes

  subroutine copy_flux_from_neighbor(boxes, id, nb, iflux)
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)        :: id
    integer, intent(in)        :: nb
    integer, intent(in)        :: iflux
    integer                    :: nc, nb_id

    nc    = boxes(id)%n_cell
    nb_id = boxes(id)%neighbors(nb)

    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
       boxes(id)%fc(1, 1:nc, 1, iflux) = &
            boxes(nb_id)%fc(nc+1, 1:nc, 1, iflux)
    case (af_neighb_highx)
       boxes(id)%fc(nc+1, 1:nc, 1, iflux) = &
            boxes(nb_id)%fc(1, 1:nc, 1, iflux)
    case (af_neighb_lowy)
       boxes(id)%fc(1:nc, 1, 2, iflux) = &
            boxes(nb_id)%fc(1:nc, nc+1, 2, iflux)
    case (af_neighb_highy)
       boxes(id)%fc(1:nc, nc+1, 2, iflux) = &
            boxes(nb_id)%fc(1:nc, 1, 2, iflux)
#elif NDIM == 3
       error stop
#endif

    end select
  end subroutine copy_flux_from_neighbor


end module m_dielectric
