program test_photoionization
  use m_afivo_2d
  use m_lookup_table
  use m_random
  use m_photons

  implicit none

  integer, parameter  :: dp           = kind(0.0d0)
  integer, parameter  :: box_size     = 8
  real(dp), parameter :: frac_O2      = 0.2_dp
  real(dp), parameter :: gas_pressure = 1.0e0_dp
  real(dp), parameter :: domain_len   = 8e-3_dp
  real(dp), parameter :: dr           = domain_len / box_size
  integer, parameter  :: num_photons  = 10*1000

  integer, parameter  :: i_src        = 1
  integer, parameter  :: i_pho        = 2
  integer, parameter  :: i_sol        = 3

  real(dp), parameter :: grid_factor  = 0.75_dp
  logical, parameter  :: use_const_dx = .false.

  type(a2_t)          :: tree
  type(ref_info_t)    :: ref_info
  type(PH_tbl_t)      :: photoi_tbl
  type(RNG_t)         :: sim_rng       ! Random number generator

  integer             :: n, id
  integer             :: ix_list(2, 1) ! Spatial indices of initial boxes
  integer             :: nb_list(4, 1) ! Neighbors of initial boxes
  real(dp)            :: sum_pho

  character(len=100) :: fname
  character(len=10)  :: cc_names(3) = &
       [character(len=10) :: "src", "pho", "sol"]

  call PH_get_tbl_air(photoi_tbl, frac_O2 * gas_pressure, 2 * domain_len)

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, &
       dr=dr, n_boxes = 25*1000)

  ! Set up geometry
  id             = 1          ! One box ...
  ix_list(:, id) = [1,1]      ! With index 1,1 ...
  nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)

  do n = 1, 20
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  call set_photoionization(tree, num_photons)
  write(fname, "(A,I0,A,L1,I0,A)") "pho_", nint(1e6_dp * gas_pressure), "_", &
       use_const_dx, nint(100 * grid_factor), ".silo"
  call a2_write_silo(tree, fname, &
       cc_names, 0, 0.0_dp)
  call a2_tree_sum_cc(tree, i_pho, sum_pho, .true.)
  print *, "Sum photoionization", sum_pho
  call a2_tree_sum_cc(tree, i_sol, sum_pho, .true.)
  print *, "Sum solution", sum_pho

contains

  ! Refinement function
  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)

    if (boxes(id)%dr > 1.0e-3_dp * domain_len) then
       ref_flags(id) = a5_do_ref
    end if

  end subroutine set_ref_flags

  subroutine set_photoionization(tree, num_photons)
    use m_photons
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: num_photons

    call a2_loop_box_arg(tree, set_photoi_rate, [1.0_dp], .true.)
    call PH_set_src_2d(tree, photoi_tbl, sim_rng, num_photons, &
         i_src, i_pho, grid_factor, use_const_dx)

  end subroutine set_photoionization

  subroutine set_photoi_rate(box, coeff)
    use m_geom
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), xy_rel(2), r
    real(dp), parameter         :: pi = acos(-1.0_dp)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          xy_rel = xy - 0.5_dp * domain_len
          r = norm2(xy_rel)
          if (r < box%dr) then
             box%cc(i, j, i_src) = 0.25_dp * coeff(1) / box%dr**2
          else
             box%cc(i, j, i_src) = 0.0_dp
          end if
          box%cc(i, j, i_sol) = coeff(1) / (2 * pi * r) * &
               PH_absfunc_air(r, gas_pressure * frac_O2)
       end do
    end do
  end subroutine set_photoi_rate

end program test_photoionization
