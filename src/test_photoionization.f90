program test_photoionization
  use m_a2_t
  use m_a2_core
  use m_a2_io
  use m_a2_utils
  use m_lookup_table
  use m_random
  use m_photons

  implicit none

  integer, parameter  :: box_size     = 8
  real(dp), parameter :: frac_O2      = 0.2_dp
  real(dp), parameter :: domain_len   = 8e-3_dp
  real(dp), parameter :: dr           = domain_len / box_size
  integer, parameter  :: num_photons  = 50*1000

  integer, parameter  :: i_src        = 1
  integer, parameter  :: i_pho        = 2
  integer, parameter  :: i_sol        = 3

  real(dp) :: gas_pressure = 1.0e-3_dp
  real(dp) :: grid_factor  = 0.3_dp
  logical  :: use_const_dx = .false.
  logical  :: use_cyl      = .false.
  integer  :: rng_seed(4)  = [234, 45, 843, 234]

  type(a2_t)          :: tree
  type(ref_info_t)    :: ref_info
  type(PH_tbl_t)      :: photoi_tbl
  type(RNG_t)         :: sim_rng       ! Random number generator

  integer             :: n, id
  integer             :: ix_list(2, 1) ! Spatial indices of initial boxes
  integer             :: nb_list(4, 1) ! Neighbors of initial boxes
  real(dp)            :: sum_pho

  character(len=100) :: fname, tmp
  character(len=10)  :: cc_names(3) = &
       [character(len=10) :: "src", "pho", "sol"]

  if (command_argument_count() < 4) &
       stop "Need 4 arguments: pressure grid_fac const_dx cylindrical"
  call get_command_argument(1, tmp)
  read(tmp, *) gas_pressure
  call get_command_argument(2, tmp)
  read(tmp, *) grid_factor
  call get_command_argument(3, tmp)
  read(tmp, *) use_const_dx
  call get_command_argument(4, tmp)
  read(tmp, *) use_cyl

  if (command_argument_count() > 4) then
     call get_command_argument(5, tmp)
     read(tmp, *) rng_seed(1)
  end if

  print *, "Gas pressure (bar):", gas_pressure
  print *, "Fraction oxygen:   ", frac_O2
  
  call PH_get_tbl_air(photoi_tbl, frac_O2 * gas_pressure, 2 * domain_len)

  ! Initialize tree
  if (use_cyl) then
     call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, &
          dr=dr, n_boxes=25*1000, coord=a5_cyl)
  else
     call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, &
          dr=dr, n_boxes=25*1000)
  end if

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

  call sim_rng%set_seed(rng_seed)
  call set_photoionization(tree, num_photons)
  write(fname, "(A,I0,A,L1,L1,I0,A)") "pho_", nint(1e3_dp * gas_pressure), "_", &
       use_const_dx, use_cyl, nint(100 * grid_factor), ".silo"
  call a2_write_silo(tree, fname, &
       cc_names, 0, 0.0_dp)
  call a2_tree_sum_cc(tree, i_pho, sum_pho)
  print *, "Sum photoionization", sum_pho
  call a2_tree_sum_cc(tree, i_sol, sum_pho)
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
         i_src, i_pho, grid_factor, use_const_dx, use_cyl, 0.05e-3_dp)

  end subroutine set_photoionization

  subroutine set_photoi_rate(box, coeff)
    use m_geom
    use m_a2_utils
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), xy_rel(2), r
    real(dp), parameter         :: pi = acos(-1.0_dp)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          xy_rel = xy - [0.0_dp, 0.5_dp * domain_len]
          r = norm2(xy_rel)

          if (use_cyl) then
             box%cc(i, j, i_sol) = coeff(1) / (4 * pi * r**2) * &
                  PH_absfunc_air(r, gas_pressure * frac_O2)

             if (r < box%dr) then
                box%cc(i, j, i_src) = 0.5_dp * coeff(1) / &
                     (2 * pi * (0.5_dp * box%dr) * box%dr**2)
             else
                box%cc(i, j, i_src) = 0.0_dp
             end if
          else
             box%cc(i, j, i_sol) = coeff(1) / (2 * pi * r) * &
                  PH_absfunc_air(r, gas_pressure * frac_O2)

             if (r < box%dr) then
                box%cc(i, j, i_src) = 0.5_dp * coeff(1) / box%dr**2
             else
                box%cc(i, j, i_src) = 0.0_dp
             end if
          end if
       end do
    end do
  end subroutine set_photoi_rate

end program test_photoionization
