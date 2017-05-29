program test_photoionization
  use m_a2_types
  use m_a2_core
  use m_a2_output
  use m_a2_utils
  use m_lookup_table
  use m_random
  use m_photons

  implicit none

  integer, parameter  :: box_size     = 8
  integer, parameter :: i8 = selected_int_kind(18)
  real(dp), parameter :: frac_O2      = 0.2_dp
  real(dp), parameter :: domain_len   = 8e-3_dp
  real(dp), parameter :: dr           = domain_len / box_size

  integer, parameter  :: i_src        = 1
  integer, parameter  :: i_photo      = 2
  integer, parameter  :: i_sol        = 3

  real(dp) :: gas_pressure
  real(dp) :: grid_factor
  logical  :: use_const_dx
  logical  :: use_cyl
  integer(i8) :: rng_seed(2)  = [189023890123_i8, 348978834_i8]
  integer  :: num_photons

  type(a2_t)          :: tree
  type(ref_info_t)    :: ref_info
  type(photoi_tbl_t)  :: photoi_tbl
  type(RNG_t)         :: sim_rng       ! Random number generator

  integer             :: n, id
  integer             :: ix_list(2, 1) ! Spatial indices of initial boxes
  real(dp)            :: sum_pho

  character(len=100) :: fname, tmp

  if (command_argument_count() < 4) then
     print *, "Need >5 arguments: num_photons pressure grid_fac", &
          " const_dx cylindrical [seed]"
     stop
  end if

  call get_command_argument(1, tmp)
  read(tmp, *) num_photons
  call get_command_argument(2, tmp)
  read(tmp, *) gas_pressure
  call get_command_argument(3, tmp)
  read(tmp, *) grid_factor
  call get_command_argument(4, tmp)
  read(tmp, *) use_const_dx
  call get_command_argument(5, tmp)
  read(tmp, *) use_cyl

  if (command_argument_count() > 5) then
     call get_command_argument(6, tmp)
     read(tmp, *) rng_seed(1)
  end if

  print *, "Gas pressure (bar):", gas_pressure
  print *, "Fraction oxygen:   ", frac_O2

  call photoi_get_table_air(photoi_tbl, frac_O2 * gas_pressure, 2 * domain_len)

  ! Initialize tree
  if (use_cyl) then
     call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, &
          dr=dr, n_boxes=25*1000, coord=af_cyl, cc_names=["src", "pho", "sol"])
  else
     call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, &
          dr=dr, n_boxes=25*1000, cc_names=["src", "pho", "sol"])
  end if

  ! Set up geometry
  id             = 1          ! One box ...
  ix_list(:, id) = [1,1]      ! With index 1,1 ...

  ! Create the base mesh
  call a2_set_base(tree, 1, ix_list)

  do n = 1, 20
     call a2_adjust_refinement(tree, refine_routine, ref_info, 0)
     if (ref_info%n_add == 0) exit
  end do

  call sim_rng%set_seed(rng_seed)
  call set_photoionization(tree, num_photons)
  write(fname, "(A,I0,A,L1,L1,I0,A)") "pho_", nint(1e3_dp * gas_pressure), "_", &
       use_const_dx, use_cyl, nint(100 * grid_factor), ".silo"
  call a2_write_silo(tree, fname)
  call a2_tree_sum_cc(tree, i_photo, sum_pho)
  print *, "Sum photoionization", sum_pho
  call a2_tree_sum_cc(tree, i_sol, sum_pho)
  print *, "Sum solution", sum_pho

contains

  ! Refinement function
  subroutine refine_routine(box, cell_flags)
    type(box2_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)

    if (box%dr > 1.0e-3_dp * domain_len) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  subroutine set_photoionization(tree, num_photons)
    use m_photons
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: num_photons

    call a2_loop_box_arg(tree, set_photoionization_rate, [1.0_dp], .true.)
    call photoi_set_src_2d(tree, photoi_tbl, sim_rng, num_photons, &
         i_src, i_photo, grid_factor, use_const_dx, use_cyl, 0.05e-3_dp)

  end subroutine set_photoionization

  subroutine set_photoionization_rate(box, coeff)
    use m_geometry
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
                  photoi_absorption_func_air(r, gas_pressure * frac_O2)

             if (r < box%dr) then
                box%cc(i, j, i_src) = 0.5_dp * coeff(1) / &
                     (2 * pi * (0.5_dp * box%dr) * box%dr**2)
             else
                box%cc(i, j, i_src) = 0.0_dp
             end if
          else
             box%cc(i, j, i_sol) = coeff(1) / (2 * pi * r) * &
                  photoi_absorption_func_air(r, gas_pressure * frac_O2)

             if (r < box%dr) then
                box%cc(i, j, i_src) = 0.5_dp * coeff(1) / box%dr**2
             else
                box%cc(i, j, i_src) = 0.0_dp
             end if
          end if
       end do
    end do
  end subroutine set_photoionization_rate

end program test_photoionization
