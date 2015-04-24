program test_photoionization
  use m_afivo_2d
  use m_lookup_table
  use m_random
  use m_photons

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: box_size = 8
  real(dp), parameter :: domain_len = 4e-3_dp
  real(dp), parameter :: dr = domain_len / box_size
  real(dp), parameter :: eta = 0.05_dp
  real(dp), parameter :: frac_O2 = 0.2_dp
  real(dp), parameter :: gas_pressure = 1.0_dp
  integer, parameter :: num_photons = 100*1000

  integer, parameter :: i_src = 1
  integer, parameter :: i_pho = 2
  integer, parameter :: i_sol = 3

  type(a2_t) :: tree
  type(ref_info_t) :: ref_info
  type(LT_table_t) :: photoi_tbl
  type(RNG_t)        :: sim_rng ! Random number generator

  integer :: n, id
  integer    :: ix_list(2, 1) ! Spatial indices of initial boxes
  integer    :: nb_list(4, 1) ! Neighbors of initial boxes
  real(dp) :: sum_pho

  character(len=10)  :: cc_names(3) = &
       [character(len=10) :: "src", "pho", "sol"]

  photoi_tbl = PH_get_tbl_air(frac_O2 * gas_pressure)

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, dr=dr)

  ! Set up geometry
  id             = 1          ! One box ...
  ix_list(:, id) = [1,1]      ! With index 1,1 ...
  nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)

  do n = 1, 20
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
  end do

  call set_photoionization(tree, eta, gas_pressure, num_photons)
  call a2_write_silo(tree, "test_photoionization.silo", &
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

    if (boxes(id)%dr > 1.0e-5_dp) then
       ref_flags(id) = a5_do_ref
    end if

  end subroutine set_ref_flags

  subroutine set_photoionization(tree, eta, gas_pressure, num_photons)
    use m_photons
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta, gas_pressure
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30.0D0 * UC_torr_to_bar
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (gas_pressure + p_quench)
    print *, "Coeff. photoionization", eta * quench_fac

    call a2_loop_box_arg(tree, set_photoi_rate, [eta * quench_fac], .true.)
    call PH_set_src_2d(tree, photoi_tbl, sim_rng, 0.25_dp, &
         num_photons, i_src, i_pho)

  end subroutine set_photoionization

  subroutine set_photoi_rate(box, coeff)
    use m_geom
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), xy_rel(2), r

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          xy_rel = xy - 0.5_dp * domain_len
          r = norm2(xy_rel)
          if (r < box%dr) then
             box%cc(i, j, i_src) = 0.25_dp * coeff(1) / box%dr**2
             box%cc(i, j, i_sol) = 0.0_dp
          else
             box%cc(i, j, i_src) = 0.0_dp
             box%cc(i, j, i_sol) = coeff(1) * &
               PH_absfunc_air(r, frac_O2) / (2 * acos(-1.0_dp) * r)
          end if
       end do
    end do
  end subroutine set_photoi_rate

end program test_photoionization
