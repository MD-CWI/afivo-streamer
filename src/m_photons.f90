module m_photons

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: PH_get_tbl_air
  public :: PH_do_absorp
  public :: PH_absfunc_air
  public :: PH_set_src_2d
  public :: PH_set_src_3d

contains

  function PH_get_tbl_air(p_O2) result(tbl)
    use m_lookup_table
    real(dp), intent(in) :: p_O2     !< Partial pressure of oxygen (bar)
    type(LT_table_t)     :: tbl      !< The lookup table

    integer, parameter   :: tbl_size  = 1000
    integer              :: n
    real(dp), parameter  :: huge_dist = 1e9_dp, sixth = 1/6.0_dp
    real(dp), allocatable :: fsum(:), dist(:)
    real(dp)             :: dr, r_prev, incr

    allocate(fsum(tbl_size))
    allocate(dist(tbl_size))
    tbl     = LT_create(0.0_dp, 1.0_dp, tbl_size, 1)
    dr      = 1e-10_dp
    dist(1) = 0
    fsum(1) = 0
    r_prev = epsilon(1.0_dp)

    ! Use Simpson's rule for integration
    do n = 2, tbl_size
       incr = dr * sixth * (PH_absfunc_air(r_prev, p_O2) + &
            4 * PH_absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
            PH_absfunc_air(r_prev+dr, p_O2))

       do while (incr > 2.0_dp / tbl_size)
          dr = dr * 0.5_dp
          incr = dr * sixth * (PH_absfunc_air(r_prev, p_O2) + &
               4 * PH_absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
               PH_absfunc_air(r_prev+dr, p_O2))
       end do

       do while (incr < 1.0_dp / (tbl_size-1) .and. dr < huge_dist)
          dr = dr * 2
          incr = dr * sixth * (PH_absfunc_air(r_prev, p_O2) + &
               4 * PH_absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
               PH_absfunc_air(r_prev+dr, p_O2))
       end do

       dist(n) = r_prev + dr
       fsum(n) = fsum(n-1) + incr

       ! The integral will not be completely accurate, and small errors can lead
       ! to big problems (dr -> infinity). So, stop when dr gets too large.
       if ((n > 2 .and. dr > 10 * r_prev) .or. fsum(n) > 1)  then
          fsum(n:) = 1
          dist(n:) = dist(n-1)
          exit
       end if

       r_prev = r_prev + dr
    end do

    call LT_set_col(tbl, 1, fsum, dist)
  end function PH_get_tbl_air

  real(dp) function PH_absfunc_air(dist, p_O2)
    use m_units_constants
    real(dp), intent(in) :: dist, p_O2
    real(dp), parameter  :: chi_min  = 3.5_dp / UC_torr_to_bar
    real(dp), parameter  :: chi_max  = 200 / UC_torr_to_bar

    PH_absfunc_air = (exp(-chi_min * p_O2 * dist) - &
         exp(-chi_max * p_O2 * dist)) / (dist * log(chi_max/chi_min))
  end function PH_absfunc_air

  integer function get_lvl_length(dr_base, length)
    real(dp), intent(in) :: dr_base, length
    real(dp), parameter :: invlog2 = 1 / log(2.0_dp)
    real(dp) :: ratio

    ratio = dr_base / length
    if (ratio <= 1) then
       get_lvl_length = 1
    else
       get_lvl_length = 1 + ceiling(log(ratio) * invlog2)
    end if
  end function get_lvl_length

  integer function get_rlvl_length(dr_base, length, rng)
    use m_random
    real(dp), intent(in) :: dr_base, length
    real(dp), parameter :: invlog2 = 1 / log(2.0_dp)
    type(RNG_t), intent(inout) :: rng
    real(dp) :: ratio, tmp

    ratio = dr_base / length
    if (ratio <= 1) then
       get_rlvl_length = 1
    else
       tmp = log(ratio) * invlog2
       get_rlvl_length = floor(tmp)
       if (rng%uni_01() < tmp - get_rlvl_length) &
            get_rlvl_length = get_rlvl_length + 1
    end if
  end function get_rlvl_length

  subroutine PH_do_absorp(xyz_in, xyz_out, n_photons, tbl, rng)
    use m_lookup_table
    use m_random
    use omp_lib
    integer, intent(in)          :: n_photons
    real(dp), intent(in)         :: xyz_in(3, n_photons)
    real(dp), intent(out)        :: xyz_out(3, n_photons)
    type(LT_table_t), intent(in) :: tbl
    type(RNG_t), intent(inout)   :: rng
    integer                      :: n, n_procs, proc_id
    real(dp)                     :: rr, dist
    type(PRNG_t)                 :: prng

    !$omp parallel private(n, rr, dist, proc_id)
    !$omp single
    n_procs = omp_get_num_threads()
    call prng%init(n_procs, rng)
    !$omp end single

    proc_id = 1+omp_get_thread_num()

    !$omp do
    do n = 1, n_photons
       rr = prng%rngs(proc_id)%uni_01()
       dist = LT_get_col(tbl, 1, rr)
       xyz_out(:, n) =  xyz_in(:, n) + prng%rngs(proc_id)%sphere(dist)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine PH_do_absorp

  subroutine PH_set_src_2d(tree, pi_tbl, rng, rel_length, &
       num_photons, i_src, i_pho)
    use m_random
    use m_afivo_2d
    use m_lookup_table
    use omp_lib

    type(a2_t), intent(inout)   :: tree   !< Tree
    type(LT_table_t)            :: pi_tbl !< Table to sample abs. lenghts
    type(RNG_t), intent(inout)  :: rng    !< Random number generator
    !> Use dx smaller than pi_tbl at this value
    real(dp), intent(in) :: rel_length
    !> How many discrete photons to use
    integer, intent(in)         :: num_photons
    !> Input variable that contains photon production per cell
    integer, intent(in)         :: i_src
    !> Output variable that contains photoionization source rate
    integer, intent(in)         :: i_pho

    integer                     :: lvl, ix, id, nc
    integer                     :: i, j, n, n_create, n_used, i_ph
    integer                     :: proc_id, n_procs
    integer                     :: pho_lvl
    real(dp)                    :: r_create, dr, fac, dist
    real(dp)                    :: sum_production
    real(dp), allocatable       :: xyz_src(:, :)
    real(dp), allocatable       :: xyz_dst(:, :)
    type(PRNG_t)                :: prng
    type(a2_loc_t), allocatable :: ph_loc(:)

    nc = tree%n_cell

    ! Allocate a bit more space because of stochastic production
    allocate(xyz_src(3, nint(1.2_dp * num_photons + 1000)))

    ! Compute the sum of photon production
    call a2_tree_sum_cc(tree, i_src, sum_production, .true.)

    ! Create approximately num_photons
    fac    = num_photons / max(sum_production, epsilon(1.0_dp))
    n_used = 0
    ! print *, "num photons / sum", num_photons, sum_production

    ! Now loop over all leaves and create photons using random numbers

    !$omp parallel private(lvl, ix, id, i, j, dr, i_ph, proc_id, r_create, n_create)

    !$omp single
    n_procs = omp_get_num_threads()
    call prng%init(n_procs, rng)
    !$omp end single

    proc_id = 1+omp_get_thread_num()

    do lvl = 1, tree%max_lvl
       dr = a2_lvl_dr(tree, lvl)
       !$omp do
       do ix = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(ix)

          do j = 1, nc
             do i = 1, nc
                r_create = fac * tree%boxes(id)%cc(i, j, i_src) * dr**2
                n_create = floor(r_create)

                if (prng%rngs(proc_id)%uni_01() < r_create - n_create) &
                     n_create = n_create + 1

                if (n_create > 0) then
                   !$omp critical
                   i_ph = n_used
                   n_used = n_used + n_create
                   !$omp end critical

                   do n = 1, n_create
                      xyz_src(1:2, i_ph+n) = a2_r_cc(tree%boxes(id), [i, j])
                      xyz_src(3, i_ph+n) = 0
                   end do
                end if
             end do
          end do
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    allocate(xyz_dst(3, n_used))
    allocate(ph_loc(n_used))

    ! Get location of absorbption
    call PH_do_absorp(xyz_src, xyz_dst, n_used, pi_tbl, rng)

    ! TODO
    !x$omp parallel do private(dist, lvl)
    do n = 1, n_used
       dist = norm2(xyz_dst(1:2, n) - xyz_src(1:2, n))
       lvl = get_rlvl_length(tree%dr_base, 0.6_dp * dist, rng)
       ph_loc(n) = a2_get_loc(tree, xyz_dst(1:2, n), lvl)
    end do
    !x$omp end parallel do

    ! Clear variable i_pho, in which we will store the photoionization source term

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a2_box_clear_cc(tree%boxes(id), i_pho)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Add photons to production rate. Currently, this is done sequentially.
    do n = 1, n_used
       id = ph_loc(n)%id
       if (id > a5_no_box) then
          i = ph_loc(n)%ix(1)
          j = ph_loc(n)%ix(2)
          dr = tree%boxes(id)%dr
          tree%boxes(id)%cc(i, j, i_pho) = &
               tree%boxes(id)%cc(i, j, i_pho) + 1/(fac * dr**2)
       end if
    end do

    ! Set ghost cells on highest level with photon source

    !$omp parallel private(lvl, i, id)

    ! Prolong to finer grids
    do lvl = 1, tree%max_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_gc_box_sides(tree%boxes, id, i_pho, &
               a2_sides_interp, a2_bc_neumann)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_prolong1_from(tree%boxes, id, i_pho, f_add=1.0_dp)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine PH_set_src_2d

  subroutine PH_set_src_3d(tree, pi_tbl, rng, num_photons, i_src, i_pho)
    use m_random
    use m_afivo_3d
    use m_lookup_table
    use omp_lib

    type(a3_t), intent(inout)   :: tree   !< Tree
    type(LT_table_t)            :: pi_tbl !< Table to sample abs. lenghts
    type(RNG_t), intent(inout)  :: rng    !< Random number generator
    !> How many discrete photons to use
    integer, intent(in)         :: num_photons
    !> Input variable that contains photon production per cell
    integer, intent(in)         :: i_src
    !> Output variable that contains photoionization source rate
    integer, intent(in)         :: i_pho

    integer                     :: lvl, ix, id, nc
    integer                     :: i, j, k, n, n_create, n_used, i_ph
    integer                     :: proc_id, n_procs
    integer                     :: pho_lvl
    real(dp)                    :: r_create, dr, fac, sum_production, pi_lengthscale
    real(dp), allocatable       :: xyz_src(:, :)
    real(dp), allocatable       :: xyz_dst(:, :)
    type(PRNG_t)                :: prng
    type(a3_loc_t), allocatable :: ph_loc(:)
    real(dp) :: tmp

    nc = tree%n_cell

    ! Get a typical length scale for the absorption of photons
    pi_lengthscale = LT_get_col(pi_tbl, 1, 0.5_dp)

    ! Determine at which level we estimate the photoionization source term. This
    ! depends on the typical lenght scale for absorption.
    pho_lvl = get_lvl_length(tree%dr_base, pi_lengthscale)

    ! Allocate a bit more space because of stochastic production
    allocate(xyz_src(3, nint(1.2_dp * num_photons + 1000)))

    ! Compute the sum of photon production
    call a3_tree_sum_cc(tree, i_src, sum_production, .true.)

    ! Create approximately num_photons
    fac    = num_photons / max(sum_production, epsilon(1.0_dp))
    n_used = 0

    ! Now loop over all leaves and create photons using random numbers

    !$omp parallel private(lvl, ix, id, i, j, dr, i_ph, proc_id, r_create, n_create)
    !$omp single
    n_procs = omp_get_num_threads()
    call prng%init(n_procs, rng)
    !$omp end single

    proc_id = 1+omp_get_thread_num()
    tmp = 0
    do lvl = 1, tree%max_lvl
       dr = a3_lvl_dr(tree, lvl)
       !$omp do
       do ix = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(ix)

          do k = 1, nc
             do j = 1, nc
                do i = 1, nc
                   tmp = tmp + fac * tree%boxes(id)%cc(i, j, k, i_src) * dr**3
                   r_create = fac * tree%boxes(id)%cc(i, j, k, i_src) * dr**3
                   n_create = floor(r_create)

                   if (prng%rngs(proc_id)%uni_01() < r_create - n_create) &
                        n_create = n_create + 1

                   if (n_create > 0) then
                      !$omp critical
                      i_ph = n_used
                      n_used = n_used + n_create
                      !$omp end critical

                      do n = 1, n_create
                         xyz_src(:, i_ph+n) = &
                              a3_r_cc(tree%boxes(id), [i, j, k])
                      end do
                   end if
                end do
             end do
          end do
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    allocate(xyz_dst(3, n_used))
    allocate(ph_loc(n_used))

    ! Get location of absorbption
    call PH_do_absorp(xyz_src, xyz_dst, n_used, pi_tbl, rng)

    !$omp parallel do
    do n = 1, n_used
       ph_loc(n) = a3_get_loc(tree, xyz_dst(1:2, n), pho_lvl)
    end do
    !$omp end parallel do

    ! Clear variable i_pho, in which we will store the photoionization source term

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a3_box_clear_cc(tree%boxes(id), i_pho)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Add photons to production rate. Currently, this is done sequentially.
    do n = 1, n_used
       id = ph_loc(n)%id
       if (id > a5_no_box) then
          i = ph_loc(n)%ix(1)
          j = ph_loc(n)%ix(2)
          k = ph_loc(n)%ix(3)
          dr = tree%boxes(id)%dr
          tree%boxes(id)%cc(i, j, k, i_pho) = &
               tree%boxes(id)%cc(i, j, k, i_pho) + 1/(fac * dr**3)
       end if
    end do

    ! Set ghost cells on highest level with photon source


    !$omp parallel private(lvl, i, id)

    ! Prolong to finer grids
    do lvl = pho_lvl, tree%max_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a3_gc_box_sides(tree%boxes, id, i_pho, &
               a3_sides_interp, a3_bc_neumann)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a3_prolong1_from(tree%boxes, id, i_pho)
       end do
       !$omp end do
    end do
    !$omp end parallel


  end subroutine PH_set_src_3d

end module m_photons
