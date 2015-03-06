module m_photons

  implicit none
  private

  integer, parameter                     :: dp = kind(0.0d0)

  ! Public methods
  public :: PH_get_tbl_air
  public :: PH_do_absorp
  public :: PH_set_src

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
       incr = dr * sixth * (absfunc_air(r_prev, p_O2) + &
            4 * absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
            absfunc_air(r_prev+dr, p_O2))

       do while (incr > 2.0_dp / tbl_size)
          dr = dr * 0.5_dp
          incr = dr * sixth * (absfunc_air(r_prev, p_O2) + &
               4 * absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
               absfunc_air(r_prev+dr, p_O2))
       end do

       do while (incr < 1.0_dp / (tbl_size-1) .and. dr < huge_dist)
          dr = dr * 2
          incr = dr * sixth * (absfunc_air(r_prev, p_O2) + &
               4 * absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
               absfunc_air(r_prev+dr, p_O2))
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

  real(dp) function absfunc_air(dist, p_O2)
    use m_units_constants
    real(dp), intent(in) :: dist, p_O2
    real(dp), parameter  :: chi_min  = 3.5_dp / UC_torr_to_bar
    real(dp), parameter  :: chi_max  = 200 / UC_torr_to_bar

    absfunc_air = (exp(-chi_min * p_O2 * dist) - &
         exp(-chi_max * p_O2 * dist)) / (dist * log(chi_max/chi_min))
  end function absfunc_air

  subroutine PH_do_absorp(xyz_in, xyz_out, n_photons, tbl, rng)
    use m_lookup_table
    use m_random
    integer, intent(in)          :: n_photons
    real(dp), intent(in)         :: xyz_in(3, n_photons)
    real(dp), intent(out)        :: xyz_out(3, n_photons)
    type(LT_table_t), intent(in) :: tbl
    type(RNG_t), intent(inout)   :: rng
    integer                      :: n
    real(dp)                     :: rr, dist

    !!$omp parallel do private(rr, dist) firstprivate(rng)
    do n = 1, n_photons
       rr = rng%uni_01()
       dist = LT_get_col(tbl, 1, rr)
       xyz_out(:, n) =  xyz_in(:, n) + rng%sphere(dist)
    end do
    !!$omp end parallel do
  end subroutine PH_do_absorp

  subroutine PH_set_src(tree, pi_tbl, num_photons, i_pho)
    use m_random
    use m_afivo_2d
    use m_lookup_table

    type(a2_t), intent(inout) :: tree   !< Tree
    type(LT_table_t)          :: pi_tbl !< Table to sample abs. lenghts
    !> How many discrete photons to use
    integer, intent(in)       :: num_photons
    !> Index of variable that contains photon production per cell
    integer, intent(in)       :: i_pho

    integer :: lvl, ix, id, nc
    integer :: i, j, n, n_create, n_used, i_ph
    integer :: pho_lvl
    real(dp) :: r_create, dr, fac, sum_production, pi_lengthscale
    real(dp), allocatable :: xyz_src(:, :)
    real(dp), allocatable :: xyz_dst(:, :)
    type(RNG_t) :: rng
    type(a2_loc_t), allocatable :: ph_loc(:)

    nc = tree%n_cell

    ! Get a typical length scale for the absorption of photons
    pi_lengthscale = LT_get_col(pi_tbl, 1, 0.5_dp)

    ! Determine at which level we estimate the photoionization source term. This
    ! depends on the typical lenght scale for absorption.
    do lvl = 1, tree%lvls_max
       if (a2_lvl_dr(tree, lvl) < pi_lengthscale) exit
    end do
    pho_lvl = lvl

    ! Allocate a bit more space because of stochastic production
    allocate(xyz_src(3, nint(1.2_dp * num_photons + 1000)))

    ! Compute the sum of photon production
    call a2_tree_sum_cc(tree, i_pho, sum_production)

    ! Create approximately num_photons
    fac = num_photons/sum_production
    n_used = 0
    ! print *, "num_photons", num_photons
    ! print *, "sum_production", sum_production

    ! Now loop over all leaves and create photons using random numbers

    !$omp parallel private(lvl, ix, id, i, j, i_ph, r_create, n_create) &
    !$omp & firstprivate(rng)
    do lvl = 1, tree%max_lvl
       !$omp do
       do ix = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(ix)

          do j = 1, nc
             do i = 1, nc
                r_create = fac * tree%boxes(id)%cc(i, j, i_pho)
                n_create = floor(r_create)

                if (rng%uni_01() < r_create - n_create) &
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

    !$omp parallel do
    do n = 1, n_used
       ph_loc(n) = a2_get_loc(tree, xyz_dst(1:2, n), pho_lvl)
    end do
    !$omp end parallel do

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
    ! !$omp do
    ! do i = 1, size(tree%lvls(pho_lvl)%ids)
    !    id = tree%lvls(pho_lvl)%ids(i)
    !    call a2_gc_box_sides(tree%boxes, id, i_pho, &
    !         a2_sides_extrap, sides_neumann)
    ! end do
    ! !$omp end do

    ! Prolong to finer grids
    do lvl = pho_lvl, tree%max_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_prolong0_from(tree%boxes, id, i_pho, .false.)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine PH_set_src

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine sides_neumann(boxes, id, nb, iv)
    use m_afivo_2d
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       boxes(id)%cc(1:nc, 0, iv) = boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       boxes(id)%cc(1:nc, nc+1, iv) = boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_neumann

end module m_photons
