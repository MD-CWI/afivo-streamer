#include "../afivo/src/cpp_macros_$Dd.h"
!> Top-module for photoionization, which can make use of different methods
module m_photoi_$Dd
  use m_photoi_mc
  use m_photoi_helmh_$Dd
  use m_a$D_all
  use m_streamer

  implicit none
  private

  ! Whether photoionization is enabled
  logical, protected, public :: photoi_enabled = .false.

  ! Which photoionization method to use (helmholtz, montecarlo)
  character(len=ST_slen) :: photoi_method = 'helmholtz'

  ! Photoionization efficiency factor, typically around 0.05-0.1, not for Helmholtz-Luque should be 1.0
  real(dp) :: photoi_eta = 0.1_dp

  ! Update photoionization every N time step
  integer, protected, public :: photoi_per_steps = 10
  public :: photoi_initialize
  public :: photoi_set_src
  public :: define_surf_photo
  public :: photo_absorb
  
  ! Imported from Helmholtz module
  public :: photoi_helmh_bc

contains

  !> Initialize photoionization parameters
  subroutine photoi_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg          !< The configuration for the simulation

    call CFG_add_get(cfg, "photoi%enabled", photoi_enabled, &
         "Whether photoionization is enabled")
    call CFG_add_get(cfg, "photoi%per_steps", photoi_per_steps, &
         "Update photoionization every N time step")
    call CFG_add_get(cfg, "photoi%method", photoi_method, &
         "Which photoionization method to use (helmholtz, montecarlo)")
    call CFG_add_get(cfg, "photoi%eta", photoi_eta, &
         "Photoionization efficiency factor, typically around 0.05-0.1")
         
    if (photoi_eta <= 0.0_dp) error stop "photoi%eta <= 0.0"
    if (photoi_eta > 1.0_dp) error stop "photoi%eta > 1.0"


    if (photoi_enabled) then
       i_photo        = ST_add_cc_variable("photo", .true.)
       i_photo_low    = ST_add_cc_variable("photo_low", .true.)
       surf_photo     = ST_add_fc_variable("surf_photo", .true.)
    end if

    select case (photoi_method)
       case ("helmholtz")
          call photoi_helmh_initialize(cfg, .true.)
          call phmc_initialize(cfg, .false.)
       case ("montecarlo")
          call photoi_helmh_initialize(cfg, .false.)
          call phmc_initialize(cfg, .true.)
       case default
          print *, "Unknown photoi_method: ", trim(photoi_method)
          error stop
    end select
  end subroutine photoi_initialize
  

  !> Sets the photoionization
  subroutine photoi_set_src(tree, dt)
    use m_units_constants

    type(a$D_t), intent(inout)     :: tree
    real(dp), intent(in), optional :: dt
    real(dp), parameter            :: p_quench = 40.0e-3_dp
    real(dp)                       :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a$D_loop_box_arg(tree, set_photoionization_rate, &
         [photoi_eta * quench_fac], .true.)

    select case (photoi_method)
    case ("helmholtz")
       ! Use Helmholtz approximation
       call photoi_helmh_compute(tree)
    case ("montecarlo")
       if (phmc_physical_photons) then
#if $D == 2
          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, &
               i_photo, ST_cylindrical, .false., dt)
          call phmc_set_src_$Dd(tree, ST_rng, i_pos_ion_old, &
               i_photo_low, ST_cylindrical, .true., dt)              
!#elif $D == 3
!          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, i_photo, .false., dt)
 !         call phmc_set_src_$Dd(tree, ST_rng, i_pos_ion_old, i_photo_low, .true., dt)
#endif
       else
#if $D == 2
          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, &
               i_photo, ST_cylindrical, .false.)
          call phmc_set_src_$Dd(tree, ST_rng, i_pos_ion_old, &
               i_photo_low, ST_cylindrical, .true.)     
!#elif $D == 3
!          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, i_photo, .false.)
!          call phmc_set_src_$Dd(tree, ST_rng, i_pos_ion_old, i_photo_low, .true.)
#endif
       end if

    end select

  end subroutine photoi_set_src

  !> Sets the photoionization_rate
  subroutine set_photoionization_rate(box, coeff)
    use m_units_constants
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)         :: coeff(:)
    integer                      :: IJK, nc
    real(dp)                     :: fld, alpha, mobility, k_N2, tmp, freq_0, tau
    type(LT_loc_t)               :: loc

    nc = box%n_cell

    do KJI_DO(1,nc)
       fld      = box%cc(IJK, i_electric_fld)
       loc      = LT_get_loc(ST_td_tbl, fld)
       alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
       mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)


       tmp = fld * mobility * alpha * box%cc(IJK, i_electron) * coeff(1)
       if (tmp < 0) tmp = 0
       box%cc(IJK, i_electron_old) = tmp     
       k_N2          = 0.0_dp!LT_get_col_at_loc(ST_ex_tbl, i_ex_N2, loc)
       freq_0        = 1.0_dp/UC_N2_tau0
       tau           = 1.0_dp/(freq_0 + UC_N2_N2_qfr + UC_N2_O2_qfr)
       box%cc(IJK, i_pos_ion_old) = k_N2 * box%cc(IJK, i_electron) * 2.5e25_dp * tau * freq_0
    end do; CLOSE_DO
  end subroutine set_photoionization_rate
  
  function photo_DI_absorption_box(box) result(phot)
    type(box$D_t), intent(in)       :: box
    real(dp)                        :: phot(2)
    integer                         :: lvl, li, id, nc, dir, IJK, ix($D)
    
    phot(:) = 0.0_dp
    nc  = box%n_cell
    if (maxval(box%cc(DTIMES(:), i_eps)) > 1.0_dp) then
      do KJI_DO(1, nc)
        if (box%cc(IJK, i_eps) == ST_epsilon_die) then
          phot(1) = phot(1) + (box%cc(IJK, i_photo) + &
                    box%cc(IJK, i_photo_low) ) * box%dr**$D
          do dir = 1, a$D_num_neighbors
            ix = a$D_neighb_dix(:, dir)
#if $D == 2
            if (box%cc(i+ix(1), j+ix(2), i_eps) < ST_epsilon_die) then
              phot(2) = phot(2) + box%cc(IJK, i_photo)
              exit
            end if
#endif
          end do
        end if
      end do; CLOSE_DO
      
    end if

  end function photo_DI_absorption_box

  
  
  subroutine photo_absorb(tree, out_1, out_2)
    type(a$D_t), intent(inout) :: tree    
    real(dp), intent(out)      :: out_1, out_2
    real(dp)                   :: my_val(2), f_dummy(2)
    integer                    :: i, id, lvl

    call a$D_tree_clear_fc(tree, surf_photo)
    if (.not. tree%ready) stop "Tree not ready"
    my_val(:) = 0.0_dp
    out_1 = 0.0_dp
    out_2 = 0.0_dp

    !$omp parallel private(lvl, i, id, f_dummy) firstprivate(my_val)
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          f_dummy = photo_DI_absorption_box(tree%boxes(id))
          my_val(:) = my_val(:) + f_dummy(:)
       end do
       !$omp end do
    end do

    !$omp critical
    out_1 =  out_1 + my_val(1)
    out_2 =  out_2 + my_val(2)
    !$omp end critical
    !$omp end parallel
  end subroutine photo_absorb
  
  subroutine define_surf_photo(box, photo)
    type(box$D_t), intent(inout)   :: box
    real(dp), intent(in)           :: photo(:)
    integer                        :: IJK, nc, s_count, ix($D), dir, i_count
    real(dp)                       :: dist_amount, i_fac
      
    nc  = box%n_cell 
    if (maxval(box%cc(DTIMES(:), i_eps)) > 1.0_dp) then    
      do KJI_DO(1, nc)
        if (box%cc(IJK, i_eps) == ST_epsilon_die) then
          s_count = 0
          i_count = 0
          do dir = 1, a$D_num_neighbors
            ix = a$D_neighb_dix(:, dir)
#if $D == 2
            if (box%cc(i+ix(1), j+ix(2), i_eps) < ST_epsilon_die) then
              s_count = s_count + 1
            end if
#endif
          end do
!            dist_amount = 2.0_dp * box%cc(IJK, i_photo) * ST_phe_yield_high 

 !           do dir = 1, a$D_num_neighbors ! To smoothen the surface distribution we insert at surf_photo, using neighbour cells
  !            ix = a$D_neighb_dix(:, dir)
!#if $D == 2
!              if (box%cc(i+ix(1), j+ix(2), i_eps) == ST_epsilon_die) then
!                dist_amount = dist_amount + box%cc(i+ix(1), j+ix(2), i_photo) * ST_phe_yield_high
!                i_count = i_count + 1
 !             end if
!#endif
 !           end do
  !        i_fac = 1.0_dp/(2 + i_count)
   !       dist_amount = dist_amount * i_fac
          
          do dir = 1, a$D_num_neighbors
            ix = a$D_neighb_dix(:, dir)
#if $D == 2

            if (box%cc(i+ix(1), j+ix(2), i_eps) < ST_epsilon_die .and. s_count > 0) then

              box%fc(max(i, i+ix(1)), max(j, j+ix(2)), int((dir-1)/2)+1, surf_photo) = box%cc(IJK, i_photo) * photo(1) / &
                                                                        (epsilon(1.0_dp) + photo(2) * s_count * box%dr**2)
            end if
#endif
          end do 
        end if
      end do; CLOSE_DO
    end if
    
  end subroutine define_surf_photo
  
  

end module m_photoi_$Dd