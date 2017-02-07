module m_init_cond_$Dd
  use m_streamer
  use m_a$D_all

  implicit none
  private

  ! Type to store initial conditions in
  type initcnd_t
     real(dp)              :: background_density
     integer               :: n_cond
     real(dp), allocatable :: seed_r0(:, :)
     real(dp), allocatable :: seed_r1(:, :)
     real(dp), allocatable :: seed_density(:)
     integer, allocatable  :: seed_charge_type(:)
     real(dp), allocatable :: seed_width(:)
     integer, allocatable  :: seed_falloff(:)
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected :: ST_init_cond

  public :: field_bc_select
  public :: field_from_potential

contains

  subroutine init_cond_initialize(cfg)
    type(CFG_t), intent(inout)  :: cfg !< Settings

    call CFG_add(cfg, "background_density", 1.0d12, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "seed_density", [5.0d19], &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed_charge_type", [0], &
         "Type of seed: neutral (0), ions (1) or electrons (-1)", .true.)
    call CFG_add(cfg, "seed_rel_r0", [0.5d0, 0.4d0], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", [0.5d0, 0.6d0], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_width", [0.5d-3], &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed_falloff", [1], &
         "Fallof type for seed, see m_geometry.f90", .true.)

  end subroutine init_cond_initialize

end module m_init_cond_$Dd
