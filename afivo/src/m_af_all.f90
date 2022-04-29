!> Module which contains all Afivo modules, so that a user does not have to
!> include them separately.
module m_af_all
  use m_af_core
  use m_af_ghostcell
  use m_af_interp
  use m_af_multigrid
  use m_af_output
  use m_af_prolong
  use m_af_restrict
  use m_af_types
  use m_af_utils
  use m_af_particles
  use m_af_flux_schemes
  use m_af_advance
  use m_af_stencil
  use m_af_surface
  use m_coarse_solver

  implicit none
  public

end module m_af_all
