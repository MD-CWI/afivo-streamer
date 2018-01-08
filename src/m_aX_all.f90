!> Module which contains all Afivo modules, so that a user does not have to
!> include them separately.
module m_a$D_all
  use m_a$D_core
  use m_a$D_ghostcell
  use m_a$D_interp
  use m_a$D_multigrid
  use m_a$D_output
  use m_a$D_prolong
  use m_a$D_restrict
  use m_a$D_types
  use m_a$D_utils
  use m_a$D_particles

  implicit none
  public

end module m_a$D_all
