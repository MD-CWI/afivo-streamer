program test_m_morton
  use m_morton

  implicit none

  integer(morton_k) :: m_ix
  integer           :: ix(2)

  ix = (/15, 2047/)
  print *, "index       ", ix
  m_ix = morton_from_ix2(ix)
  print *, "morton index", m_ix
  ix = morton_to_ix2(m_ix)
  print *, "index again ", ix
  call print_bits(ix(1))
  call print_bits(ix(2))
  call print_bits_morton(m_ix)
end program test_m_morton