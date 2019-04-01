program test_m_morton
  use m_morton

  implicit none

  integer(morton_k) :: m_ix
  integer           :: ix(NDIM)

#if NDIM == 2
  ix = (/15, 2047/)
#elif NDIM == 3
  ix = (/15, 2047, 2047/)
#endif
  print *, "index       ", ix
#if NDIM == 2
  m_ix = morton_from_ix2(ix)
#elif NDIM == 3
  m_ix = morton_from_ix3(ix)
#endif
  print *, "morton index", m_ix
#if NDIM == 2
  ix = morton_to_ix2(m_ix)
#elif NDIM == 3
  ix = morton_to_ix3(m_ix)
#endif
  print *, "index again ", ix
  call print_bits(ix(1))
  call print_bits(ix(2))
  call print_bits_morton(m_ix)
end program test_m_morton
