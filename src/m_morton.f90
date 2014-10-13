! This module contains methods to convert indices to morton numbers.
!
! Because fortran does not support unsigned integers, you can only use these
! routines for positive integers.
module m_morton

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: morton_k = selected_int_kind(15)

  ! Public types
  public :: morton_k

  ! Public methods
  public :: morton_from_ix2
  public :: morton_from_ix3
  public :: morton_to_ix2
  public :: morton_to_ix3
  public :: morton_bsearch
  public :: morton_rank
  public :: print_bits
  public :: print_bits_morton

contains

  function morton_from_ix2(ix) result(m_ix)
    integer, intent(in) :: ix(2)
    integer(morton_k)      :: m_ix
    integer(morton_k) :: a, b
    a = bit_space_1(ix(1))
    b = bit_space_1(ix(2))
    m_ix = ior(a, ishft(b, 1))
  end function morton_from_ix2

  function morton_to_ix2(m_ix) result(ix)
    integer(morton_k), intent(in) :: m_ix
    integer                       :: ix(2)
    ix(1) = bit_despace_1(m_ix)
    ix(2) = bit_despace_1(ishft(m_ix, -1))
  end function morton_to_ix2

  function morton_from_ix3(ix) result(m_ix)
    integer, intent(in) :: ix(3)
    integer(morton_k)      :: m_ix
    integer(morton_k) :: a, b, c
    a = bit_space_2(ix(1))
    b = bit_space_2(ix(2))
    c = bit_space_2(ix(3))
    m_ix = ior(ior(a, ishft(b, 1)), ishft(c, 2))
  end function morton_from_ix3

  function morton_to_ix3(m_ix) result(ix)
    integer(morton_k), intent(in) :: m_ix
    integer                       :: ix(3)
    ix(1) = bit_despace_2(m_ix)
    ix(2) = bit_despace_2(ishft(m_ix, -1))
    ix(3) = bit_despace_2(ishft(m_ix, -2))
  end function morton_to_ix3

  function morton_bsearch(list, val) result(ix)
    integer(morton_k), intent(in) :: list(:)
    integer(morton_k), intent(in) :: val
    integer                       :: ix, i_min, i_max, i_middle

    i_min = 1
    i_max = size(list)

    do while (i_min < i_max)
       i_middle = i_min + (i_max - i_min) / 2

       if (val <= list(i_middle)) then
          i_max = i_middle
       else
          i_min = i_middle + 1
       end if
    end do

    ix = i_min
    if (val > list(ix)) ix = -1
  end function morton_bsearch

  ! Add two "zero" bits between each bit of the input. Because the result has 64
  ! bits available, only the first 21 bits from the input are spaced.
  function bit_space_2(a) result(x)
    integer, intent(in) :: a
    integer(morton_k)   :: x

    ! We only look at the first 21 bits
    x = iand(int(a, morton_k),     int(z'1fffff', morton_k))
    x = iand(ior(x, ishft(x, 32)), int(z'1f00000000ffff', morton_k))
    x = iand(ior(x, ishft(x, 16)), int(z'1f0000ff0000ff', morton_k))
    x = iand(ior(x, ishft(x, 8)),  int(z'100f00f00f00f00f', morton_k))
    x = iand(ior(x, ishft(x, 4)),  int(z'10c30c30c30c30c3', morton_k))
    x = iand(ior(x, ishft(x, 2)),  int(z'1249249249249249', morton_k))
  end function bit_space_2

  ! Invert bit_space_2
  function bit_despace_2(a) result(y)
    integer(morton_k), intent(in) :: a
    integer(morton_k)             :: x
    integer                  :: y

    x = iand(a,                     int(z'1249249249249249', morton_k))
    x = iand(ior(x, ishft(x, -2)),  int(z'10c30c30c30c30c3', morton_k))
    x = iand(ior(x, ishft(x, -4)),  int(z'100f00f00f00f00f', morton_k))
    x = iand(ior(x, ishft(x, -8)),  int(z'1f0000ff0000ff', morton_k))
    x = iand(ior(x, ishft(x, -16)), int(z'1f00000000ffff', morton_k))
    x = iand(ior(x, ishft(x, -32)), int(z'1fffff', morton_k))
    y = int(x)
  end function bit_despace_2

  ! Add one "zero" bit between each bit of the input.
  function bit_space_1(a) result(x)
    integer, intent(in) :: a
    integer(morton_k)   :: x

    x = a
    x = iand(ior(x, ishft(x, 32)), int(z'00000000ffffffff', morton_k))
    x = iand(ior(x, ishft(x, 16)), int(z'0000ffff0000ffff', morton_k))
    x = iand(ior(x, ishft(x, 8)),  int(z'00ff00ff00ff00ff', morton_k))
    x = iand(ior(x, ishft(x, 4)),  int(z'0f0f0f0f0f0f0f0f', morton_k))
    x = iand(ior(x, ishft(x, 2)),  int(z'3333333333333333', morton_k))
    x = iand(ior(x, ishft(x, 1)),  int(z'5555555555555555', morton_k))
  end function bit_space_1

  ! Invert bit_space_1
  function bit_despace_1(a) result(y)
    integer(morton_k), intent(in) :: a
    integer(morton_k)             :: x
    integer                  :: y

    x = a
    x = iand(x,                     int(z'5555555555555555', morton_k))
    x = iand(ior(x, ishft(x, -1)),  int(z'3333333333333333', morton_k))
    x = iand(ior(x, ishft(x, -2)),  int(z'0f0f0f0f0f0f0f0f', morton_k))
    x = iand(ior(x, ishft(x, -4)),  int(z'00ff00ff00ff00ff', morton_k))
    x = iand(ior(x, ishft(x, -8)),  int(z'0000ffff0000ffff', morton_k))
    x = iand(ior(x, ishft(x, -16)), int(z'00000000ffffffff', morton_k))
    y = int(x)
  end function bit_despace_1

  ! Jannis: This is taken from ORDERPACK and modified for morton_k
  subroutine morton_rank (XDONT, IRNGT)
    ! __________________________________________________________
    !   MRGRNK = Merge-sort ranking of an array
    !   For performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    Integer(morton_k), Dimension (:), Intent (In)  :: XDONT
    Integer, Dimension (:), Intent (Out) :: IRNGT
    ! __________________________________________________________
    Integer(morton_k) :: XVALA, XVALB
    !
    Integer, Dimension (SIZE(IRNGT)) :: JWRKT
    Integer :: LMTNA, LMTNC, IRNG1, IRNG2
    Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
    NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
    Select Case (NVAL)
    Case (:0)
       Return
    Case (1)
       IRNGT (1) = 1
       Return
    Case Default
       Continue
    End Select
    !
    !  Fill-in the index array, creating ordered couples
    !
    Do IIND = 2, NVAL, 2
       If (XDONT(IIND-1) <= XDONT(IIND)) Then
          IRNGT (IIND-1) = IIND - 1
          IRNGT (IIND) = IIND
       Else
          IRNGT (IIND-1) = IIND
          IRNGT (IIND) = IIND - 1
       End If
    End Do
    If (Modulo(NVAL, 2) /= 0) Then
       IRNGT (NVAL) = NVAL
    End If
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    LMTNA = 2
    LMTNC = 4
    !
    !  First iteration. The length of the ordered subsets goes from 2 to 4
    !
    Do
       If (NVAL <= 2) Exit
       !
       !   Loop on merges of A and B into C
       !
       Do IWRKD = 0, NVAL - 1, 4
          If ((IWRKD+4) > NVAL) Then
             If ((IWRKD+2) >= NVAL) Exit
             !
             !   1 2 3
             !
             If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
             !
             !   1 3 2
             !
             If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                IRNG2 = IRNGT (IWRKD+2)
                IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                IRNGT (IWRKD+3) = IRNG2
                !
                !   3 1 2
                !
             Else
                IRNG1 = IRNGT (IWRKD+1)
                IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                IRNGT (IWRKD+2) = IRNG1
             End If
             Exit
          End If
          !
          !   1 2 3 4
          !
          If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
          !
          !   1 3 x x
          !
          If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
             IRNG2 = IRNGT (IWRKD+2)
             IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
             If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
                !   1 3 2 4
                IRNGT (IWRKD+3) = IRNG2
             Else
                !   1 3 4 2
                IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                IRNGT (IWRKD+4) = IRNG2
             End If
             !
             !   3 x x x
             !
          Else
             IRNG1 = IRNGT (IWRKD+1)
             IRNG2 = IRNGT (IWRKD+2)
             IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
             If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                IRNGT (IWRKD+2) = IRNG1
                If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
                   !   3 1 2 4
                   IRNGT (IWRKD+3) = IRNG2
                Else
                   !   3 1 4 2
                   IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                   IRNGT (IWRKD+4) = IRNG2
                End If
             Else
                !   3 4 1 2
                IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                IRNGT (IWRKD+3) = IRNG1
                IRNGT (IWRKD+4) = IRNG2
             End If
          End If
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 4
       Exit
    End Do
    !
    !  Iteration loop. Each time, the length of the ordered subsets
    !  is doubled.
    !
    Do
       If (LMTNA >= NVAL) Exit
       IWRKF = 0
       LMTNC = 2 * LMTNC
       !
       !   Loop on merges of A and B into C
       !
       Do
          IWRK = IWRKF
          IWRKD = IWRKF + 1
          JINDA = IWRKF + LMTNA
          IWRKF = IWRKF + LMTNC
          If (IWRKF >= NVAL) Then
             If (JINDA >= NVAL) Exit
             IWRKF = NVAL
          End If
          IINDA = 1
          IINDB = JINDA + 1
          !
          !   Shortcut for the case when the max of A is smaller
          !   than the min of B. This line may be activated when the
          !   initial set is already close to sorted.
          !
          !          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
          !
          !  One steps in the C subset, that we build in the final rank array
          !
          !  Make a copy of the rank array for the merge iteration
          !
          JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
          !
          XVALA = XDONT (JWRKT(IINDA))
          XVALB = XDONT (IRNGT(IINDB))
          !
          Do
             IWRK = IWRK + 1
             !
             !  We still have unprocessed values in both A and B
             !
             If (XVALA > XVALB) Then
                IRNGT (IWRK) = IRNGT (IINDB)
                IINDB = IINDB + 1
                If (IINDB > IWRKF) Then
                   !  Only A still with unprocessed values
                   IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                   Exit
                End If
                XVALB = XDONT (IRNGT(IINDB))
             Else
                IRNGT (IWRK) = JWRKT (IINDA)
                IINDA = IINDA + 1
                If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                XVALA = XDONT (JWRKT(IINDA))
             End If
             !
          End Do
       End Do
       !
       !  The Cs become As and Bs
       !
       LMTNA = 2 * LMTNA
    End Do
    !
    Return
    !
  end subroutine morton_rank

  subroutine print_bits(x)
    integer, intent(in) :: x
    integer :: n
    do n = 0, bit_size(x)-1
       if (btest(x, n)) then
          write(*, "(A)", advance="NO") "1"
       else
          write(*, "(A)", advance="NO") "0"
       end if
    end do
    write(*, *) ""
  end subroutine print_bits

  subroutine print_bits_morton(x)
    integer(morton_k), intent(in) :: x
    integer :: n
    do n = 0, bit_size(x)-1
       if (btest(x, n)) then
          write(*, "(A)", advance="NO") "1"
       else
          write(*, "(A)", advance="NO") "0"
       end if
    end do
    write(*, *) ""
  end subroutine print_bits_morton

end module m_morton