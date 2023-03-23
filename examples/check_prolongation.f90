#include "../src/cpp_macros.h"
!> \example check_prolongation
!>
!> This example checks whether prolongation behaves as expected
program check_prolongation
  use m_af_all
  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_runs   = 50
  integer            :: i_phi
  real(dp)           :: gradient(NDIM)
  logical, parameter :: write_silo = .false.

  write(*,'(A,I0,A)') 'program check_prolongation_', NDIM, "d"
  print *, "Number of threads", af_get_max_threads()

  ! Check whether prolongation is exact for a linear gradient
  call test_gradient(n_runs, af_prolong_linear, "linear")
  call test_gradient(n_runs, af_prolong_sparse, "sparse")
  call test_gradient(n_runs, af_prolong_limit, "limit")

  ! Check whether values lie in a valid range after prolongation
  call test_valid_range(n_runs, af_prolong_linear, "linear")
  call test_valid_range(n_runs, af_prolong_sparse, "sparse")
  call test_valid_range(n_runs, af_prolong_zeroth, "zeroth")
  call test_valid_range(n_runs, af_prolong_limit, "limit")

  print *, "Tests passed"

contains

  subroutine test_valid_range(max_iter, prolongation_method, name)
    integer, intent(in)          :: max_iter
    procedure(af_subr_prolong)   :: prolongation_method
    character(len=*), intent(in) :: name
    type(af_t)                   :: tree
    integer                      :: iter
    type(ref_info_t)             :: ref_info
    character(len=100)           :: fname

    print *, "Testing range of values (", name, ")"

    call af_add_cc_variable(tree, "phi", ix=i_phi)
    call af_set_cc_methods(tree, i_phi, af_bc_dirichlet_zero, &
         prolong=prolongation_method)

    call af_init(tree, box_size, [DTIMES(1.0_dp)], &
         [DTIMES(box_size)], periodic=[DTIMES(.true.)], &
         mem_limit_gb=0.5_dp)

    do iter = 1, max_iter
       call af_loop_box(tree, set_random_state)
       call af_restrict_tree(tree, [i_phi])
       call af_gc_tree(tree, [i_phi])

       call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=0)

       call af_loop_tree(tree, check_range_box, leaves_only=.true.)

       if (write_silo) then
          write(fname, "(A,I0)") "output/check_prolongation_" // name &
               // "_" // DIMNAME // "_", iter
          call af_write_silo(tree, trim(fname), n_cycle=iter)
       end if
    end do

  end subroutine test_valid_range

  subroutine test_gradient(max_iter, prolongation_method, name)
    integer, intent(in)          :: max_iter
    procedure(af_subr_prolong)   :: prolongation_method
    character(len=*), intent(in) :: name
    type(af_t)                   :: tree
    integer                      :: iter
    type(ref_info_t)             :: ref_info
    character(len=100)           :: fname

    print *, "Testing linear gradient (", name, ")"

    call af_add_cc_variable(tree, "phi", ix=i_phi)
    call af_set_cc_methods(tree, i_phi, bc_gradient, &
         prolong=prolongation_method)

    call af_init(tree, box_size, [DTIMES(1.0_dp)], &
         [DTIMES(box_size)], periodic=[DTIMES(.false.)], &
         mem_limit_gb=0.5_dp)

    do iter = 1, max_iter
       call random_number(gradient)
       call af_loop_box(tree, set_solution_gradient)
       call af_restrict_tree(tree, [i_phi])
       call af_gc_tree(tree, [i_phi])

       call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=0)

       call af_loop_tree(tree, check_gradient_box, leaves_only=.true.)

       if (write_silo) then
          write(fname, "(A,I0)") "output/check_prolongation_gradient_" &
               // name // "_" // DIMNAME // "_", iter
          call af_write_silo(tree, trim(fname), n_cycle=iter)
       end if
    end do
  end subroutine test_gradient

  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: rr

    call random_number(rr)

    if (rr < 0.5_dp**NDIM .and. box%lvl < 5) then
       cell_flags = af_do_ref
    else
       cell_flags = af_rm_ref
    end if
  end subroutine ref_routine

  subroutine set_random_state(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    nc = box%n_cell
    call random_number(box%cc(DTIMES(1:nc), i_phi))
  end subroutine set_random_state

  subroutine set_solution_gradient(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc, IJK
    real(dp)                   :: r(NDIM)

    nc = box%n_cell

    do KJI_DO(1,nc)
       r = af_r_cc(box, [IJK])
       box%cc(IJK, i_phi) = sum(gradient * r)
    end do; CLOSE_DO
  end subroutine set_solution_gradient

  subroutine check_range_box(tree, id)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer                   :: nc
    real(dp)                  :: min_phi, max_phi
    logical                   :: mask(DTIMES(0:tree%n_cell+1))
    integer                   :: minloc_phi(NDIM), maxloc_phi(NDIM)

    nc = tree%n_cell

    ! Set the mask for the interior and the sides (excluding corners)
    mask                    = .false.
    mask(DTIMES(1:nc))      = .true.
#if NDIM == 1
    mask(:0)                = .true.
    mask(nc+1:)             = .true.
#elif NDIM == 2
    mask(:0, 1:nc)          = .true.
    mask(nc+1:, 1:nc)       = .true.
    mask(1:nc, :0)          = .true.
    mask(1:nc, nc+1:)       = .true.
#elif NDIM == 3
    mask(:0, 1:nc, 1:nc)    = .true.
    mask(nc+1:, 1:nc, 1:nc) = .true.
    mask(1:nc, :0, 1:nc)    = .true.
    mask(1:nc, nc+1:, 1:nc) = .true.
    mask(1:nc, 1:nc, :0)    = .true.
    mask(1:nc, 1:nc, nc+1:) = .true.
#endif

    associate (cc => tree%boxes(id)%cc)
      min_phi = minval(cc(DTIMES(:), 1), mask=mask)
      minloc_phi = minloc(cc(DTIMES(:), 1), mask=mask)
      max_phi = maxval(cc(DTIMES(:), 1), mask=mask)
      maxloc_phi = maxloc(cc(DTIMES(:), 1), mask=mask)

      if (min_phi < 0 .or. max_phi > 1) then
         print *, "minimum of phi: ", min_phi, minloc_phi - 1
         print *, "maximum of phi: ", max_phi, maxloc_phi - 1
         error stop "min/max not in range [0, 1]"
      end if
    end associate
  end subroutine check_range_box

  subroutine check_gradient_box(tree, id)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer                   :: IJK, nc
    logical                   :: mask(DTIMES(0:tree%n_cell+1))
    real(dp)                  :: diff(DTIMES(0:tree%n_cell+1))
    real(dp)                  :: max_deviation, r(NDIM)
    real(dp), parameter       :: threshold = 1e-13_dp

    nc = tree%n_cell

    ! Set the mask for the interior and the sides (excluding corners)
    mask                    = .false.
    mask(DTIMES(1:nc))      = .true.
#if NDIM == 1
    mask(:0)                = .true.
    mask(nc+1:)             = .true.
#elif NDIM == 2
    mask(:0, 1:nc)          = .true.
    mask(nc+1:, 1:nc)       = .true.
    mask(1:nc, :0)          = .true.
    mask(1:nc, nc+1:)       = .true.
#elif NDIM == 3
    mask(:0, 1:nc, 1:nc)    = .true.
    mask(nc+1:, 1:nc, 1:nc) = .true.
    mask(1:nc, :0, 1:nc)    = .true.
    mask(1:nc, nc+1:, 1:nc) = .true.
    mask(1:nc, 1:nc, :0)    = .true.
    mask(1:nc, 1:nc, nc+1:) = .true.
#endif

    associate (cc => tree%boxes(id)%cc)
      ! Set solution
      diff = 0.0_dp
      do KJI_DO(0, nc+1)
         if (mask(IJK)) then
            r = af_r_cc(tree%boxes(id), [IJK])
            diff(IJK) = abs(sum(gradient * r) - cc(IJK, 1))
         end if
      end do; CLOSE_DO

      max_deviation = maxval(diff)

      if (max_deviation > threshold) then
         print *, "max deviation from linear solution: ", max_deviation
         print *, "location: ", maxloc(diff) - 1
         error stop "max deviation larger than threshold"
      end if
    end associate
  end subroutine check_gradient_box

  !> Set boundary conditions for a linear gradient solution
  subroutine bc_gradient(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: nc

    nc = box%n_cell
    bc_type = af_bc_dirichlet
    bc_val = matmul(gradient, coords)
  end subroutine bc_gradient

end program
