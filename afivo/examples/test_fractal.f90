!> \example test_fractal_2d.f90
!> This example shows the basic functionality of m_afivo_2d.
program test_fractal
  use m_afivo_2d

  implicit none

  integer, parameter   :: dp           = kind(0.0d0)
  type(a2_t)           :: tree
  integer              :: i
  integer, parameter   :: n_boxes_base = 1
  integer              :: ix_list(2, n_boxes_base)
  integer              :: nb_list(4, n_boxes_base)
  integer, parameter   :: box_size     = 2
  integer, parameter   :: i_phi        = 1, i_mrtn = 2
  integer, parameter   :: n_var_cell   = 1
  integer, parameter   :: n_var_face   = 0
  type(ref_info_t)     :: ref_info
  character(len=40)    :: var_names(1) = ["phi"]
  real(dp)             :: dr
  character(len=40)    :: fname

  dr = 0.25_dp / box_size

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! Number of cells per coordinate in a box
       n_var_cell, &   ! Number of face-centered variables
       n_var_face, &   ! Number of cell-centered variables
       dr, &           ! Distance between cells on base level
       r_min = [-0.235_dp, 0.93_dp])

  print *, tree%r_base

  ! Set up geometry
  ix_list(:, 1) = [1,1] ! One box at 1,1

  ! Set neighbors for box one
  nb_list(:, 1) = -1

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)

  ! Set variables on base
  call a2_loop_box(tree, set_phi_box)
  call a2_gc_sides(tree, i_phi, a2_sides_interp, a2_bc_neumann)

  do i = 1, 12
     print *, "i = ", i, "max_id", tree%max_id

     write(fname, "(A,I0)") "test_fractal_2d_", i
     call a2_write_vtk(tree, trim(fname), var_names, i, i * 1.0_dp)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     call prolong_to_new_children(tree, ref_info)
     call a2_tidy_up(tree, max_frac_used=0.75_dp, goal_frac_used=0.5_dp, &
          n_clean_min=10000, reorder=.true.)
  end do

  call a2_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer :: nc
    real(dp)                 :: pmax, pmin, diff

    nc = boxes(id)%n_cell
    pmax = maxval(boxes(id)%cc(1:nc, 1:nc, i_phi))
    pmin = minval(boxes(id)%cc(1:nc, 1:nc, i_phi))

    diff = max( &
         maxval(abs(boxes(id)%cc(1:nc+1, 1:nc, i_phi) - &
         boxes(id)%cc(0:nc, 1:nc, i_phi))), &
         maxval(abs(boxes(id)%cc(1:nc, 1:nc+1, i_phi) - &
         boxes(id)%cc(1:nc, 0:nc, i_phi))))

    if (boxes(id)%lvl < 3) then
       ref_flags(id) = a5_do_ref
    else if (diff * boxes(id)%dr > 1.0e-2_dp .and. boxes(id)%lvl < 14) then
       ref_flags(id) = a5_do_ref
    end if
  end subroutine set_ref_flags

  subroutine set_phi_box(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)
    complex(dp) :: z0

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          z0 = cmplx(xy(1), xy(2), dp)
          box%cc(i, j, i_phi) = get_num_its(z0)
       end do
    end do
  end subroutine set_phi_box

  real(dp) function get_num_its(z0)
    complex(dp), intent(in) :: z0
    complex(dp)             :: zz
    real(dp), parameter :: radius = 2.0_dp
    integer, parameter  :: max_its = 10*1000

    zz          = 0
    get_num_its = 0
    do
       if (abs(zz) > radius .or. get_num_its >= max_its) exit
       zz          = zz**2 + z0
       get_num_its = get_num_its + 1
    end do
  end function get_num_its

  subroutine prolong_to_new_children(tree, ref_info)
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       !$omp parallel do private(id)
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call set_phi_box(tree%boxes(id))
       end do
       !$omp end parallel do

       !$omp parallel do private(id)
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box_sides(tree%boxes, id, i_phi, &
               a2_sides_interp, a2_bc_neumann)
       end do
       !$omp end parallel do
    end do
  end subroutine prolong_to_new_children

  subroutine set_morton_variable(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    boxes(id)%cc(:,:,i_mrtn) = id
  end subroutine set_morton_variable

end program test_fractal
