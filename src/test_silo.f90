program test_silo
  use m_write_silo

  implicit none

  integer, parameter             :: dp        = kind(0.0d0)
  integer, parameter             :: NN        = 16
  character(len=*), parameter    :: filename  = "test.silo"
  character(len=*), parameter    :: grid_name = "gg"
  character(len=*), parameter    :: var_name  = "vv"
  character(len=*), parameter    :: amr_name  = "amr"
  integer                        :: ix, n_grids
  character(len=40), allocatable :: mesh_name_list(:), var_name_list(:)

  type grid_t
     real(dp) :: r_min(3)
     real(dp) :: dr(3)
     real(dp) :: var(NN, NN, NN)
  end type grid_t

  type(grid_t) :: grid_list(5)
  type(grid_t) :: grid

  ! Create 5 blocks; 1 coarse with 4 fine next to it
  grid_list(1)%r_min = 0
  grid_list(1)%dr = 1.0_dp/(NN-1)
  grid_list(2)%dr = 0.5_dp/(NN-1)
  grid_list(2)%r_min = (/1.0_dp, 0.0_dp, 0.0_dp/)
  grid_list(3)%dr = 0.5_dp/(NN-1)
  grid_list(3)%r_min = (/1.0_dp, 0.5_dp, 0.0_dp/)
  grid_list(4)%dr = 0.5_dp/(NN-1)
  grid_list(4)%r_min = (/1.0_dp, 0.0_dp, 0.5_dp/)
  grid_list(5)%dr = 0.5_dp/(NN-1)
  grid_list(5)%r_min = (/1.0_dp, 0.5_dp, 0.5_dp/)

  n_grids = size(grid_list)
  do ix = 1, n_grids
     call fill_grid(grid_list(ix))
  end do

  call SILO_create_file(filename)
  allocate(mesh_name_list(n_grids))
  allocate(var_name_list(n_grids))

  do ix = 1, n_grids
     grid = grid_list(ix)
     write(mesh_name_list(ix), "(A,I0)") grid_name, ix
     call SILO_add_grid(filename, mesh_name_list(ix), 3, &
          (/NN, NN, NN/), grid%r_min, grid%dr)
     write(var_name_list(ix), "(A,I0)") trim(var_name) // "_", ix
     call SILO_add_var(filename, var_name_list(ix), mesh_name_list(ix), &
          pack(grid%var, .true.), (/NN, NN, NN/), "UNIT")
  end do

  call SILO_set_mmesh_grid(filename, amr_name, mesh_name_list, 1, 0.0d0)
  call SILO_set_mmesh_var(filename, trim(var_name), amr_name, &
       var_name_list(:), 1, 0.0d0)

contains

  subroutine fill_grid(grid)
    type(grid_t), intent(inout) :: grid
    integer :: i, j, k
    real(dp) :: x, y, z
    do k = 1, NN
       z = grid%r_min(3) + (k-1) * grid%dr(3)
       do j = 1, NN
          y = grid%r_min(2) + (j-1) * grid%dr(2)
          do i = 1, NN
             x = grid%r_min(1) + (i-1) * grid%dr(1)

             grid%var(i, j, k) = sin(x**2 + y**2 + z**2)
          end do
       end do
    end do
  end subroutine fill_grid

end program test_silo