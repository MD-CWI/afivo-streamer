module m_io

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: io_write_2d

contains

  subroutine io_write_2d(tree, filename, cc_names, cc_units, n_cycle, time)
    use m_write_silo
    type(a5_2d_t), intent(in) :: tree
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: cc_names(:)
    character(len=*), intent(in) :: cc_units(:)
    integer, intent(in) :: n_cycle
    real(dp), intent(in) :: time

    character(len=*), parameter :: mesh_name = "mesh"
    character(len=80), allocatable :: mesh_names(:), var_names(:,:)
    integer :: lvl, i, ig, iv
    integer :: n_grids, n_vars

    n_vars = size(cc_names)

    n_grids = 0
    do lvl = 1, tree%n_levels
       n_grids = n_grids + size(tree%levels(lvl)%ids)
    end do

    allocate(mesh_names(n_grids))
    allocate(var_names(n_vars, n_grids))

    call SILO_create_file(filename)
    ig = 0

    do lvl = 1, tree%n_levels
       do i = 1, size(tree%levels(lvl)%ids)
          ig = ig + 1
          write(mesh_names(ig), "(A,I0)") mesh_name, ig
          call SILO_add_grid(filename, mesh_names(ig), 2, &
               (/tree%n_points, tree%n_points/), &
               a2_get_rmin(b_ix), tree%dr(lvl))
          do iv = 1, n_vars
             write(var_names(iv, ig), "(A,I0)") trim(cc_names(iv)) // "_", ig
             call SILO_add_var(filename, var_names(iv, ig), mesh_names(ig), &
                  pack(grid%var, .true.), (/tree%n_points, tree%n_points/), &
                  trim(cc_units(iv)))
          end do
       end do
    end do

    call SILO_set_mmesh_grid(filename, mesh_name, mesh_names, n_cycle, time)
    do iv = 1, n_vars
       call SILO_set_mmesh_var(filename, trim(cc_names(iv)), mesh_name, &
         var_names, n_cycle, time)
    end do
  end subroutine io_write_2d

end module m_io