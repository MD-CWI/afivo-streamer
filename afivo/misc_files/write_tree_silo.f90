subroutine a$D_write_tree(tree, filename, cc_names, cc_units, n_cycle, time)
  use m_write_silo
  type(a$D_t), intent(in)          :: tree
  character(len=*)                :: filename, cc_names(:), cc_units(:)
  integer, intent(in)             :: n_cycle
  real(dp), intent(in)            :: time
  character(len=*), parameter     :: grid_name = "gg", var_name  = "vv"
  character(len=*), parameter     :: amr_name  = "amr"
  character(len=100), allocatable :: grid_list(:), var_list(:, :)
  integer                         :: lvl, i, id, ig, iv, bs, n_grids, dbix

  bs = tree%box_cells
  n_grids = 0
  do lvl = 1, tree%n_lvls
     n_grids = n_grids + size(tree%lvls(lvl)%ids)
  end do

  allocate(grid_list(n_grids))
  allocate(var_list(tree%n_var_cell, n_grids))

  call SILO_create_file(filename, dbix)
  ig = 0

  do lvl = 1, tree%n_lvls
     do i = 1, size(tree%lvls(lvl)%ids)
        id = tree%lvls(lvl)%ids(i)
        ig = ig + 1
        write(grid_list(ig), "(A,I0)") grid_name, ig
        call SILO_add_grid(dbix, grid_list(ig), 2, &
             [bs+1, bs+1], tree%boxes(id)%r_min, tree%boxes(id)%dr)
        print *, id, tree%boxes(id)%r_min, tree%boxes(id)%dr
        do iv = 1, tree%n_var_cell
           write(var_list(iv, ig), "(A,I0)") trim(cc_names(iv)) // "_", ig
           call SILO_add_var(dbix, var_list(iv, ig), grid_list(ig), &
                pack(tree%boxes(id)%cc(1:bs, 1:bs, iv), .true.), [bs, bs], &
                trim(cc_units(iv)))
        end do
     end do
  end do

  call SILO_set_mmesh_grid(dbix, amr_name, grid_list, n_cycle, time)
  do iv = 1, tree%n_var_cell
     call SILO_set_mmesh_var(dbix, trim(cc_names(iv)), amr_name, &
          var_list(iv, :), n_cycle, time)
  end do
  call SILO_close_file(dbix)
  print *, "Number of grids", ig
end subroutine a$D_write_tree
