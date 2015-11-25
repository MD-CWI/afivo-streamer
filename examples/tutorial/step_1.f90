
program tutorial
  use m_afivo_2d

  implicit none
  integer, parameter  :: dp           = kind(0.0d0)
  integer, parameter  :: box_size     = 8
  integer, parameter  :: n_var_cell   = 1
  integer, parameter  :: n_var_face   = 1

  real(dp), parameter :: domain_len   = 1.0_dp

  type(a2_t)          :: tree

  print *, "This is an Afivo tutorial for solving the heat equation."

  call initialize_tree(tree)

contains

  subroutine initialize_tree(tree)
    type(a2_t), intent(inout) :: tree

    integer, parameter  :: n_boxes_base = 3
    real(dp)            :: dr
    integer             :: ix_list(2, n_boxes_base)
    integer             :: nb_list(4, n_boxes_base)

    dr = domain_len / box_size

    ! Initialize the tree
    call a2_init(tree, & ! Tree to initialize
         box_size, &     ! Number of cells per coordinate in a box
         n_var_cell, &   ! Number of face-centered variables
         n_var_face, &   ! Number of cell-centered variables
         dr)             ! Distance between cells on base level

    ! Set the indices for boxes 1 to 3
    ix_list(:, 1) = [1, 1]
    ix_list(:, 2) = [2, 1]
    ix_list(:, 3) = [2, 2]

    ! Negative values are used for boundary conditions
    nb_list(1:4, 1) = -1
    nb_list(1:4, 2) = -1
    nb_list(1:4, 3) = -1

    ! Connect the two boxes, nb_hx stands for "neighbor in high-x direction".
    ! Connections only have to be specified from one side.
    nb_list(a2_nb_hx, 1) = 2
    nb_list(a2_nb_hy, 2) = 3

    ! Construct the base (level 1) mesh
    call a2_set_base(tree, ix_list, nb_list)

    call a2_print_info(tree)
  end subroutine initialize_tree

end program
