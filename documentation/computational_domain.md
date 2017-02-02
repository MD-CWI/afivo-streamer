# Defining the computational domain

In Afivo, the computational domain is defined by the coarse grid. 

## Step 1: calling a2_init()

First a user has to call `m_a2_core::a2_init()`, for example as shown below.

```{f90}
    integer    :: n_cell   = 4      ! Boxes contain n_cell^dim cells
    real(dp)   :: dr       = 1.0_dp ! Coarse grid spacing
    integer    :: n_var_cc = 1      ! Number of cell-centered variables
    integer    :: n_var_fc = 0      ! Number of face-centered variables
    type(a2_t) :: tree             ! Will contain the quad/octree grid

    call a2_init(tree, n_cell, n_var_cc, n_var_fc, dr)
```

The coarse boxes on this quadtree contain \f$4 \times 4\f$ cells, and these
cells have spacing of `dr = 1.0`.

## Step 2: Setting the box spatial indices

The user can now choose how many boxes to add to the computational domain. For
each box, a *spatial index* has to be given, which will determine where the box
is located. The lowest index a box can have is [1, 1]. A box with index [2, 1]
lies to the right (+x direction) of the box at [1, 1], etc. One could for
example do the following:

```{f90}
    integer, parameter  :: n_boxes = 2
    integer             :: ix_list(2, n_boxes)

    ! Two boxes along x-direction
    ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
    ix_list(:, 2) = [2, 1]        ! Box 2 at [2, 1]
```
The coordinates of these boxes would then be

box | x,y min | x,y max
---|---|---
1 | 0, 0 | 4, 4
2 | 4, 0 | 8, 4

To specify two boxes that touch in the y-direction, one would similarly do:

```{f90}
    ! Two boxes along y-direction
    ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
    ix_list(:, 2) = [1, 2]        ! Box 2 at [1, 2]
```

## Step 3: Optionally setting the neighbor information

Now the user can call `m_a2_core::a2_set_base` to create the computational
domain:

    call a2_set_base(tree, n_boxes, ix_list)

Afivo will automatically resolve the connectivity between neighbors. This
routine accepts an optional argument with the connectivity, which can for
example be used to specify periodic boundary conditions. An example is shown
below:

```{f90}
    integer :: nb_list(a$D_num_neighbors, n_boxes)

    nb_list(:, :) = af_no_box      ! Start with default value
    nb_list(a2_neighb_lowx, 1) = 2 ! Box 1's lower x-neighbor is box 2
    nb_list(a2_neighb_lowy, 1) = 1 ! Box 1 is periodic in y-direction
    nb_list(a2_neighb_lowy, 2) = 2 ! Box 2 is periodic in y-direction

    call a2_set_base(tree, n_boxes, ix_list, nb_list)
```

By specifying the default value `af_no_box`, Afivo will automatically handle the
internal connectivity between the boxes. Note also that the periodic boundary
only needs to be specified from one side, the other sides are handled
automatically.

## Examples

* @ref computational_domain_2d.f90
* @ref computational_domain_3d.f90
