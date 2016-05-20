\section Contents
<ul>
	<li> \ref sect_intro </li>
	<li> \ref sect_example1
		<ul>
			<li> \ref sect_parallel </li>
		</ul>
	</li>

</ul>

\section sect_intro Introduction 

In this tutorial we describe several examples.

\subsection sect_example1 Example1: Afivo for solving the Poisson equation in two- and three dimensions

In this section we present two Poisson test problems
to demonstrate the multigrid behavior on a partially refined mesh.
We use the **method of manufactured solutions**: from an analytic solution the right-hand side and boundary
conditions are computed. Two test problems are considered, a
constant-coefficient
<a class="code" href="poisson__basic__2d_8f90.html#aef02b53cac21b72a47afc1c34286c443">two-</a> and
<a class="code" href="poisson__basic__3d_8f90.html#a1b0ee6adabd66bb15df1120e11776ead">three-dimensional</a>
Poisson equation shown in <a href="examples.html">Examples</a>
<a name="eq_mg-example-lpl-1" >(11)</a>
\f[
  \nabla^2 u = \nabla \cdot (\nabla u) = \rho ,
\f]
on a two-dimensional rectangular domain \f$[0,1] \times [0,1]\f$ and
on a three-dimensional cubic domain \f$[0,1] \times [0,1] \times [0,1]\f$, respectively.
We pick the following solution for \f$u\f$
\f[
  u(r) = \exp(|{\vec{r}-\vec{r}_1}|/\sigma) + \exp(|{\vec{r}-\vec{r}_2}|/\sigma),
\f]
where \f$\vec{r_1} = (0.25, 0.25)\f$, \f$\vec{r_2} = (0.75, 0.75)\f$, in the 2D case
  and \f$\vec{r_1} = (0.25, 0.25, 0.25)\f$, \f$\vec{r_2} = (0.75, 0.75, 0.75)\f$, in the 3D case
  and \f$\sigma = 0.04\f$.
An analytic expression for the right-hand side \f$\rho\f$ is obtained by plugging
the solution in <a href="#eq_mg-example-lpl-1">equation(11)</a>.

In this tutorial we concentrate on the 2D case. The 3D example will be used for timing results.
<code>gaussian</code> of type <a class="el" href="structm__gaussians_1_1gauss__t.html">gauss_t</a>:

\verbatim
  !> A type to store a collection of gaussians in
  type gauss_t
     integer :: n_gauss                !< Number of gaussians
     integer :: n_dim                  !< Dimensionality
     real(dp), allocatable :: ampl(:)  !< Amplitudes
     real(dp), allocatable :: sigma(:) !< Widths
     real(dp), allocatable :: r0(:,:)  !< Centers
  end type gauss_t
\endverbatim
defined in module <code>m_gaussians</code>. This module also contains the routine 
<a class="el" href="namespacem__gaussians.html#a050ea2228d5aa753d13e1f81996045d7">gauss_init</a>,
which is called to store variable <code>gaussian</code>:

\verbatim
! The manufactured solution exists Gaussians, which are stored in gaussian
  if (n_gaussian == 2) then
     ! Amplitudes:  [1.0_dp, 1.0_dp]
     ! Sigmas    :  [0.04_dp, 0.04_dp]
     ! Locations :  [[0.25_dp, 0.25_dp], [0.75_dp, 0.75_dp]]
     call gauss_init(gaussian, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
          reshape([0.25_dp, 0.25_dp, 0.75_dp, 0.75_dp], [2,2]))
  end if
\endverbatim

Defining <code>box_size = 8</code>  denotes that each box has \f$ 8 \times 8\f$ cells, and,
since the problem is solved on a rectangular domain of \f$[0,1] \times [0,1]\f$,
implies that variable <code>dr = 1 / 8 </code>.
The variable <code>tree</code> of type
<a class="el" href="structm__a2__types_1_1a2__t.html">a2_t</a>
is initialized by a call to 
<a class="el" href="namespacem__a2__core.html#aa3e7687c41f7b3915de71b96f4511de3">a2_init</a>

\verbatim
! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! Number of cells per coordinate in a box
       n_var_cell, &   ! Number of face-centered variables
       n_var_face, &   ! Number of cell-centered variables
       dr)             ! Distance between cells on base level
\endverbatim
Note that <a class="el" href="namespacem__a2__core.html#aa3e7687c41f7b3915de71b96f4511de3">a2_init</a>
also has some optional variables. In example 
<a class="code" href="poisson__basic__2d_8f90.html#aef02b53cac21b72a47afc1c34286c443">poisson_basic_2d</a>
the variable <code>n_var_cell</code>, the number of cell centered values, has been set
to 4 (reserved for <code>i_phi, i_rhs, i_err</code> and <code>i_tmp</code>, see below),
whereas the variable <code>n_var_face</code>, the number of face centered values,
has been initialized to 0.

In Afivo the coarsest mesh, which covers the full computational domain, is not
supposed to change. To create this mesh there is a routine
<a class="el" href="namespacem__a2__core.html#ab7007734c1625a6057a63f4676ce7244">a2_set_base</a>,
which takes as input the spatial indices of the coarse boxes and their neighbors.
Below, a 2D example is shown for creating a single coarse box at index \f$(1,1)\f$. 
Physical (non-periodic) boundaries are indicated by a
negative index for the neighbor. By adjusting the neighbors one can specify
different geometries, the possibilities include meshes that contain a hole, or
meshes that consist of two isolated parts. The treatment of boundary conditions
is discussed in section \ref sect_fill-ghost-cell.

\verbatim
  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of box 1

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1            ! Dirichlet zero -> -1

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)
\endverbatim

The routine
<a class="el" href="namespacem__a2__types.html#a42dbfb05abcda7baa560a475849d2bbb">a2_print_info</a>
shows some values of the tree. At this position in the example, it lists:
\verbatim
 *** a2_print_info(tree) ***
 Current maximum level:           1
 Maximum allowed level:          30
 Number of boxes used:            3
 Memory limit(boxes):       4793490
 Memory limit(GByte):      0.16E+02
 Highest id in box list:          3
 Size of box list:             1000
 Box size (cells):                8
 Number of cc variables:          4
 Number of fc variables:          0
 Type of coordinates:             0
 min. coords:          0.0000E+00  0.0000E+00
 dx at min/max level   0.1250E+00  0.1250E+00
\endverbatim

Note that the number of boxes is 3, although the maximum multigrid level  is just 1. 
This arises from the fact that for multigrid cycling also a level 0 of box size \f$ 4 \times 4 \f$,
and a level -1 have been defined of size  \f$ 2 \times 2 \f$.
The memory limit has been set on 16 GB (default value) which results in at most 4793490 boxes, based on the box size and the number of ghost cells in each direction (i.e., 2) and
the numbers of cell and face centered variables. 
The length of the box list is 1000, but if more boxes are required, the length will be extended
and boxes which are not longer used will be cleaned up.
Since <code>box_size = 8</code>  this implies that each box has \f$ 8 \times 8\f$ cells, excluding
ghost cells, and, since the problem is solved on a rectangular domain of \f$[0,1] \times [0,1]\f$,
(<code> min. coords:          0.0000E+00  0.0000E+00</code>, the cell size <code>dx</code>
at the coarsest level is <code>dr = 1 / 8 </code>.

Next, we start refining the mesh by means of subroutine
<td class="memItemRight" valign="bottom"><a class="el" href="poisson__basic__2d_8f90.html#a1c44f0de37167b4b60895cf38c39e430">refine_routine</a>
which has been added to program 
<a class="code" href="poisson__basic__2d_8f90.html#aef02b53cac21b72a47afc1c34286c443">poisson_basic_2d</a>

\verbatim
   ! Return the refinement flag for boxes(id)
     subroutine refine_routine(boxes, id, refine_flag)
       type(box2_t), intent(in) :: boxes(:)
       integer, intent(in)      :: id
       integer, intent(inout)   :: refine_flag
       integer                  :: i, j, nc
       real(dp)                 :: xy(2), dr2, drhs
   
       nc = boxes(id)%n_cell
       dr2 = boxes(id)%dr**2
   
       outer: do j = nc
          do i = nc
             xy = a2_r_cc(boxes(id), [i, j])
   
             ! This is an estimate of the truncation error in the right-hand side,
             ! which is related to the fourth derivative of the solution.
             drhs = dr2 * gauss_4th(gaussian, xy) / 
   
             if (abs(drhs) > 0.05_dp) then
                refine_flag = af_do_ref
                exit outer
             end if
          end do
       end do outer
     end subroutine refine_routine
\endverbatim
Here variable <code>id</code> denotes the index number of the box.
<a class="el" href="namespacem__a2__types.html#a36d16420768ccb9874afb2cfe392ed01">a2_r_cc</a>
is one of those useful features from module <code>m_a2_types</code> to calculate the center of a cell.
The refinement is based on the fourth derivative of the solution.
This example only shows when a box should be refined or not, the derefinement does not play a role here.
Besides, value <code>af_do_ref</code>, indicating you want to refine a box,
the value <code>af_keep_ref</code>, indicating you want to keep a box's refinement, and
<code>af_rm_ref</code>, indicating you want to derefine a box, are available.
The routine  
<a class="el" href="namespacem__a2__core.html#a3f0dd0f48f710e79b6b2085177b19dd6">a2_adjust_refinement</a>
ensures the 2:1 balance is ensured, so that there is never a jump of more than
one refinement level between neighboring boxes.
These constraints are automatically handled, so that the user-defined refinement
function does not need to impose them.

When boxes are added or removed in the refinement procedure, their connectivity
is automatically updated.
References to a removed box are removed from its parent and neighbors.
When a new box is added, its neighbors are found through its parent.
Three scenarios can occur: the neighbor can be one of the other children of the
parent, the neighbor can be a child from the neighbor of the parent, or the
neighbor does not exist.
In the latter case, there is a refinement boundary, which is indicated by the
special value
<a class="el" href="namespacem__af__t.html#aa0395be4a37c035fae9b5663b2d41698">af_no_box</a>.

In the next loop, first the routine 
<a class="el" href="namespacem__a2__utils.html#a7209ee5c6acc86eb86c66a9198f135df">a2_loop_box</a>
is called for each box, with as argument the user-defined routine
\verbatim
  ! This routine sets the initial conditions for each box
     subroutine set_initial_condition(box)
       type(box2_t), intent(inout) :: box
       integer                     :: i, j, nc
       real(dp)                    :: xy(2)
   
       nc = box%n_cell
   
       do j = nc
          do i = nc
             ! Get the coordinate of the cell center at i,j
             xy = a2_r_cc(box, [i,j])
   
             ! And set the rhs values
             box%cc(i, j, i_rhs) = gauss_laplacian(gaussian, xy)
          end do
       end do
     end subroutine set_initial_condition
\endverbatim
which calls  <code>gauss_laplacian</code> from <code>m_gaussian</code> corresponding the problem
described above. The routine <code>a2_r_cc</code> computes the cell center of the cells in the box.
In this example, each box has 4 cell centered matrices. Here <code>box%cc(:,:,i_rhs)</code> is
initialized
with the right hand side values. 

The following <code>do loop</code> all right hand side field of the boxed used are initialized and
boxes are refined in accordance with the refinement routine to obtain an adaptive mesh.
\verbatim
  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_initial_condition)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a2_adjust_refinement(tree, refine_routine, refine_info)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
\endverbatim
If no further refinement is required, i.e., <code>refine_info\%n_add == 0</code>, the do loop stops.

The routine
<a class="el" href="namespacem__a2__types.html#a42dbfb05abcda7baa560a475849d2bbb">a2_print_info</a>
shows after the refinement loop:
\verbatim
*** a2_print_info(tree) ***
 Current maximum level:          10
 Maximum allowed level:          30
 Number of boxes used:         5487
 Memory limit(boxes):       4793490
 Memory limit(GByte):      0.16E+02
 Highest id in box list:       5487
 Size of box list:             8030
 Box size (cells):                8
 Number of cc variables:          4
 Number of fc variables:          0
 Type of coordinates:             0
 min. coords:          0.0000E+00  0.0000E+00
 dx at min/max level   0.1250E+00  0.2441E-03
\endverbatim
So ten levels of refinement are required leading to 5487 boxes and a cell size of the highest level
of \f$1 / 4096\f$. 

By means of a call to <code>a2_write_vtk</code> from module <code>m_a2_output</code>
a <code>vtu</code> file can be produced to show the adapted mesh, but also the solution.
\verbatim
  ! This writes a VTK output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  write(fname, "(A,I0)") "poisson_basic_2d_", 0
  call a2_write_vtk(tree, trim(fname), dir="output")
\endverbatim
The resulting file "poisson_basic_2d_0.vtu" in directory output can be viewed by e.g., visit.
See Figure <a href="#fig_start_poisson_basic_2d">Figure 1</a>.

<a name="fig_start_poisson_basic_2d" />
a) | b) 
------------- | -------------
<img src="../../figures/poisson_basic_2d_rhs.png" width=400px /> | <img src="../../figures/poisson_basic_2d_0.png" width=400px />

**Figure 1**. a) initialization of the right hand side b) adapted mesh

Now the grid has been constructed to match the RHS function and in accordance with the refinement function
<code>refine_routine</code> showed above, we continue to solve the Poisson problem using multigrid.
Therefore we first initialize the multigrid options.
The routine <a class="el" href="structm__a2__multigrid_1_1mg2__t.html">mg2_t</a>
performs some basic checks and sets default values where necessary in variable <code>mg</code> of type 
<a class="el" href="structm__a2__multigrid_1_1mg2__t.html">mg2_t</a>.

\verbatim
  !> Type to store multigrid options in
  type, public :: mg2_t
     integer :: i_phi        = -1 !< Variable holding solution
     integer :: i_rhs        = -1 !< Variable holding right-hand side
     integer :: i_tmp        = -1 !< Internal variable (holding prev. solution)

     integer :: i_eps        = -1 !< Optional variable (diel. permittivity)
     integer :: i_lsf        = -1 !< Optional variable for level set function
     integer :: i_bval       = -1 !< Optional variable for boundary value

     integer :: n_cycle_down = -1 !< Number of relaxation cycles in downward sweep
     integer :: n_cycle_up   = -1 !< Number of relaxation cycles in upward sweep
     integer :: n_cycle_base = -1 !< Number of relaxation cycles at bottom level

     logical :: initialized  = .false.

     !> Routine to call for filling ghost cells near physical boundaries
     procedure(a2_subr_bc), pointer, nopass   :: sides_bc => null()

     !> Routine to call for filling ghost cells near refinement boundaries
     procedure(a2_subr_rb), pointer, nopass   :: sides_rb => null()

     !> Subroutine that performs the (non)linear operator
     procedure(mg2_box_op), pointer, nopass   :: box_op => null()

     !> Subroutine that performs Gauss-Seidel relaxation on a box
     procedure(mg2_box_gsrb), pointer, nopass :: box_gsrb => null()

     !> Subroutine that corrects the children of a box
     procedure(mg2_box_corr), pointer, nopass :: box_corr => null()

     !> Subroutine for restriction
     procedure(mg2_box_rstr), pointer, nopass :: box_rstr => null()
  end type mg2_t
\endverbatim
To this end, we perform a number of multi-grid iterations. Within each loop we call routine
<a class="el" href="namespacem__a2__multigrid.html#a9a06b7943954767e5a207ed1c1a1d5bb">mg2_fas_fmg</a>, performing a full approximation scheme using full multigrid.
At the first call of
<a class="el" href="namespacem__a2__multigrid.html#a9a06b7943954767e5a207ed1c1a1d5bb">mg2_fas_fmg</a>sets the initial guess for <code>phi</code> and
restricts the right hand side fields <code>rhs</code> from the highest level down to the
coarsest level.
After the first call of
<a class="el" href="namespacem__a2__multigrid.html#a9a06b7943954767e5a207ed1c1a1d5bb">mg2_fas_fmg</a>
it sets the right hand side field on coarser grids and
restricts the field with the holding solution of <code>phi</code>.
If required the ghost cells are filled.

Next the user-defined routine <code>set_error</code> is called to computed the error
compared to the analytic solution.
\verbatim
 ! Set the error compared to the analytic solution
  subroutine set_error(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - gauss_value(gaussian, xy)
       end do
    end do
  end subroutine set_error
\endverbatim

The results of the first 10 iterations steps are listed below:
\verbatim
 Multigrid iteration | max residual | max error
       1                1.00579E+01   5.13811E-05
       2                4.54153E-01   5.13903E-05
       3                2.26078E-02   5.13902E-05
       4                1.14302E-03   5.13902E-05
       5                5.79126E-05   5.13902E-05
       6                2.94218E-06   5.13902E-05
       7                1.50202E-07   5.13902E-05
       8                2.38767E-08   5.13902E-05
       9                2.36430E-08   5.13902E-05
      10                2.01674E-08   5.13902E-05
\endverbatim
Note that the residual is still decreasing, whereas the maximum error reaches
its minimum after about three iterations.

After 10 multigrid iterations we obtain 
the resulting file "poisson_basic_2d_10.vtu" in directory output, which can be viewed by e.g., visit.
See <a href="#fig_end_poisson_basic_2d">Figure 3</a>.

<a name="fig_end_poisson_basic_2d" />
a) | b)
------------- | -------------
<img src="../../figures/poisson_basic_2d_10_phi.png" width=400px /> | <img src="../../figures/poisson_basic_2d_10_err.png" width=400px />

**Figure 3**. a) solution \f$ \phi \f$ b) error \f$ \epsilon \f$ after 10 multigrid iteration steps.


\subsection sect_parallel Parallel results
Most operations in Afivo loop over a number of boxes, for example the leaves at a certain refinement
level, as described in the Afivo Manual. All such loops have been parallelized by adding OpenMP
statements around them, for example as in Figure 2.

\verbatim
   do lvl = 1, tree%max_lvl
      !$omp parallel do private(id)
      do i = 1, size(tree%lvls(lvl)%ids)
         id = tree%lvls(lvl)%ids(i)
         call my_method(tree%boxes(id))
      end do
      !$omp end parallel do
   end do
\endverbatim
Figure 2. Fortran code fragment that shows how to call my_method for all the boxes in a tree,
from level 1 to the maximum level. Within each level, the routine is called in parallel using OpenMP.

The parallel speedup that one can get depends on the cost of the algorithm that one is using.
The communication cost (updating ghost cells) is always about the same, so that an expensive
algorithm will show a better speedup. Furthermore, on a shared memory system, it is not unlikely
for an algorithm to be memory-bound instead of CPU-bound.

The info below is coming from the corresponding 3d problem for boxes of size \f$8 \times 8 \times 8\f$, excluding ghost cells.
\verbatim
 *** a3_print_info(tree) ***
 Current maximum level:           9
 Maximum allowed level:          30
 Number of boxes used:         8811
 Memory limit(boxes):        525314
 Memory limit(GByte):      0.16E+02
 Highest id in box list:       8811
 Size of box list:            12502
 Box size (cells):                8
 Number of cc variables:          4
 Number of fc variables:          0
 Type of coordinates:             0
 min. coords:          0.0000E+00  0.0000E+00  0.0000E+00
 dx at min/max level   0.1250E+00  0.4883E-03
\endverbatim
On the coarsest level the variable \f$dx = 1 / 8\f$, the size of a cell, whereas on the
highest level \f$dx = 1 / 2048\f$.
We have compared, for this 3D case, the wall clock time, obtained by a calculation with just one
thread and with multiple threads. This shows that the use of OpenMP is very succesful.
We distinguish the wall clock time for both the creation of the adaptive mesh refinement part
of the program and the multigrid process.
For the  adaptive mesh refinement part, we obtain a factor of 5.3, comparing the wall clock time
for 16 and just one thread. 
For the multigrid process we achieve a speedup factor of nearly a factor 9.
For a multigrid process, which is hard to parallelize, this is a very good result.

<img src="../../figures/WCT_basic_3d_lisa_new.png" width=300px /> 
