\section sect_intro Introduction 

\section sect_photoi Photoionization (11.3.3)

Photoionization can play an important role in electrical discharges in air and
other gases. 

\verbatim
[ .. ]
\endverbatim

If we consider models for streamer discharges, then two typical challenges are
solving Poisson’s equation to obtain the electrostatic potential, and computing
the photoionization profile. Both challenges involve a non-local process or
interaction – at least when the speed of light is assumed to be effectively
infinite.
However, the required accuracy for these two problems is quite different:
discharges are much less sensitive t  the amount of photoionization around them
than to the electric field. The method presented here was thus designed to
be fast and flexible, but not necessarily highly accurate.

\verbatim
[ .. ]
\endverbatim

\subsection sect_descr Description of the method

In our Monte Carlo method for photoionization, we assume that photon scattering
can be neglected, and that the direction of ionizing photons is isotropic. The
method can be divided in three parts:
	-  Determine the coordinates at which photons are produced, and
store these coordinates in a list \f$L_{src}\f$.
	- Determine the coordinates at which photons are absorbed, and
store these coordinates in a list \f$L_{dst}\f$.
	- Compute the resulting photoionization profile on a mesh.

The implementation of these steps is described below.

\subsection sect_ionizing The source of ionizing photons

In a plasma fluid simulation, the computational domain is divided into cells. We
assume that for each cell the production of ionizing photons \f$I\f$ is known within
a given time step \f$\Delta t\f$. In our Monte Carlo method, this information is converted
to a list  \f$L_{src}\f$ of approximately \f$N\f$ discrete photons in the following way.
1. Determine the total photon production \f$I_{\Sigma}\f$ within the time step \f$\Delta t\f$,
by summing over all the cells.
2. For each cell \f$n_{\gamma} = N I / I_{\Sigma}\f$ photons should be produced.
To convert \f$n_{\gamma}\f$ to an integer, first draw a uniform random number
\f$U (0 , 1)\f$.
If \f$n_{\gamma} - \lfloor{ n_{\gamma}} \rfloor > U\f$ round up,
else round down.
3. For each produced photon, add the coordinate of the cell center to the list
\f$L_{src}\f$ of photons.

Some additional remarks: In principle the number of photons N can be chosen
adaptively, for example to create discrete ‘super-photons’ with a given weight.
It is also possible to assign different weights to different photons by modifying
the above procedure, see section 11.3.4.
Instead of rounding \f$n_{\gamma}\f$ to an integer as described above, it can be more
realistic to sample from the Poisson distribution with parameter \f$n_{\gamma}\f$.
For \f$n_{\gamma} \ll 1\f$ the result would be almost identical, but for larger
values there are differences.
However, most of the random fluctuations in our method are due to the stochastic
photon direction and absorption, as discussed in the next section. Therefore our
rounding with uniform random numbers should be fine for most applications.
When grid cells are large compared to typical photoionization length scales,
one could determine a ‘subgrid’ production location in the cell, instead of using
the cell center. However, this scenario should generally be avoided, because the
resulting photoionization profile would itself be insufficiently resolved.

\subsection sect_absorption Absorption of ionizing photons

Now that we know where photons are produced, we need to determine where
they are absorbed. This is done in the following way. First, we determine
the absorption distance \f$r\f$ for a photon. Given an absorption function \f$f(r)\f$,
the cumulative absorption function \f$F(r) =\int_0^r f(r')dr'\f$ can be computed, either
numerically or analytically. Then a so-called lookup table is constructed with
columns \f$F(r)\f$ and \f$r\f$. Now we can do inversion sampling: given a random number
\f$U(0, 1)\f$, the corresponding distance \f$r\f$ is obtained by linear interpolation from the
table. The procedure is illustrated in figure 11.1. (In the special case where the
inverse of \f$F(r)\f$ is known, one could directly compute
\f$r = F^{-1}(U)\f$, but this will often be slower than using a lookup table.)
Then a random orientation \f$\hat{r}\f$ for the photon is determined, using a procedure
for picking a point on a sphere proposed by Marsaglia [205]:
1. Get two random numbers \f$U_1(-1,1)\f$ and \f$U_2(-1,1)\f$. If
\f$U_1^2 + U_2^2 \leq 1\f$ accept them, otherwise draw new numbers.
2. Let \f$a = U_1^2 + U_2^2\f$ and \f$b = 2 \sqrt{1-a}\f$.
The isotropic random direction is
\f$\hat{r} = (b U_1, b U_2, 1 - 2a)\f$ in Cartesian coordinates.
In the special case of a 2D Cartesian \f$(x, y)\f$ coordinate system, we should not
pick points on a sphere but points on a circle. Then the second step is replaced
by
\f[
\hat{r} = ( \frac{U_1^2-U_2^2}{U_1^2+U_2^2},\frac{2 U_1 U_2}{U_1^2+U_2^2}) .
\f]
Given the direction and the distance, the location of absorption
\f$\mathbf{r} = r \hat{\mathbf{r}}\f$ is known,
which is added to the list \f$L_{dst}\f$.
Sometimes, the typical absorption length of photons is much larger than the
domain size. Since most photons will not be absorbed within the domain, the
Monte Carlo approximation becomes less accurate. This can be resolved by

[nog niet klaar]

\subsection sect_config Usage of configuration files

The streamer programs discussed here are using configuration files
and set variables accordingly. For each simulation program a configuration
file is available in directory
<a class="el" href="dir_53bba52d758a816fa6b9c79aab7dc56b.html">ConfigureFiles</a>.
After compiling the streamer programs the executables can be run as follows:
* streamer_2d ConfigureFiles/streamer_2d.cfg
* streamer_3d ConfigureFiles/streamer_3d.cfg
* streamer_cyl ConfigureFiles/streamer_cyl.cfg

In <code>streamer_2d.cfg</code> values like the name of the simulation,
the directory where the output should be written, the kind of gas,
the end time of the simulation are listed. Values corresponding to the
default values used by <code>streamer_2d.f90</code> as initialized by calling
subroutine
<a class="el" href="namespacem__streamer.html#ae57f8b89a83407090092d2ca5e0133af">st_create_config</a>
do not have to be included.
After the configuration file has been read, a new configuration file will be
created and written in directory 'output' (this directory should be present!!),
named 'output/streamer_2d_output.cfg'.
Besides the variables initialized by reading <code>streamer_2d.cfg</code>
it contains default values for the variables. 
Analogously, the configuration files <code>streamer_3d.cfg</code> and
<code>streamer_cyl.cfg</code> are read and new files
'output/streamer_3d_output.cfg' and 'output/streamer_cyl_output.cfg'
are created and written.
To be able to compare the original configuration files *.cfg the new created
'_output' files have been added to the directory 
<a class="el" href="dir_53bba52d758a816fa6b9c79aab7dc56b.html">ConfigureFiles</a>.

\subsection sect_lookup Lookup tables (LT)

This streamer package makes use of lookup tables. These tables can be used to
efficiently interpolate one or more values.
Therefore a special lookup table type <code>lookup_table_t</code> has been designed
\verbatim
  !> The lookup table type. There can be one or more columns, for which values
  !> can be looked up for a given 'x-coordinate'.
  type lookup_table_t
     private
     integer  :: n_rows !< The number of rows
     integer  :: n_cols !< The number of columns
     real(dp) :: x_min  !< The minimum lookup coordinate
     real(dp) :: dx     !< The x-spacing in the lookup coordinate
     real(dp) :: inv_dx !< The inverse x-spacing

     !The table is stored in two ways, to speed up different types of lookups.
     real(dp), allocatable :: cols_rows(:, :) !< The table in column-major order
     real(dp), allocatable :: rows_cols(:, :) !< The table in row-major order
  end type lookup_table_t
\endverbatim
Besides the type <code>lookup_table_t</code> another special type <code>LT_loc_t</code> has been created,
which can be used to speed up multiple lookups of different columns
\verbatim
  type LT_loc_t
     private
     integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
     real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                          !< as a real number between 0 and 1.
  end type LT_loc_t
\endverbatim
In order to optimize the use of the lookup tables, the following public functions and
subroutines are presented in module <code>m_lookup_table</code>:
\verbatim
   ! Returns a new lookup table
   LT_create
   ! Returns the x-coordinates of the lookup table
   LT_get_xdata
   ! Linearly interpolate the (x, y) input data to the new_x coordinates
   LT_get_spaced_data
   ! Fill the column with index col_ix using the linearly interpolated (x, y) data
   LT_set_col
   ! Add a new column by linearly interpolating the (x, y) data
   LT_add_to_col
   ! Add the (x,y) data to a given column
   LT_add_col
   ! Get a location in the lookup table
   LT_get_loc
   ! Get the values of all columns at x
   LT_get_col
   ! Get the value of a single column at x
   LT_get_mcol
   ! Get the values of all columns at a location
   LT_get_col_at_loc
   ! Get the value of a single column at a location
   LT_get_mcol_at_loc
   ! Return the number of rows
   LT_get_num_rows
   ! Return the number of columns
   LT_get_num_cols
   ! Get the x-coordinates and the columns of the lookup table
   LT_get_data
   ! Compute by use of linear interpolation the value in the middle of
   ! domain D = [x_list(1) , x_list(size(x_list))].
   ! If x_value is left of domain  D, then the value becomes the value
   ! at the left side of D,
   ! if x_value is right of domain D, then the value becomes the value
   ! at the right side of D
   LT_lin_interp_list
   ! Write the lookup table to file (in binary, potentially unportable)
   LT_to_file
   ! Read the lookup table from file (in binary, potentially unportable)
   LT_from_file
\endverbatim

\subsection sect_TD Transport data (TD)

The module <code>m_transport_data</code> provides a routine <code>TD_get_td_from_file</code>
for reading in arbitrary transport data. <code>TD_get_td_from_file</code> reads in transport
data from a file. Searches 'file_name' for transport 'data_name' concerning 'gas_name'
\verbatim
! Look for collision processes with the correct gas name in the file,
! which should contains entries like below:

!     Efield[V/m]_vs_energy[eV]     [description of the type of transport data]
!     AIR                           [the gas (mixture) name]
!     COMMENT: Compiled by xxxxx    [possibly comments, these are optional]
!     UPDATED: 2010-06-24 15:04:36
!     ------------------            [at least 5 dashes]
!     xxx       xxx                 [transport data in two column format]
!     ...       ...
!     xxx       xxx
!     ------------------
\endverbatim
The directory 'input' contains five transport data files for cross sections:
\verbatim
   td_air_siglo_swarm.txt   
   td_air_phelps_bolsig.txt 
   td_air_props.txt         
   td_n2_siglo_133mbar.txt  
   td_example.txt      
\endverbatim
The electron cross sections were retrieved from
<a href="www.lxcat.net/SIGLO">SIGLO(N2,O2)</a> and
<a href="http://jilawww.colorado.edu/~avp/">Phelps(NO)</a> databases and
<a href="http://fr.lxcat.net/solvers/BOLSIG+/">BOLSIG</a>.
The main utility of the databases is to obtain electron transport coefficients and collision rate coefficients
from more fundamental cross section data, which can then be used as input for fluid models.
In the configuration file (see above) as shown in 
<a class="el" href="dir_53bba52d758a816fa6b9c79aab7dc56b.html">ConfigureFiles</a>
one can specify which transport data file one prefers, e.g.,
<code>input_file = input/td_air_siglo_swarm.txt</code>.

\subsection sect_LL Methods to work properly with linked lists (LL)

Two different types for linked lists are defined:
\verbatim
   type LL_int_t
      private
      integer                 :: x
      type(LL_int_t), pointer :: next
   end type LL_int_t

   type LL_int_head_t
      private
      type(LL_int_t), pointer :: ptr => null()
      integer                 :: n_values = 0
   end type LL_int_head_t
\endverbatim
The module <code>m_linked_list</code> contains four subroutines to work properly with linked lists:
\verbatim
   LL_Clear    ! Clears an linked list of type LL_int_head_t
   LL_add      ! Adds an element to a linked list of type LL_int_head_t at position head
   LL_pop      ! TODO
   LL_get_size ! Computes the number of elements of a linked list of type LL_int_head_t
\endverbatim

\subsection sect_MRGRNK Merge-sort ranking

The subroutines from this subsection are not used elsewhere by the package

\subsection sect_FI Finding special indices in monotonically increasing sorted lists (FI)

The module <code>m_find_index</code> contains several auxiliary subroutines for finding special indices in
monotonically increasing sorted lists. 

\subsection sect_GM Geometric operations and calculations (GM)

This subsection provides routines for geometric operations and calculations, e.g.,
subroutine <code>GM_density_line</code> for computing .......

\subsection sect_UC Physical constants (UC)

This subsection contains several physical constants:
\verbatim
  UC_eps0             = 8.8541878176d-12        ! permitivity of vacuum (SI)
  UC_elem_charge      = 1.6022d-19              ! the elementary charge in Coulombs
  UC_elec_charge      = -1.6022d-19             ! the electron charge in Coulombs
  UC_elec_volt        = 1.6022d-19              ! the eV in joules
  UC_elec_mass        = 9.10938189d-31          ! the electron mass in kg
  UC_atomic_mass      = 1.66053886d-27          ! the atomic mass unit in kg
  UC_N2_mass          = 28.0D0 * UC_atomic_mass ! The mass of a N2 molecule
  UC_lightspeed       = 299792458d0             ! the speed of light in m/s
  UC_boltzmann_const  = 1.3806503d-23           ! the Boltzmann constant
  UC_bohr_radius      = 5.29d-11                ! the Bohr radius (m)
  UC_torr_to_bar      = 133.322368 * 1.0d-5     ! one Torr in units of bar
  UC_elec_q_over_eps0 = UC_elec_charge / UC_eps0
  UC_elec_q_over_m    = UC_elec_charge / UC_elec_mass
\endverbatim
All constants are of type <code>real(dp)</code>, where <code>dp = kind(0.0d0)</code>.


A Lookup table with the absorption locations of the photons has been created.
These coordinates need to be mapped to a mesh to get the photoionization profile
Two options can be considered: the photons can be absorbed on a single mesh with a constant grid spacing, or they can be absorbed at different levels in a refined mesh. 
In both cases the nearest grid point (NGP) is used to map the photons to densities,
which means that photons are absorbed within a single cell.
With linear (and higher order) interpolation the resulting density profiles are smoother,
but it becomes harder to handle refinement boundaries.

\subsection sect_similarities Similarities of streamer examples

Since the three streamer examples, <code>streamer_2d.f90</code>,
<code>streamer_3d.f90</code>, and <code>streamer_cyl.f90</code> show great similarities, 
a module called <code>m_streamer.f90</code> contains a couple of subroutines,
like the subroutines for handling the configuration files,
getting the electric field, setting the voltages at a given time or to a fixed value,
a new implementation of the Koren limiter, 
but also many variables like the names of centered variables, photoionization
quantities, time steps and so on.

\section sect_2d streamer_2d: a 2-dimensional streamer program

As described above, the program starts with creating, loading and writing configuration files, such that the variables are initialized in a proper and clear way.
Next a lookup table will be created of type 
<a class="el" href="structm__lookup__table_1_1lookup__table__t.html">lookup_table_t</a>
\verbatim
type lookup_table_t
   private
   integer               :: n_rows
   integer               :: n_cols
   real(dp)              :: x_min, inv_dx
   real(dp), allocatable :: cols_rows(:, :)
   real(dp), allocatable :: rows_cols(:, :)
   end type lookup_table_t
\endverbatim
to initialize the transport coefficients. See subroutine
<a class="el" href="namespacem__streamer.html#ac544e682ffced5108a9e5f33b7c6c1ba">st_load_transport_data</a>. 

The initial conditions will be stored in a variable of type
<a class="el" href="structm__streamer_1_1initcnd__t.html">initcnd_t</a>
\verbatim
! Type to store initial conditions in
   type initcnd_t
      real(dp)              :: background_density
      integer               :: n_cond
      real(dp), allocatable :: seed_r0(:, :)
      real(dp), allocatable :: seed_r1(:, :)
      real(dp), allocatable :: seed_density(:)
      real(dp), allocatable :: seed_width(:)
      integer, allocatable  :: seed_falloff(:)
   end type initcnd_t
\endverbatim

Next the tree, from the Afivo package, will be initialized. The tree contains
all the mesh information. For more details how to use Afivo,
see the tutorial?? of Afivo.

For a streamer model the electric field must be computed on a tree.
This can be performed by
<a class="el" href="streamer__2d_8f90.html#a3078ef0e5f2dbb937b8b9d54baccaa6c">compute_electric_field</a>,
a user defined subroutine, which uses multigrid to compute the electric field
for each level of the mesh.
The mesh is sufficiently refined when no more new grids are created based on criteria
defined in the refinement subroutine
<a class="el" href="streamer__2d_8f90.html#a6b6ddafe30bfd36b772e2bc7888bf5fb">refine_routine</a>,
defined by the user, too.

Further the user may decide whether photoionization will be added.
Also here the user has to enclose a subroutine how the photoionization
is defined.
The photon production rate per cell is proportional to the ionization rate
and will be computed per cell and stored in the <code>i_photo</code> part
for each box.
For more details see
<a class="el" href="streamer__2d_8f90.html#abc619b5d02ebb950749f874f4f9f6dcb">set_photoionization</a>
as an example. 

\subsection sect_timestepping Time Stepping

The maximum allowed time step is restricted by e.g., the CFL criteria:
\verbatim
   ! CFL condition
   dt_cfl = dr_min / (mobility * max_electric_fld) ! Factor ~ sqrt(0.5)
  
   ! Diffusion condition
   dt_dif = 0.25_dp * dr_min**2 / diff_coeff
  
   ! Dielectric relaxation time
   dt_drt = uc_eps0 / (uc_elem_charge * max_mobility * max_dns)
  
   ! Ionization limit
   dt_alpha = 1 / max(mobility * max_electric_fld * alpha, epsilon(1.0_dp))
\endverbatim
where the minimal cell distance <code>dr_min</code> can be computed from
<code>tree</code>,
the values of the variables <code>mobility</code>, and
<code>max_mobility</code> can be found in the lookup table,
see \ref sect_lookup, the diffusion coefficient <code>diff_coeff</code>
as well as <code>alpha</code> are also included in that table. 
The parameters
<code>uc_eps0</code> and <code>uc_elem_charge</code> have the following values:

	* <code>uc_eps0 = 8.8541878176d-12</code> and
	* <code>uc_elem_charge = 1.6022d-19</code>.

This leads to a maximum time step \f$dt_{amr-max}\f$ :
\f[
    dt_{amr-max} = \frac{1}{2} \; min(\frac{1}{(\frac{1}{dt_{cfl}} + \frac{1}{dt_{dif}})},
    \; dt_{alpha}, \; dt_{max}) 
\f]
where \f$ dt_{max}\f$ arises from the configuration file <code>streamer_2d.cfg</code>.
The maximum time step is calculated by the user defined function
<a class="el" href="streamer__2d_8f90.html#a2eb820779f8ab73aa61e8665e0c13473">get_max_dt</a>. 

Inside the time stepping loop, two forward Euler iteration steps are performed.
In an Euler iteration the following steps are completed:

	* First calculate the fluxes
		* In order to obtain consistent fluxes on the refinement
boundaries, the fluxes on children grids are restricted to parent grids. 
	* Update the solution
		* Advance solution over time step based on the fluxes and source
term, using forward Euler. 
	* Restrict the electron and ion densities to lower levels
	* Fill ghost cells for the electron and ion densities
		* Fill ghost cells on the sides of all boxes.
Near refinement boundaries we interpolate between fine points and
coarse neighbors to fill ghost cells.
The ghost values are less than twice the coarse values. 
Near physical boundaries Neumann boundary conditions are applied.
	* Compute new field on first iteration
		* After the first iteration the electric field should be
recalculated.

\subsection sect_mesh New mesh

The user may decide to compute a new mesh after each time step, but he may
also choose to keep the mesh a fixed number of time steps unchanged. This
can be handled by the variable
<a class="el" href="namespacem__streamer.html#a27ad38736d3b7d879bcfcc8826411d98">st_ref_per_steps</a>.
Its default value is 2.

In case
<a class="el" href="namespacem__streamer.html#a27ad38736d3b7d879bcfcc8826411d98">st_ref_per_steps</a>
<code>= 1</code>, each time step the subroutine
<a class="el" href="streamer__2d_8f90.html#a6b6ddafe30bfd36b772e2bc7888bf5fb">refine_routine</a>
will be called to recompute the mesh.
There is a 'new' mesh if new grids are added or if children grids are deleted. 
See manual.md????, how grids can be deleted. In case new grids are created,
the values of the parent box (like <code>i_photo, i_pos_ion</code> and so on) will be prolongated
to its children by using linear interpolation. 
2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) are applied, which do not need corner ghost cells.
The new mesh needs a recalculation of the electric field by calling 
<a class="el" href="streamer__2d_8f90.html#a3078ef0e5f2dbb937b8b9d54baccaa6c">compute_electric_field</a>.

