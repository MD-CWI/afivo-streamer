\section sect_intro Introduction 

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

\subsection sect_lookup Lookup tables
TODO

\subsection sect_similarities Similarities of streamer examples

Since the three streamer examples, <code>streamer_2d.f90</code>,
<code>streamer_3d.f90</code>, and <code>streamer_cyl.f90</code> show great similarities, 
a module called <code>m_streamer.f90</code> contains a couple of subroutines,
like the subroutines for handling the configuration files,
getting the electric field, setting the voltages at a given time or to a fixed value,
a new implemenation of the Koren limiter, 
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
      real(dp)              :: bg_dens
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
defined in the refinement subroutine <code>refine_routine</code>, defined by the user, too.

Further the user may decide whether photoionization will be added.
Also here the user has to enclose a subroutine how the photoionization
is defined.
See <code>set_photoionization</code> as an example. 
The photon production rate per cell is proportional to the ionization rate
and will be computed per cell and stored in the <code>i_photo</code> part
for each box.
For more details see <code>set_photoionization_rate</code> as an example. 

\subsection sect_timestepping Time Stepping

The maximum time step is restricted by e.g., the CFL criteria:
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
can be handled by the variable <code>st_ref_per_steps</code>. Its default value is 2.

In case <code>st_ref_per_steps = 1</code>, each time step the subroutine
<code>a2_adjust_refinement</code> will be called to recompute the mesh.
There is a 'new' mesh if new grids are added or if children grids are deleted. 
See manual.md????, how grids can be deleted. In case new grids are created,
the values of the parent box (like i_photo, i_pos_ion and so on) will be prolongated
to its children by using linear interpolation. 
2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) are applied, which
do not need corner ghost cells.
The new mesh needs a recalculation of the electric field by calling 
<code>compute_electric_field</code>

