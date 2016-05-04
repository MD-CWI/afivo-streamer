### Design considerations

#### Why do not you fill corner ghost cells by default?

Filling these is relatively easy in 2D, but quite a pain in 3D. There you have to
consider 8 corner points and 12 edges between these corners. These edges can be
shared by multiple neighbors, and filling them in a consistent way is quite
difficult.

#### Why use Fortran (2003+)?

Because it is one of the more convenient languages for scientific computing.

#### Why do not you use MPI?

There are a couple of reasons for this:

* There are already many frameworks out there aimed at "big" simulations.
* I think most people dont need "big" simulations, running on more than say 10
  cores.
* Most of the complexity of the current frameworks is in the communication, this
  is much simpler for AFiVO. There is much less code, and it is probably easier
  to make changes in a project if you can read all the data from each core, so
  that you do not have to think about MPI. (Although getting good OpenMP
  performance is also quite tricky, I admit).
* When your code is more efficient, you can use a smaller system to do the same
  type of simulations. This is what I hope to achieve.
* Most parallel codes do not scale so well, especially if there is a lot of grid
  refinement. The work is then harder to distribute, and more communication is
  required.
* If your simulation fits in memory, then you can also consider running 5
  different test cases on a big system instead of one bigger one. The bigger one
  will almost always be less efficient.

#### Why use one fixed ghost cell?

I have considered a couple of options, which are listed below with some remarks:

* Variable number of ghost cells (depending on the variable)

	* Suitable for all ghost cell requirements, flexible.
	* Having more than 1 ghost cells is not very memory efficient. For example,
	  for a 2D 8x8 block with 2 ghost cells per side, you would have to store
	  12x12 cells (144). So the flexibility of having more than one is not really
	  all that useful.
	* It is harder to write code for a variable number of ghost cells, for
      example, when copying data, should we copy the ghost cells? And writing
      code for filling more than one ghost cell is also hard.
	* Storing variables is annoying, because they cannot be stored in the same
      array (since they have different shapes). This makes indexing harder.

* One ghost cell

	* Restricted, not flexible.
	* Simple to implement because all variables are the same.
	* One ghost cells does not cost too much memory.
	* One ghost cells is convenient for 2nd order schemes, such as the multigrid
      examples, which are quite common.
	* When you need more than one ghost cell, you can simply access the data on
      a neighbor directly. See for example the drift-diffusion test.
	* Perhaps I will add variables without ghost cells in the future, to not
      waste memory for them.

* No ghost cells

	* Perhaps the most general / elegant idea: do not waste memory on ghost cells
      but just look the values up at the neighbors.
	* Hard to write efficient code: typically you would work on an enlarged copy
      of the current box that includes neighbor data. Copying data takes time,
      and it is hard to write elegant routines for this. For example, to get a
      corner ghost cell, you typically want to use the "side" ghost cell of a
      neighbor, but if these are not stored, they have to be recomputed each
      time.
	* If you do not work on an enlarged copy of the box, indexing is really
      annoying.

NEWNEWNEWNEW

\subsection sect_motivation-history Motivation: a brief history

Given the fact that there are already several frameworks available, why develop
another one? The main reason was that a *simple* or *basic* framework
seemed to be missing. Our motivation came from work on time-dependent
simulations of streamer discharges. These discharge have a multiscale nature,
and require a fine mesh in the region where they grow. Furthermore, at every
time step Poisson''s equation has to be solved. A streamer model that uses a
uniform Cartesian grid is therefore computationally expensive, especially for 3D
simulations.

In Xcite{Pancheshnyi_2008}, Paramesh was used for streamer simulations. The main
bottleneck in this implementation was however the Poisson solver. Other streamer
models (see e.g., Xcite{Montijn_2006, Li_hybrid_ii_2012, Luque_2012}) had the
same problem, because the non-local nature of Poisson''s equation makes an
efficient parallel solution difficult, especially on an adaptively refined grid.
An attractive solution method to get around this is geometric multigrid,
discussed in section \ref sect_afivo-multigrid.

We first considered implementing multigrid in Paramesh Xcite{Macneice_2000},
which already includes an *alpha* version of a multigrid solver with the
following comment Xcite{www_paramesh_mg}:
> ##
>  This is an ALPHA version of this feature.
>  You should be aware that it may be **buggy**.
>  Also, construction of multigrid algorithms and AMR is much less
>  straightforward than incorporating AMR into finite-volume hydro codes.

Because Paramesh does not seem to be actively maintained, we decided not to move
forward with it after experiencing several problems (writing and visualizing
output, having a fixed number of blocks).

NEWNEWNEWNEW
Next, we considered Boxlib Xcite{www_boxlib}, an actively maintained framework
which is also used in Chombo Xcite{chombo}. Boxlib contains a significant amount
of multigrid code, including several examples that demonstrate how a solver can
be set up and used. After spending some time getting familiar with the
framework, we tried to modify the multigrid solver to our needs. This involves
operations like: get the coarse grid values next to refinement boundaries to
perform a special type of ghost cell filling (see section
\ref sect_mg-ghost-cells). An operation like this can definitely be implemented
in Boxlib, but we realized that this task can be made easier in several ways.

TODO: explain clearly and positively!

they are not trivial to implement. There are a number of reasons for
this. First of all, Boxlib is a quite flexible framework used in diverse
applications, which increases the number of supported features and thereby the
size of the code. Furthermore, block-structured grids are used, whose
connectivitty is not as simple as that of orthtrees. This is the , MPI is used
for parallelization, and .

In our experience, a large number of scientific simulations fit into the memory
of a desktop machine or cluster node, which nowadays often have 16 or 32
gigabytes of RAM. A practical reason for this is that for larger problems, the
visualization of the results becomes quite challenging. In general, writing
parallel code with good scaling takes considerable effort, for which the
manpower and resources are often lacking. For the application we had in mind,
efficient large scale parallelism would anyway be hard, due to the non-locality
of Poisson''s equation. This inspired us to develop a framework that uses
shared-memory parallelism, which makes many operations much simpler, because all
data can directly be accessed. The goal was to create a relatively simple
framework that could easily be modified, to provide an option in between the
**advanced** distributed-memory codes of table Xref{tab:amr-frameworks and simple
uniform grid computations}.


