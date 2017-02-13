# Writing and viewing output

Afivo can write output in
the [Silo](https://wci.llnl.gov/simulation/computer-codes/silo) and
the [VTK unstructured](www.vtk.org/VTK/img/file-formats.pdf) format. The Silo format has some advantages:

* Large files can be viewed more efficiently (important in 3D)
* Ghost cells are included in the output, so that contours and surface plots
don't 'break' near refinement boundaries

At the moment, only cell-centered variables are included in the output.

# Writing Silo/VTK files

The routines for writing Silo and VTK files are m_a2_output::a2_write_silo and
m_a2_output::a2_write_vtk, with equivalents in 3D. These routines can be called
as follows:

    ! To write "output/test.silo"
    call a2_write_silo(tree, "test", dir="output")

    ! To write "output/test.vtu"
    call a2_write_vtk(tree, "test", dir="output")

They write the full mesh structure, and by default include all cell-centered
variables. Because files can get pretty big in 3D, there is an optional argument
`ixs_cc`. By specifying e.g. `ixs_cc = [1,2]` only the first two variables will
be included in the output. There are also optional arguments `n_cycle` and
`time`, which correspond to the simulation cycle (number) and the simulation
time.

# Writing a plane

Sometimes it is convenient to write the solution in some region with a uniform
resolution. For this purpose there is the routine m_a2_output::a2_write_plane,
which writes a VTK file.

# Visualizing results

The Silo and VTK files can be visualized
with [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads).
Consult
the [Visit manual](https://wci.llnl.gov/simulation/computer-codes/visit/manuals)
or one of
the [tutorials](http://www.visitusers.org/index.php?title=Short_Tutorial) for
more details.



