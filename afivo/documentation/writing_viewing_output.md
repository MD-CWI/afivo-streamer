# Writing and viewing output

Afivo can write output in
the [Silo](https://wci.llnl.gov/simulation/computer-codes/silo) and
the [VTK unstructured](www.vtk.org/VTK/img/file-formats.pdf) format. The Silo format has some advantages:

* Large files can be viewed more efficiently (important in 3D)
* Ghost cells are included in the output, so that contours and surface plots
don't 'break' near refinement boundaries

At the moment, only cell-centered variables are included in the output.

# Different output formats

* `m_af_output::af_write_silo` Write Silo files
* `m_af_output::af_write_vtk` Write VTK files
* `m_af_output::af_write_plane` Write data on a plane
* `m_af_output::af_write_line` Write data along a line
* `m_af_output::af_write_tree` Write the full mesh in binary format (for restarting)
* `m_af_output::af_read_tree` Read the full mesh in binary format (for restarting)

See `m_af_output` for more details.

## Silo output

Silo is the preferred output format. With this format, ghost cells are included so that no clear artifacts should be visible near refinement boundaries.

Some details about this format:
* Only data from the AMR leaves (so the highest visible refinement levels) is included
* In 2D and 3D, multiple quadtree/octree boxes is converted into rectangular regions for efficiency. These are written as so-called *quadmeshes* into a Silo *multimesh*.
* In 1D, output is written in the *curve* format

## VTK output

Write the grid as an unstructured grid in the VTK XML format. Not very efficient, especially in 3D. Does not properly handle the refinement boundaries.

# Visualizing results

The Silo and VTK files can be visualized
with [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads).
Consult
the [Visit manual](https://wci.llnl.gov/simulation/computer-codes/visit/manuals)
or one of
the [tutorials](http://www.visitusers.org/index.php?title=Short_Tutorial) for
more details.
