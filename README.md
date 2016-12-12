# config_fortran

A configuration file parser for Fortran. The intended usage is as follows:

1. You create your configuration variables, by providing a default value and
   a description.
2. You read in a text file in which new values are specified for (some of) the
   variables.
3. You use the updated values in your program, so that there is no need to recompile.

Steps 1 and 2 can also be reversed, so that you read in the configuration files
before specifying the variables. Variables can be of type integer, real,
logical/bool, or string, and they can also be an array of such types.

## Example

Suppose you want to use a grid of size `n_grid`, then you could do:

    integer     :: n_grid
    type(CFG_t) :: my_cfg

    call CFG_add(my_cfg, "grid_size", 1024, "Size of the grid")
    call CFG_read_file(my_cfg, "my_input_file.txt")
    call CFG_get(my_cfg, "grid_size", n_grid)

Here, the default grid size will be 1024. If the file `my_input_file.txt` contains a line

    grid_size = 512

the actual grid size used in your program will be 512. This can also be achieved by combining the `add` and the `get` like this:

    integer     :: n_grid = 1024
    type(CFG_t) :: my_cfg

    call CFG_read_file(my_cfg, "my_input_file.txt")
    call CFG_add_get(my_cfg, "grid_size", n_grid, "Size of the grid")

When parsing the input file, the variable `n_grid` will be stored as plain text,
since its type is not yet known. The call `CFG_add_get` converts it to the right
type. See `example_1.f90` and `example_2.f90`for more usage examples.

## Configuration file syntax

There are three types of lines:

1. Blank lines, or lines only containing a comment (`# ...`), which are ignored.
2. Lines indicating the start of a category: `[category_name]`
3. Lines with an `=`-sign. If they are part of a user-defined category, they
   should start with an indent.

An example of a configuration file is shown below

    age = 29
    name = John

    [weather]
        temperature = 25.2
        humidity = 23.5

    happy = .true.

    weather%temperature = 23.9

Note that `temperature` and `humidity` are indented, and that `happy` is not,
which means that `happy` is not part of weather (it is in the default unnamed
category). Any space or tab counts as indentation. Outside an indented
`[weather]` group, you can directly refer to its members by using e.g.
`weather%temperature`, as is done on the last line. To place variables in a
category, you add them like this:

    call CFG_add(my_cfg, "weather%temperature", 25.0_dp, "The temperature")

Variables can also be arrays:

    name_of_variable = value1 [value2 value3 ...] # Optional comment

The extra values `[value2 value3 ...]` are omitted for a scalar variable. You
can create variables of varying array size, by specifying `dynamic_size=.true.`
when creating a config variable:

    call CFG_add(my_cfg, "numbers", [1, 2], "Comment", dynamic_size=.true.)

## Methods

* `CFG_add`: Add a variable to the configuration
* `CFG_get`: Get the value of a variable
* `CFG_add_get`: First `CFG_add`, then `CFG_get`
* `CFG_check`: Check whether all variables have been defined
* `CFG_get_size`: Get the array size of a variable
* `CFG_get_type`: Get the type of a variable
* `CFG_sort`: Sort the configuration (for faster lookup when there are many variables)
* `CFG_write`: Write the configuration to a standard text/config file, which can
  be read in again
* `CFG_write_markdown`: Write the configuration to a file in markdown format
* `CFG_read_file`: Read in a configuration file
* `CFG_update_from_arguments`: Read in the program's arguments as configuration files.

## Requirements

A modern Fortran compiler that supports Fortran 2008. The included `Makefile` was written for `gfortran`.

## TODO

* Write tests

## Alternatives

* [libconfig](http://www.hyperrealm.com/libconfig/) (C/C++)
* [config4*](http://www.config4star.org/) (C/C++)
* [KRACKEN](http://www.urbanjost.altervista.org/LIBRARY/libCLI/arguments/src2015/krackenhelp.html) (Fortran argument parser)
* [FLAP](https://github.com/szaghi/FLAP) (Fortran 2003+ argument parser)
* [FiNeR](https://github.com/szaghi/FiNeR) (Fortran 2003+ config file parser)
