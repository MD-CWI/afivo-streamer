# Profiling

With [gperftools](https://github.com/gperftools/gperftools) installed, profiling can be done by first recompiling:

    make clean
    make PROF=gperftools

And then executing the following commands, where `my_executable` could be one of the included examples:

    CPUPROFILE=ls.prof ./my_executable

    # Get performance data in callgrind format
    pprof --callgrind ./my_executable ls.prof > callgrind.out

    # Inspect results with the kcachegrind GUI
    kcachegrind
