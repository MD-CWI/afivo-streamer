# Parallelization with OpenMP

Most operations in Afivo loop over a number of boxes, for example the leaves at
a certain refinement level. All such loops have been parallelized by adding
OpenMP statements around them, for example as shown below:

```{.f90}
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call my_procedure(tree%boxes(id))
       end do
       !$omp end do
    end do
    !$omp end parallel
```

The parallel speedup that one can get depends on the cost of the algorithm that
one is using. The communication cost (updating ghost cells) is always about the
same, so that an expensive algorithm will show a better speedup. Furthermore, on
a shared memory system, it is not unlikely for an algorithm to be memory-bound
instead of CPU-bound. An example of the scaling can be found
