## fortran_linked_list

A linked list (see <a href="https://en.wikipedia.org/wiki/Linked_list">Wikipedia</a>)
is a linear collection of data elements, called nodes, pointing to the next node by means
of a pointer. It is a data structure consisting of a group of nodes which together
represent a sequence. Under the simplest form, each node is composed of data and
a reference (in other words, a link) to the next node in the sequence.
This structure allows for efficient insertion or removal of elements from any position
in the sequence during iteration. More complex variants add additional links, allowing
efficient insertion or removal from arbitrary element references.

## Example

The package `fortran_linked_list` can easily be used in fortran programs
where linked lists are useful.
It has been developed as a by-product of `a5_streamer`,
a package to simulate streamers. However, the program shown here is completely
separate from `a5_streamer`.

We define a linked list by the following types

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

Besides these types we present four subroutines operating on a linked list: 

	LL_clear    to clear a linked list of type `LL_int_t`
	LL_pop      to add an element to a linked list of type  `LL_int_head_t` at position head
	LL_add      ?????
	LL_get_size to compute the number of elements of a linked list of type `LL_int_head_t`


A test example `test_m_linked_list` has been added to show how a list can be created, added,
read and cleared. If you have a sufficiently recent `gfortran` compiler, you can run the test with

    $ make
    $ ./test_m_linked_list

## Output

The output of this program is:

    Lists do not need to be initialized
    Initial size:           0
    Adding           10  values to list
    Popping values until the list is empty
              10
               9
               8
               7
               6
               5
               4
               3
               2
               1


## More examples

If you are interested on how this package is used by `a5_streamer`
to handle with linked lists, see subroutine: wordt NIET gebruikt

