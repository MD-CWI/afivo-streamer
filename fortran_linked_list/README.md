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
It has been developed as a by-product of <b>`a5_streamer`</b>,
a package to simulate streamers. However, the program shown here is completely
separate from <b>`a5_streamer`</b>.

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

Besides these types we present for subroutines operating on a linked list: 

	LL_clear    to clear a linked list of type `LL_int_t`
	LL_pop      to add an element to a linked list of type  `LL_int_head_t` at position head
	LL_add      ?????
	LL_get_size to compute the number of elements of a linked list of type
`LL_int_head_t`


An example how to add data is shown in test program
<a class="el" href="test__lookup__table_8f90_source.html">test_lookup_table</a>.
The result of the first part of this program
is:


\subsection sect_example More examples

If you are interested on how this package is used by  `a5_streamer`
to create a lookup table with transport data, see subroutine 
<a class="el" href="namespacem__streamer.html#ac543e682ffced5108a9e5f33b7c6c1ba">st_load_transport_data</a>

