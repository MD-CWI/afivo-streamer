program test_m_config
  use m_config

  integer, parameter :: dp = kind(0.0d0)
  type(CFG_t) :: my_cfg
  integer :: n_reals, v_type
  real(dp), allocatable :: my_reals(:)
  logical :: my_logic
  integer :: my_int
  print *, "Testing m_config.f90 implementation"

  ! Create some parameters
  call CFG_add(my_cfg, "first_name", "jannis", "First name of the author of this code")
  call CFG_add(my_cfg, "full_name", (/"Jannis   ", "Teunissen"/), "Full name of the author")
  call CFG_add(my_cfg, "my_age", 25, "Age of the author of this code")
  call CFG_add(my_cfg, "my_fav_reals", (/1.337d0, 13.37d0, 3.1337d0/), "My favorite numbers", dyn_size=.true.)
  call CFG_add(my_cfg, "lots_of_work", .true., "Whether I have a lot of work to do")

  print *, "Original values"
  call CFG_write(my_cfg, "stdout")                      ! Write to stdout (only when given the filename "stdout")

  call CFG_read_file(my_cfg, "test_m_config_input.txt") ! Update values with file

  print *, "Updated values"
  call CFG_write(my_cfg, "stdout")                      ! Write to stdout
  call CFG_write(my_cfg, "test_m_config_output.txt")    ! Write to file

  print *, "We can get values in different ways:"
  call CFG_get(my_cfg, "lots_of_work", my_logic)
  print *, "lots_of_work:", my_logic

  my_logic = CFG_fget_logic(my_cfg, "lots_of_work")
  print *, "lots_of_work:", my_logic

  call CFG_get(my_cfg, "my_age", my_int)
  print *, "my_age:", my_int

  my_int = CFG_fget_int(my_cfg, "my_age")
  print *, "my_age:", my_int

  call CFG_get_size(my_cfg, "my_fav_reals", n_reals)
  print *, "number of favourite numbers:", n_reals

  n_reals = CFG_fget_size(my_cfg, "my_fav_reals")
  print *, "number of favourite numbers:", n_reals

  allocate(my_reals(n_reals))
  call CFG_get(my_cfg, "my_fav_reals", my_reals)
  print *, "my favourite number:", my_reals(1)
  print *, "my second favourite number:", my_reals(2)
  call CFG_get_type(my_cfg, "full_name", v_type)
  print *, "type of full name:", v_type

  print *, "Done, the configuration file has been written to test_m_config_output.txt"
  print *, ""
end program test_m_config
