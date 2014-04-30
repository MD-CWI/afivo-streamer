program test_m_config
  use m_config

  integer, parameter :: dp = kind(0.0d0)
  print *, "Testing m_config.f90 implementation"

  ! Create some parameters
  call CFG_add("first_name", "jannis", "First name of the author of this code")
  call CFG_add("full_name", (/"Jannis   ", "Teunissen"/), "Full name of the author")
  call CFG_add("my_age", 25, "Age of the author of this code")
  call CFG_add("my_fav_reals", (/1.337d0, 13.37d0, 3.1337d0/), "My favorite numbers", varSize=.true.)
  call CFG_add("lots_of_work", .true., "Whether I have a lot of work to do")

  print *, "Original values"
  call CFG_write("stdout")                      ! Write to stdout (only when given the filename "stdout")

  call CFG_read_file("test_m_config_input.txt") ! Update values with file

  print *, "Updated values"
  call CFG_write("stdout")                      ! Write to stdout
  call CFG_write("test_m_config_output.txt")    ! Write to file

  print *, "We can get values in different ways:"
  print *, "lots_of_work:", CFG_get_logic("lots_of_work")
  print *, "my_age:", CFG_get_int("my_age")
  print *, "number of favourite numbers:", CFG_get_size("my_fav_reals")
  print *, "my second favourite number:", CFG_get_real("my_fav_reals", 2)
  print *, "type of full name:", CFG_get_type("full_name")

  call CFG_destruct()
  print *, ""
  print *, "Done, the configuration file has been written to test_m_config_output.txt"
  print *, ""

end program test_m_config
