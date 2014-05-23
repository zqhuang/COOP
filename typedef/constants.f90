module coop_constants
  implicit none
  integer,parameter::coop_integer_length = 4  
  integer,parameter::coop_real_length = 8  ! double precision
  integer,parameter::coop_string_length = 256
  integer,parameter::coop_short_string_length = 32
  integer,parameter::coop_long_string_length = 8192

  COOP_REAL, parameter:: coop_pi = asin(1.d0)*2.d0


end module coop_constants
