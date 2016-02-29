  integer,parameter::coop_short_int_length = 2
  integer,parameter::coop_integer_length = kind(1)  
  integer,parameter::coop_long_integer_length = 8
  integer,parameter::coop_real_length = kind(1.d0)  ! double precision
  integer,parameter::coop_single_real_length = kind(1.)  ! single precision
  integer,parameter::coop_single_complex_length = kind( (1., 1.) )  ! double precision
  integer,parameter::coop_complex_length = kind( (1.d0, 1.d0) )  ! double precision
  integer,parameter::coop_string_length = 1024
  integer,parameter::coop_short_string_length = 32
  integer,parameter::coop_long_string_length = 16384
  real(coop_real_length), parameter::coop_LogZero = 1.d30
    
