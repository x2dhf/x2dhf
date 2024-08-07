module sharedMemory
  use params

  integer (KIND=IPREC), dimension(:), pointer :: sorptr=>null()
  real (PREC), dimension(:), pointer :: coulombptr,exchangeptr,orbptr,coulptr,exchptr,supplptr,&
       scratchptr,scratch4lxcptr
  common /c_sorptr/ sorptr
  common /c_orbptr/ orbptr
  common /c_exchptr/ exchptr
  common /c_supplptr/ supplptr  
  common /c_scratchptr/ scratchptr
  common /c_scratch4lxcptr/ scratch4lxcptr  
  common /c_coulptr/ coulptr
  common /c_coulombptr/ coulombptr
  common /c_exchangeptr/ exchangeptr
  common /c_i4b/ i4barr(20)
  !common /c_iadex/ iadex1c,iadex2c,iadex3c
  common /c_iadex/ iadex1c,iadex2c,iadex3c,iadextc,iadnorc
  integer (KIND=IPREC), dimension(20) :: i5bc,i4ec,i5ec,i4sic,i5sic,i4ngc,i5ngc
  integer (KIND=IPREC) :: iadex1c,iadex2c,iadex3c,iadextc,iadnorc
end module sharedMemory
