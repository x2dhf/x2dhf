module solver
  use params
  integer (KIND=IPREC) :: icase,ifill,muoffs,maxsor1,maxsor2,maxsor3,nthreads,nthreads4coulexch
  integer (KIND=IPREC) :: isstart,isstop,isstep
  
  integer (KIND=IPREC),dimension(4) ::  iadint2,iadint3l,iadint3r,iadint4
  real (PREC), dimension(9,4) ::  cint2,cint3l,cint3r,cint4
  real (PREC) :: omegasf,omegasfOrb,omegasfPot,omega,omega1
end module solver
