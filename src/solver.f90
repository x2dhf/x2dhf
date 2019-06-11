module solver
  use params

  integer :: icase,ifill,isym,muoffs
  integer, dimension(4) ::  iadint2,iadint3l,iadint3r,iadint4
!  real (PREC), dimension(4) ::  dmu2t,dmu1t
  real (PREC), dimension(9,4) ::  cint2,cint3l,cint3r,cint4
  real (PREC) :: omegasf,omegasfOrb,omegasfPot,omegasf4lda,omega,omega1
end module solver
