! ### checkOrbSym ###
!
!     checks the Ci symmetry of a given orbital

subroutine checkOrbSym (nmut,orb,ihsym)
  use params
  use discret

  implicit none
  integer :: ihsym,imu,inu,nmut
  real (PREC), dimension(nni,nmut) :: orb

  imu=10
  inu=2
  if (orb(inu,imu)*orb(nni-inu+1,imu).gt.0.0_PREC) then
     ihsym= 1
  else
     ihsym=-1
  endif

end subroutine checkOrbSym

