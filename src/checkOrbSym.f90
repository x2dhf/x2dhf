! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2023  Jacek Kobus 

module checkOrbSym_m
  implicit none
contains
  ! ### checkOrbSym ###
  !
  !     checks the Ci symmetry of a given orbital
  !
  subroutine checkOrbSym (nmut,orb,ihsym)
    use params
    use discrete

    implicit none
    integer (KIND=IPREC) :: ihsym,imu,inu,nmut
    real (PREC), dimension(nni,nmut) :: orb
    ! imu=13 corresponds to the smallest grid that can be used
    imu=13
    inu=2
    if (orb(inu,imu)*orb(nni-inu+1,imu).gt.0.0_PREC) then
       ihsym= 1
    else
       ihsym=-1
    endif

  end subroutine checkOrbSym
end module checkOrbSym_m
