! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### vlpcoeff ###

!     This function uses the Horner scheme to calculate the value of a
!     polynomial defined by coefficients stored in array coeff at a
!     point r0
!
!     vlpcoeff(r0)=coeff(1)+coeff(2)*r0+...+coeff(iord)*r0**(iord-1)

module vlpcoeff_m
  implicit none
contains
  function vlpcoeff (iord,r0,coeff)
    use params

    implicit none
    integer :: i,iord
    real (PREC) :: vlpcoeff
    real (PREC) :: r0
    real (PREC), dimension(*) ::coeff

    vlpcoeff=0.0_PREC
    do i=iord,2,-1
       vlpcoeff=(vlpcoeff+coeff(i))*r0
    enddo
    vlpcoeff=vlpcoeff+coeff(1)

  end function vlpcoeff
end module vlpcoeff_m
