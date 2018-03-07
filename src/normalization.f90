! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *   Copyright (C) 2018- Susi Lehtola                                      *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************

module normalization
  implicit none
contains

  ! Normalization factor for Laguerre polynomials
  function laguerre_normalization(z,n,l) result(fn)
    use params
    use commons8
    use factor_m
    implicit none

    real(PREC), intent(in) :: z
    integer, intent(in) :: n, l
    real(PREC) :: fn

    fn=(2.0_PREC*z/dble(n))**(3.0_PREC/2.0_PREC+dble(l))*sqrt(factor(n+l)/(2.0_PREC*dble(n)*factor(n-l-1)))/factor(2*l+1)
  end function laguerre_normalization

  ! Normalization factor for spherical harmonics
  function sphharm_normalization(l,m) result(shn1)
    use params
    use factor_m
    implicit none

    integer, intent(in) :: l, m
    real(PREC) :: shn1

    if(m.eq.0) then
       shn1=1.0_PREC/sqrt(4.0_PREC*pii)*sqrt((2*l+1)*factor(l-m)/factor(l+m))
    else
       shn1=(-1.0_PREC)**dble(m)/sqrt(4.0_PREC*pii)*sqrt((2*l+1)*factor(l-m)/factor(l+m))
    end if
  end function sphharm_normalization
end module normalization


