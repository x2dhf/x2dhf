!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Copyright (C) 2008 Tomasz Dziubak <tomek@fizyka.umk.pl> (C version)   !
!   Copyright (C) 2017 Susi Lehtola (converted to Fortran 90)             !
!                                                                         !
!   This program is free software; you can redistribute it and/or modify  !
!   it under the terms of the GNU General Public License version 2 as     !
!   published by the Free Software Foundation.                            !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module kh
  use params, only: PREC
  implicit none

  real (PREC), parameter :: pii = 4.0_PREC*atan(1.0_PREC)

contains

! This function returns a value of the model potential -V0/Sqrt(a^2+r^2)
! transformed into the cylindrical coordinates, i.e. the value
!         -V0/Sqrt(a^2+z^2+s^2)
! plus a m-dependent correction resulting from the azimuthal variable
! being factored out.
!
! PARAMETERS
! z,s - cylindrical coordinates
!   m - magnetic quantum number
!  V0 - potential depth
!   a - potential core width
!	 if a=0, V0=1 -> Coulomb potential
!	 if a>0, V0>0 -> smoothed Coulomb potential
  function poth3(z, s, m, V0, a, precis)
    implicit none
    real (PREC), intent(in) :: z, s, V0, a, precis
    integer, intent(in) :: m
    real (PREC) :: popr_cylind, poth3

    if (abs(s) .lt. precis .and. abs(z) .lt. precis) then
       popr_cylind = (m**2)/(2.0*precis**2)
       poth3 = -V0 / sqrt(a**2 + precis**2) + popr_cylind
    else
       if (s .lt. precis) then
          popr_cylind=(m**2)/(2.0_PREC * precis**2)
       else
          popr_cylind=(m**2)/(2.0_PREC * s**2)
       end if
       poth3=-V0/sqrt(a**2 + z**2 + s**2) + popr_cylind;
    end if

    return
  end function poth3

! This function returns a value of the Kramers-Henneberger potential
! at a point (z,s) in cylindrical coordinates
! for the smoothed Coulomb potential plus a m-dependent correction
! correction resulting from the azimuthal variable being factored
! out.
!
! PARAMETERS
!  z,s - cylindrical coordinates
!    m - magnetic quantum number
!  eps - laser field intensity
!    w - laser cycle frequency
!   V0 - original potential depth
!    a - original potential core width (a>0)
!    N - number of intervals in the Simpson quadrature
  function potkh(z, s, m, eps, w, V0, a, N, precis)
    implicit none
    real (PREC), intent(in) :: z, s, eps, w, V0, a, precis
    integer, intent(in) :: m, N
    real (PREC) :: popr_cylind, potkh

    if (abs(s) .lt. precis) then
       popr_cylind = (m**2)/(2.0_PREC * precis**2)
       potkh = popr_cylind + AM_HK_ZSimpson(z, precis, 0.0_PREC, eps, w, V0, a, N)
    else
       popr_cylind = (m**2)/(2.0_PREC * s**2);
       potkh = popr_cylind + AM_HK_ZSimpson(z, s, 0.0_PREC, eps, w, V0, a, N)
    end if

    return
  end function potkh

  ! AUXILIARY FUNCTIONS
  !-------------------------------------------------
  ! HK in integral form
  function AM_HK(t, x, y, z, e0, w, V0, a)
    implicit none
    real (PREC), intent(in) :: t, x, y, z, e0, w, V0, a
    real (PREC) :: alpha0, bx, AM_HK

    alpha0 = e0 / w**2
    bx = x + alpha0 + alpha0 * (cos(w*t)-1.0_PREC)
    AM_HK = -V0 / ( (2.0_PREC * pii / w) * sqrt(a**2 + bx**2 + y**2 + z**2) )

    return
  end function AM_HK

  ! Numerical integration by means of the composite Simpson quadrature
  function AM_HK_ZSimpson(x, y, z, e0, w, V0, a, Nt)
    real (PREC), intent(in) :: x, y, z, e0, w, V0, a
    integer, intent(in) :: Nt
    integer :: i

    ! Interval length
    real (PREC) :: dT
    ! Result
    real (PREC) :: AM_HK_ZSimpson

    ! Calculate interval
    dT = 2*pii/(w*Nt)

    ! Values at the end points
    AM_HK_ZSimpson = AM_HK(0.0_PREC, x, y, z, e0, w, V0, a) &
         + AM_HK(Nt*dT, x, y, z, e0, w, V0, a)
    ! Values at the intermediate points
    do i=1,Nt-1
       AM_HK_ZSimpson = AM_HK_ZSimpson + 2.0_PREC * AM_HK(i*dT,x,y,z,e0,w,V0,a)
    end do
    ! Values at the mid points
    do i=0,Nt-1
       AM_HK_ZSimpson = AM_HK_ZSimpson + 4.0_PREC * AM_HK((i+0.5_PREC)*dT,x,y,z,e0,w,V0,a)
    end do
    ! Normalization
    AM_HK_ZSimpson = AM_HK_ZSimpson * dt / 6.0_PREC
    
    return
  end function AM_HK_ZSimpson
end module kh
