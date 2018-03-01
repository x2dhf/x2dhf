! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### flp ###  

!     It is assumed that function f(r) is tabulated at grid points r_i,
!     i=1,...,n which are stored in array r and the corresponding values
!     f_i are in array f.

!     The function fr0 uses the Lagrange polynomial of a given order to
!     calculate f(r0). The interpolation employs grid points adjacent to
!     r0.

!   iord - number of grid point used for interpolation 
!     n  - dimension of arrys r and f
!     r  - array of abscissas
!     f  - array of corresponding values of a given function f
!     r0 - an abscissa for which f(r0) is calculated via 8th-order 
!          Lagrange polynomial
 
function flp(iord,n,r,f,r0)
  use params

  implicit none
  integer :: i,iord,istart,k,n,nearest
  real (PREC) :: flp
  real (PREC) :: r0
  real (PREC), dimension(9) :: coeff
  real (PREC), dimension(9,9) :: coeff2
  real (PREC), dimension(n) :: r,f

  real (PREC), external :: vlpcoeff

  !     iord value for 2th order (3-point) Lagrange polynomial
  !     parameter (iord=3)
  
  !     iord value for 8th order (9-point) Lagrange polynomial
  !     parameter (iord=9)
  
  
  !     find the smallest element of arry r that is greater than
  !     r0. Whenever possible [iord/2] points to its left and right are
  !     used to construct the interpolation polynomial (the start and the
  !     end of r array are taken care of)
  
  nearest=1
  do i=1,n
     if (r(i).gt.r0) then
        nearest=i
        goto 10
     endif
  enddo
10 continue
  
  !     calculate coefficients of iord Lagrange polynomials (begining with
  !     istart+1) and store them in coeff2
  
  istart=nearest-(iord/2+1)
  if (istart.lt.1) istart=0
  if (istart+iord.gt.n) istart=n-iord
  do k=1,iord
     call lpcoeff(iord,istart,k,r,coeff)
     do i=1,iord
        coeff2(i,k)=coeff(i)
     enddo
  enddo

  !     evaluate the value of the Lagrange interpolation polynomial at r0
  
  flp=0.0_PREC
  do k=1,iord
     flp=flp+f(istart+k)*vlpcoeff(iord,r0,coeff2(1,k))
  enddo
  return
end function flp


