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
! ### lpcoeff ###  

subroutine lpcoeff(iord,istart,k,r,coeff)
  use params
  
  implicit none

  integer :: i,ic1,ic2,iord,istart,j,k
  real (PREC) :: c1,denom,r0
  real (PREC), dimension(12) :: a,b
  real (PREC), dimension(9) :: coeff
  real (PREC), dimension(*) :: r


  !     This routine calculates coefficients of the (sub)Lagrange 
  !     polynomial for a grid point k
  !     \prod_{i=1,i\ne k}^{9} {(r-r_{i}) \over (r_{k}-r_{i})}
  !
  !     r_{k}= r(istart+k), k=1,...,iord
  
  !  calculate denominator product
  
  denom=1.0_PREC
  do i=1,iord
     if (i.ne.k) denom=denom*(r(istart+k)-r(istart+i))
  enddo
  
  do i=1,12
     a(i)=0.0_PREC
     b(i)=0.0_PREC	 
  enddo
  
  !     calculate nominator product: 
  !     a(1)+a(2)*r+a(3)*r^2+...+a(9)*r^8
  !     storing coefficients in a and b
  
  a(1)=1.0_PREC
  
  !  multiply polynomial a by (r+c1), c1=-r(istart+i)
  
  ic2=1	
  do i=1,iord
     if (i.ne.k) then
        ic2=ic2+1
        c1=-r(istart+i)	  
        b(1)=a(1)*c1
        do ic1=2,ic2
           b(ic1)=a(ic1)*c1+a(ic1-1)
        enddo
        b(ic2+1)=a(ic2)	
        do j=1,ic2+1
           a(j)=b(j)
        enddo
     endif
  enddo
  
  do i=1,iord
     coeff(i)=a(i)/denom
  enddo

end subroutine lpcoeff
