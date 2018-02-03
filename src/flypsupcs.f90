! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### flypsupcs ###

!     Calculates the correlation potential of Lee, Yang and Parr (PRB 37
!     (1988) 785, eq.23, i.e. using the closed-shell formula); it is
!     returned in wk10

subroutine flypsupcs (rhot,rhotup,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
  use params
  use discret
  use commons8

  implicit none
  integer :: i
  real (PREC) :: a,b,c,cf,d,const13,const23,const53,const83, &
       f1,g0,g1,g2,g3,rh


  real (PREC), dimension(*) :: wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10, &
       rhot,rhotup,rhotdown

  parameter (const13=1.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC,const53=5.0_PREC/3.0_PREC, &
       const83=8.0_PREC/3.0_PREC)

  !     coefficients of the Colle-Salveti formula
  parameter(a=0.049180_PREC,b=0.1320_PREC,c=0.25330_PREC,d=0.3490_PREC)


  !     total density in rhot

  !     nabla^2 rho (rhotup)
  call n2f (rhot,wk3,wk4,wk5,wk6)
  call copy (mxsize,wk6,ione,rhotup,ione)

  !     nabla rho nabla rho (rhotdown)
  call nfng (rhot,rhot,rhotdown,wk3,wk4,wk5,wk6,wk7,wk8,wk9)
  call copy (mxsize,wk9,ione,rhotdown,ione)


  cf=0.30_PREC*(three*pii*pii)**const23

  do i=1,mxsize
     if (abs(rhot(i)).lt.precis) then
        wk10(i)=0.0_PREC
     else
        rh=rhot(i)**(-const13)
        f1=one/(one+d*rh)
        g0=rhot(i)**(-const53)*exp(-c*rh)
        g1=-const53+c/three*rhot(i)**(-const13)
        g2=g1+d/three*f1*rhot(i)**(-const13)
        g3=d*d*f1*f1*rhot(i)**(-const13) -d*f1-c

        wk10(i)=                                                     &

             ! See documentation for the derivation of these terms
             !        term (1)
             -a*d/three*rh*f1*f1-a*f1                                &

             !        term (2)
             -a*b*cf*f1*g0*rhot(i)**(const53)*(g2+const83)           &

             !        term (3)
             -a*b*19.0_PREC/18.0_PREC*f1*g0*rhotup(i)                        &

             !        term (4)
             -a*b*1.0_PREC/72.0_PREC*f1*g0*g2*(42.0_PREC*rhotup(i)+59.0_PREC*rhotdown(i)/rhot(i)) &

        !        term (5)
             -a*b*7.0_PREC/24.0_PREC*f1*g0/rhot(i)*rhotdown(i)*( g2*(d/three*f1*rh+g1-1.0_PREC)+g3*rh/9.0_PREC )

     endif
  enddo

  return
end subroutine flypsupcs

