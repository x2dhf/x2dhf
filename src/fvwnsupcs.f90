! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fvwnsupcs ###

!     Calculates the correlation potential of VWN using the closed-shell
!     formula; it is returned in wk10 (see also fvwncs)

subroutine fvwnsupcs (rhot,rhotup,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
  use params
  use discret
  use commons8

  implicit none
  integer :: i
  real (PREC) :: ck1,cl1,cm1,cn1,const16,const76,const56,constx,constxp,g1,g2,x,xderrhot
  real (PREC), dimension(*) :: wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,rhot,rhotup,rhotdown

  real (PREC), external :: qvwn,qvwnderx

  parameter(ck1=0.03109070_PREC,cl1=-0.104980_PREC,cm1=3.727440_PREC,cn1=12.93520_PREC, &
       const16=1.0_PREC/6.0_PREC,const76=7.0_PREC/6.0_PREC,const56=5.0_PREC/6.0_PREC)

!     total density in rhot

  constx=(three/(four*pii))**const16
  constxp=-(four*pii/three)**const56/(eight*pii)

  ! dx/d(rho)=xderrhot
  do i=1,mxsize
     if (abs(rhot(i)).lt.precis) then
        wk10(i)=0.0_PREC
     else
        x=constx*rhot(i)**(-const16)
        xderrhot=constxp*rhot(i)**(-const76)
        g1=qvwn(x,ck1,cl1,cm1,cn1)
        g2=qvwnderx(x,ck1,cl1,cm1,cn1)
        !            print *,'2',i,constx,constxp,x,xderrhot,g1,g2
        wk10(i)=g1+rhot(i)*xderrhot*g2
     endif
  enddo
end subroutine fvwnsupcs
