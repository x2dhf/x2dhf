! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fbe88sup ###

!     Calculates the Becke exchange potential for a given density and
!     returns it in grhot array (see Johnson, Gill, Pople, JCP 98 (1993)
!     p.5623) and L.Fan and T.Ziegler, JCP 94 (1991) 6057)

subroutine fbe88sup (rhot,grhot,wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,iformA6
  real (PREC) :: ash,bbeta,const13,const23,const32,const34,const43,const53,const115, &
       g1,g3,f,fm1,s,s2,t1,t2,t3

  real (PREC), dimension(*) :: wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,grhot
  real (PREC), external :: fdftpot
  parameter(iformA6=2)
  parameter (const13=1.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC,const32=3.0_PREC/2.0_PREC, &
       const34=3.0_PREC/4.0_PREC,const43=4.0_PREC/3.0_PREC,const53=5.0_PREC/3.0_PREC,       &
       const115=1.0_PREC/15.0_PREC,bbeta=0.0042)


  ! arcsinh
  ash(s)=log(s+sqrt(one+s*s))
  
  !     Calculate |(nabla rhot nabla rhot)| (grhot)
  call nfng (rhot,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhot)
  
  !     calculate gamma (see JGP, p.5623)
  do i=1,mxsize
     if (abs(rhot(i)).lt.precis) then
        wk(i)=0.0_PREC
     else
        wk(i)=sqrt(grhot(i))/rhot(i)**const43
     endif
  enddo
  
  !     nabla^2 rho
  call n2f (rhot,wk0,wk1,wk2,grhot)
  
  !     nabla gamma nabla rho (wk7)
  call nfng (wk,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  
  !     calculate g and g prime
  do i=1,mxsize
     g1=fdftpot(alphaf)*rhot(i)**const13
     
     s=wk(i)
     s2=s*s
     
     if (abs(rhot(i)).lt.precis) then
        wk1(i)=g1
     else
        fm1=(one+6.0_PREC*bbeta*s*ash(s))
        if (fm1.lt.precis) then
           f=0.0_PREC
        else
           f=one/fm1
        endif
        
        g3=s/sqrt(one+s2)
        
        t1=const43*s2*rhot(i)**const53
        !           nabla^2 rho
        t2=one+f*(one-6.0_PREC*bbeta*s*g3)
        
        !           nabla rho nabla gamma
        t3=6.0_PREC*bbeta*f*((one+two*f)*ash(s)+g3*(one/(one+s2)+two*f*(two-6.0_PREC*bbeta*s*g3)))
        
        wk1(i)=g1-bbeta*f/rhot(i)**const43*(t1-grhot(i)*t2+wk7(i)*t3)
        
     endif
  enddo
  
  call copy(mxsize,wk1,ione,grhot,ione)
  
end subroutine fbe88sup

