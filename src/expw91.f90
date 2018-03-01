! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! c ### expw91 ###

! c     Calculates exchange energy according to a formula of Parr and Wang
! c     Yue Phys. Rev. B 54, 16 533 (1996), Phys. Rev. B 45, 13 244
! c     (1992).

function expw91 (psi,wgt2,rhot,rhotup,rhotdown, &
     grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,iborb,iorb,isiorb,nmut
  real (PREC) :: expw91
  real (PREC) :: const13,const43,const115,ocdown,ocup

  real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,  &
       rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

  real (PREC), external :: expw91sup,fdften

  parameter (const13=1.0_PREC/3.0_PREC,const43=4.0_PREC/3.0_PREC,const115=1.0_PREC/15.0_PREC)
  
  
  do i=1,mxsize
     rhotup(i)  =0.0_PREC
     rhotdown(i)=0.0_PREC
  enddo
  
  !     calculate total densities due to up and down spins
  do iorb=1,norb
     if (inhyd(iorb).eq.1) goto 10
     iborb=i1b(iorb)
     isiorb=i1si(iorb)
     nmut=i1mu(iorb)
     
     call exocc (iorb,ocup,ocdown)
     
     call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
     call scal (isiorb,ocup,wk1,ione)
     
     call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
     call scal (isiorb,ocdown,wk2,ione)
     
     !        store total spin densities 
     call add(isiorb,wk1,rhotup)
     call add(isiorb,wk2,rhotdown)
10   continue
  enddo
  
  !     calculate (nabla rho nabla rho)
  ! FIXME
  call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  call copy(mxsize,wk7,ione,grhotup,ione)
  
  call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  call copy(mxsize,wk7,ione,grhotdown,ione)
  
  do i=1,mxsize
     rhotup(i)   =two*rhotup(i)
     rhotdown(i) =two*rhotdown(i)
     grhotup(i)  =four*grhotup(i)
     grhotdown(i)=four*grhotdown(i)
  enddo
  
  !     total exchange energy is calculated as
  !     Ex(rhoup,rhodown)=(1/2)Ex(2*rhoup)+(1/2)Ex(2*rhodown)
  
  expw91=fdften(alphaf)*(  half*expw91sup(wgt2,rhotup,grhotup,wk0,wk1) &
       + half*expw91sup(wgt2,rhotdown,grhotdown,wk0,wk1))
  
end function expw91
