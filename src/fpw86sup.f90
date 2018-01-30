! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fpw86sup ###

!     Calculates exchange potential according to a formula of Perdew and
!     Wang Yue (PRB 33(12) (1986) 8800) and returns it in grhot array

subroutine fpw86sup (rhot,grhot,wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,ifive
  real (PREC) :: akf,const13,const43,const115,fdftcoeff,s,t,u,w1,w2

  real (PREC), dimension(*) :: wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,grhot
  real (PREC), external :: ff,ffdp,ffp,fdftpot

  parameter (const13=1.0_PREC/3.0_PREC,const43=4.0_PREC/3.0_PREC,const115=1.0_PREC/15.0_PREC,ifive=5)


  !     Calculate |(nabla rhot nabla rhot)| (wk)
  
  call nfng (rhot,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  call copy(mxsize,wk7,ione,wk,ione)
  
  do i=1,mxsize
     wk(i)=sqrt(wk(i))
  enddo
  
  !     calculate (nabla rho nabla |nabla rho|) (grhot)
  
  call nfng (rhot,wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  call copy(mxsize,wk7,ione,grhot,ione)
  
  !     calculate nabla^2 rho (wk7)
  
  call n2f (rhot,wk0,wk1,wk2,wk3)
  call copy(mxsize,wk3,ione,wk7,ione)
  call scal(mxsize,half,wk7,ione)
  
  fdftcoeff=fdftpot(alphaf)
  
  do i=1,mxsize
     if (abs(rhot(i)).lt.precis) then
        s=0.0_PREC
     else
        !            akf=(three*pii*pii*rhot(i))**const13
        akf=(two*three*pii*pii*rhot(i))**const13
        s=wk(i)/(two*akf*rhot(i))
        t=wk7(i)*(two*akf)**(-2)/rhot(i)
        u=t/(two*akf)/rhot(i)*grhot(i)
        !           write(*,'("s,t,u:",i4,3e14.4)') i,s,t,u
     endif
     
     w1=fdftcoeff*rhot(i)**const13
     
     if (abs(s).lt.precis) then
        wk0(i)=w1*const43*ff(s)
     else
        
        !        if (abs(u).gt.threshPW.or.abs(s).lt.precis) then
        !            wk0(i)=w1*const43*ff(s)
        !        else

        ! FIXME 
        !
        !           Check functional derivative (eq.24) since its second term causes
        !           problems !!!  
        !
        !           12/2007 update 
        !
        !           The problems boils down to s^3 term which blows up
        !           the potential
        
        !            term1=const43*ff(s)
        !            term2=-t/s*ffp(s)
        !            term3=-u*ffdp(s)
        !            term4=const43*s*s*s*ffdp(s)
        
        !            w2=term1+term2+term3+term4
        !            call i2numu(i,in,imu)
        !            if (mod(i,ifive).eq.0) then
        !               write(*,'("term1..term4 ",3i5,5e14.4)') 
        !     &              i,in,imu,term1,term2,term3,term4,w2
        !            endif
        
        w2=(const43*ff(s)-t/s*ffp(s)-(u-const43*s*s*s)*ffdp(s))
        wk0(i)=w1*w2
     endif
  enddo
  
  call copy(mxsize,wk0,ione,grhot,ione)
  
end subroutine fpw86sup

