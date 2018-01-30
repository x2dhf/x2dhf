! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### expw86sup ###

!     Calculates exchange energy according to a formula of Parr and Wang
!     Yue (PRB 33 (1986) 8800)

function expw86sup (wgt2,rho,grho,wk0,wk1)
  use params
  use discret
  use commons8

  implicit none
  integer :: i
  real (PREC) :: expw86sup
  real (PREC) :: a,b,c,const,const13,const23,const43,const83,const115,fgga,s,s2
  real (PREC), dimension(*) :: wgt2,rho,grho,wk0,wk1
  real (PREC), external :: dot

  parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC,const13=1.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC, &
       const43=4.0_PREC/3.0_PREC,const83=8.0_PREC/3.0_PREC,const115=1.0_PREC/15.0_PREC)

!      grho = nabla rho  nabla rho
!      |nabla rho| = sqrt(grho)

  const=(24.0_PREC*pii*pii)**(-const23)
  do i=1,mxsize
     if (abs(rho(i)).lt.precis) then
        wk0(i)=0.0_PREC
     else
        s2=const*grho(i)/rho(i)**const83
        s=sqrt(s2)
        fgga=(one+a*s2+b*s2*s2+c*s2*s2*s2)**const115

        ! Parr & Yue 1986  
        wk0(i)=rho(i)**const43*fgga

        ! Langreth-Mehl    
        !           wk0(i)=rho(i)**const43*(one+1.521*0.0864*s2)
        
        ! GEA               
        !           wk0(i)=rho(i)**const43*(one+0.0864*s2)
        
        ! LDA
        !           wk0(i)=rho(i)**const43
     endif
  enddo
  
  !     take care of F4 factor
  call multf4(wk0)
  
  expw86sup=dot(mxsize,wgt2,ione,wk0,ione)
  
end function expw86sup


