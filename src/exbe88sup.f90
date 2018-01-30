! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### exbe88sup ###

!     Calculates exchange energy according to Becke's gradient formula
!     (PRA 88 (1988) 3098-3100)

function exbe88sup (wgt2,rho,grho,wk0,wk1)
  use params
  use discret
  use commons8

  implicit none
  real (PREC) :: exbe88sup
  integer :: i
  real (PREC) :: ash,bbeta,const13,const43,const83,const115,rho43,s,s2
  
  real (PREC), dimension(*) :: wgt2,rho,grho,wk0,wk1
  real (PREC), external :: dot

  parameter (const13=1.0_PREC/3.0_PREC,const43=4.0_PREC/3.0_PREC,const83=8.0_PREC/3.0_PREC, &
       const115=1.0_PREC/15.0_PREC,bbeta=0.0042_PREC)
  
    ash(s)=log(s+sqrt(1.0_PREC+s*s))
  
  !     grho = nabla rho  nabla rho
  !    |nabla rho| = sqrt(grho)

  do i=1,mxsize
     if (abs(rho(i)).lt.precis) then
        wk0(i)=0.0_PREC
     else
        rho43=rho(i)**const43
        s2=grho(i)/(rho43*rho43)
        s=sqrt(s2)
        !           fgga=-bbeta*s2/(one+6.0_PREC*bbeta*s*asinh(s))
        wk0(i)=-bbeta*s2/(one+6.0_PREC*bbeta*s*ash(s))*rho43
     endif
  enddo
  
  !     take care of f4 factor
  call multf4(wk0)
  
  exbe88sup=dot(mxsize,wgt2,ione,wk0,ione)
  
end function exbe88sup
    

