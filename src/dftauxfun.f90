! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### qvwn ###
!
!     Calculates q function of VWN functional (see
!     http://wild.life.nctu.edu.tw/~jsyu/molpro2002.1/doc/manual/node184.html)

function qvwn(x,a,p,c,d)
  use params
  use commons8
  implicit none
real (PREC) :: qvwn
  real (PREC) :: a,arctan,c,d,p,qcap,t1,t2,t3,xcap,x

  xcap=x*x+c*x+d
  qcap=sqrt(four*d-c*c)      
  arctan=atan(qcap/(two*x+c))
  
  t1=log(x*x/xcap)
  t2=two*c/qcap*arctan
  t3=-c*p/xcap*(log((x-p)*(x-p)/xcap) +two*(c+two*p)/qcap*arctan)
  
  qvwn=a*(t1+t2+t3)
  return
end function qvwn

! ### qvwnderx ###

!     Calculates dq/dx of the VWN functional 

function qvwnderx(x,a,p,c,d)
  use params
  use commons8

  implicit none
  real (PREC) qvwnderx
  real (PREC) :: a,c,d,p,qcap,qcap2,t1,t2,t3,t4,x,xcap,xcapder,xcapder2,xcapp

! X=xcap, dX/dx=xcapder 
  xcap=x*x+c*x+d
  xcapder=two*x+c
  xcapder2=xcapder*xcapder

  xcapp=p*p+c*p+d
  
  qcap=sqrt(four*d-c*c)      
  qcap2=(four*d-c*c)      
  
  t1=two/x-xcapder/xcap
  t2=-four*c/(qcap2+xcapder2)
  t3=four*c*p*(c+two*p)/xcapp/(qcap2+xcapder2)
  t4=-c*p/xcapp*(two/(x-p)-xcapder/xcap)
  
  qvwnderx=a*(t1+t2+t3+t4)

end function qvwnderx

! ### ff ###
!     Calculates F(s) (see fpw86sup)

function ff (s)
  use params

  implicit none
  real (PREC) :: ff
  real (PREC) ::  const115,s
  real (PREC), external :: ffbar

  parameter (const115=1.0_PREC/15.0_PREC)

  ff=ffbar(s)**const115
  
end function ff

! ### ffbar ###
!     Calculates Fbar(s) (see fpw86sup)

function ffbar (s)
  use params
  implicit none 
  real (PREC) :: ffbar
  real (PREC) ::  a,b,c,s,s2
  parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC)

  s2=s*s
  ffbar=(1.0_PREC+a*s2+b*s2*s2+c*s2*s2*s2)
  
end function ffbar

! ### ffbarp ###
!     Calculates Fbar(s) prime (see fpw86sup)

function ffbarp (s)
  use params
  implicit none
  real (PREC) :: ffbarp
  real (PREC) :: a,b,c,s,s2
  parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC)

  s2=s*s
  ffbarp=(2.0_PREC*a*s+4.0_PREC*b*s*s2+6.0_PREC*c*s*s2*s2)
  
end function ffbarp


! ### ffp ###
!     Calculates dF(s)/ds (see fpw86sup)

function ffp (s)
  use params
  implicit none
  real (PREC) :: ffp
  real (PREC) :: a,b,c,const13,const43,const115,s,s2
  real (PREC), external :: ffbar,ffbarp
  parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC,const13=1.0_PREC/3.0_PREC, &
       const43=4.0_PREC/3.0_PREC,const115=1.0_PREC/15.0_PREC)
  
  s2=s*s
  ffp=const115*ffbar(s)**(const115-1.0_PREC)*ffbarp(s)
end function ffp


! ### ffdp ###
!     Calculates d(1/s dF(s)/ds)/ds (see fpw86sup)

function ffdp (s)
  use params

  implicit none
  real (PREC) :: ffdp
  real (PREC) :: a,b,c,const13,const43,const115,s,s2
  real (PREC), external :: ffbar,ffbarp

  parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC,const13=1.0_PREC/3.0_PREC, &
       const43=4.0_PREC/3.0_PREC,const115=1.0_PREC/15.0_PREC)
  
  s2=s*s
  ffdp= const115*(8.0_PREC*b*s+24.0_PREC*c*s*s2)*ffbar(s)**(const115-1.0_PREC) &
       +const115*(const115-1.0_PREC)*ffbarp(s)**2.0_PREC/s*ffbar(s)**(const115-2.0_PREC)
  
end function ffdp


