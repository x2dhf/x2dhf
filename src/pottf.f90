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
! ### pottf ###
!
!     Calculates a Thomas-Fermi potential (a modified version of dentfa by Desclaux)

function pottf(r,ek,qk,dz,ch,slim)
  use params

  implicit none
  real (PREC) :: pottf
  real (PREC) :: ch,dz,ek,qk,r,slim,t,w

  pottf=0.0_PREC
  if((dz+ch).lt.1.e-04_PREC)  return

  if(abs(qk+ek).lt.1.e-10_PREC) return

  w=r*(qk+ek)/2.0_PREC*(dz+ch)**(1.0_PREC/3.0_PREC)
  w=sqrt(w/0.885340_PREC)
  t=w*(0.6011200_PREC*w+1.8106100_PREC)+1.0_PREC
  w=w*(w*(w*(w*(0.0479300_PREC*w+0.2146500_PREC)+0.7711200_PREC)+1.3951500_PREC)+1.8106100_PREC)+1.0_PREC
  pottf=slim*(1.0_PREC-(t/w)*(t/w))*2.0_PREC/(r*(qk+ek))

end function pottf
