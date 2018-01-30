! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1992 F.A. Parpia and I.P. Grant                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### zgsz2g ###
!
!     Evaluate the model HF potential according to the formula derived by 
!     Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z2 centre
!     with Z2 modified according to the Gauss nucleus model

function zgsz2g(i,j)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,j,izz2
  real (PREC) :: zgsz2g
  real (PREC) :: hc,ri,zz2t
  real (PREC), external :: zz2g

  ri = r*(vxi(i)-veta(j))/two
  izz2=nint(z2)
  if (izz2.eq.0) then
     zgsz2g=0.0_PREC
     return
  endif
  
  if (izz2.eq.1) then
     hc=dgsz(izz2)*(z2)**0.4_PREC
     zgsz2g=(z2)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
     return
  endif
  
  zz2t=zz2g(i,j)
  if (zz2t.eq.zero) then
     zgsz2g=zero
     return
  endif
  
  hc=dgsz(izz2)*(zz2t-one)**0.4_PREC
  zgsz2g=(zz2t-one)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
  
end function zgsz2g
