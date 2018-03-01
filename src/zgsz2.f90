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
! ### zgsz2 ###
!      Evaluate the model HF potential according to the formula derived by
!      Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z1 centre

function zgsz2(i,j)
  use params
  use discret
  use commons8

  implicit none

  integer :: i,j,izz2
  real (PREC) :: zgsz2
  real (PREC) hc,ri

  ri = r*(vxi(i)-veta(j))/two
  izz2=nint(z2)
  if (izz2.eq.0) then
     zgsz2=0.0_PREC
     return
  endif

  if (izz2.eq.1) then
     hc=dgsz(izz2)*(z2)**0.4_PREC
     zgsz2=(z2)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
     return
  endif

  hc=dgsz(izz2)*(z2-one)**0.4_PREC
  zgsz2=(z2-one)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one

end function zgsz2


