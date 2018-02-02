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
! ### zgsz1g ###
!
!     Evaluate the model HF potential according to the formula derived by
!     Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z1 centre
!     with Z1 modified according to the Gauss nucleus model

function zgsz1g(i,j)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,j,izz1
  real (PREC) :: zgsz1g
  real (PREC) :: hc,ri,zz1t
  real (PREC), external :: zz1g

  ri = r*(vxi(i)+veta(j))/two
  izz1=nint(z1)
  if (izz1.eq.0) then
     zgsz1g=0.0_PREC
     return
  endif

  if (izz1.eq.1) then
     hc=dgsz(izz1)*(z1)**0.4_PREC
     zgsz1g=(z1)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one
     return
  endif

  zz1t=zz1g(i,j)
  if (zz1t.eq.zero) then
     zgsz1g=zero
     return
  endif

  hc=dgsz(izz1)*(zz1t-one)**0.4_PREC
  zgsz1g=(zz1t-one)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one

end function zgsz1g
