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
! ### zgsz1 ###
!
!     Evaluate the model HF potential according to the formula derived by
!     Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z1 centre

module zgsz1_m
  implicit none
contains
  function zgsz1(i,j)
    use params
    use discret
    use commons8

    implicit none
    integer :: i,j,izz1
    real (PREC) :: zgsz1
    real (PREC) hc,ri

    ri = r*(vxi(i)+veta(j))/two
    izz1=nint(z1)
    if (izz1.eq.0) then
       zgsz1=0.0_PREC
       return
    endif

    if (izz1.eq.1) then
       hc=dgsz(izz1)*(z1)**0.4_PREC
       zgsz1=(z1)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one
       return
    endif

    hc=dgsz(izz1)*(z1-one)**0.4_PREC
    zgsz1=(z1-one)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one

  end function zgsz1
end module zgsz1_m
