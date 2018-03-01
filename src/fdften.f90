! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fdften ###

!     Calculates coefficients appearing in LDA energy formulea.

function fdften (alpha)
  use params
  use commons8

  implicit none

  real (PREC) :: fdften
  real (PREC) :: alpha,const13,const32,const34,const94,const98
  parameter (const13=1.0_PREC/3.0_PREC,const32=3.0_PREC/2.0_PREC,const34=3.0_PREC/4.0_PREC, &
       const94=9.0_PREC/4.0_PREC,const98=9.0_PREC/8.0_PREC)
  
  
  if (idftex.eq.1.or.idftex.eq.2) then
     !        L.Laaksonen, P.Pyykko, D.Sundholm (Int. J.Quantum Chem. 27 (1985) 601)
     !        additional factor 2**(1/3) is due to the way the total density is calculated
     !        (see fldapot)
     
     !        Becke (idftex=2)
     fdften = -alpha*const98*(three/pii)**const13*two**(const13)
  elseif (idftex.eq.3) then
     !        generalized gradient approximation (GGA) (Perdew & Wang 86 )
     fdften = -const34*(three/pii)**const13
     
  elseif (idftex.eq.4) then
     !        generalized gradient approximation (GGA) (Perdew & Wang 91 )
     fdften = -const34*(three/pii)**const13
  endif
  
end function fdften
