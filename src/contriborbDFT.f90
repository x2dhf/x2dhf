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
! ### contriborbDFT ###  
!
!     Calculates various contributions to total DFT energy due to a given
!     orbital

subroutine contriborbDFT(iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
  use params
  use commons8

  implicit none
  integer :: iorb
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3

! FIXME!!! make EaDFT into contribDFTorb with detailed printouts 
  call EaDFT (iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)

end subroutine contriborbDFT
