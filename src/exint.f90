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
! ### exint ###
!
!     Calculates the number of exchange integrals within one open nonsigma
!     shell.

subroutine exint (iorb1,ox1)
  use params
  use commons8

  implicit none
  integer :: i1,i2,iorb1,ip1
  real (PREC) :: ox1

  character*8 :: alpha,beta

  data alpha,beta /'+','-'/

  ox1=0.0_PREC
  ip1=4*(iorb1-1)

  !    interaction within nonsigma shell

  do i1=1,2
     if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 10
     do i2=3,4
        if(spin(ip1+i2).ne.alpha.and.spin(ip1+i2).ne.beta) goto 12
        if(spin(ip1+i1).eq.spin(ip1+i2)) ox1=ox1+1.0_PREC
00012   continue
     enddo
00010 continue
  enddo

end subroutine exint
