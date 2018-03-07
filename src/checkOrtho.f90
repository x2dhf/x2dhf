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
! ### checkOrtho ###

!     Checks orthogonalization of a given orbital

module checkOrtho_m
  implicit none
contains
  subroutine checkOrtho (iorb1,psi,f4,wgt2,wk0)
    use params
    use commons8
    use util
    use blas_m

    implicit none
    integer :: ibeg1,ibeg3,iorb1,iorb2,iorb3,jor,jor1,ngrid

    integer, dimension(maxorb) :: istp
    real (PREC), dimension(maxorb) :: ovla,ovla1
    real (PREC), dimension(*) :: psi,f4,wgt2,wk0

    jor=0
    jor1=0
    ibeg1 = i1b (iorb1)

    !     orbitals are stored in reverse order
    do iorb2=1,norb
       iorb3=norb-iorb2+1
       if (iorb3.le.iorb1) go to 101

       ngrid= min (i1si(iorb1),i1si(iorb3))
       if (mgx(6,iorb3).ne.mgx(6,iorb1)) go to 101
       if (ige(iorb3).ne.ige(iorb1)) go to 101

       jor=jor+1
       ibeg3= i1b (iorb3 )
       istp(jor)=ibeg3

       call prod2 (ngrid,psi(ibeg1),psi(ibeg3),wk0)
       call prod  (ngrid,f4,wk0)

       ovla(jor)=dot (ngrid,wgt2,ione,wk0,ione)

       jor1=jor1+1
       ovla1(jor1)=abs(ovla(jor))

       !        print overlap integrals
       if (iprint(20).ne.0) then
          write(*,60) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb3),bond(iorb3),gut(iorb3),ovla(jor)
60        format(4x,' <',i2,1x,a5,1x,a1,1x,'|',1x,i2,1x,a5,1x,a1,1x,'> = ',e10.2)
       endif

       !        print overlap integrals only for pairs of orbitals having
       !        non-zero off-diagonal Lagrange multipliers
       if (iprint(21).ne.0.and.nlm(iorb1,iorb3).ne.0) then
          write(*,60) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb3),bond(iorb3),gut(iorb3),ovla(jor)
       endif
101    continue
    enddo

  end subroutine checkOrtho
end module checkOrtho_m
