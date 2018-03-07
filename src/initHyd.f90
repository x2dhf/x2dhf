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
! c ### initHyd ###
! c
! c     This routine initializes molecular orbitals as linear combinations
! c     of hydrogenic functions on centres A and B.  In the case of HF or
! c     HFS calculations Coulomb (exchange) potentials are approximated as
! c     a linear combination of Thomas-Fermi (1/r) potentials of the two
! c     centres; if method OED is chosen the potential functions are set
! c     to zero.

subroutine initHyd (psi,pot,excp,f2,f4,wgt2,wk0)
  use params
  use commons8
  use blas_m

  implicit none

  integer :: i1beg,iorb,l1,l2,m1,m2,n1,n2,ngorb
  real (PREC) ez1,ez2

  real (PREC), dimension(*) :: psi,pot
  real (PREC), dimension(*) :: excp
  real (PREC), dimension(*) :: f2,f4,wgt2
  real (PREC), dimension(*) :: wk0


  !     Initialization of molecular orbitals

  print *,'... initializing orbitals from hydrogenic functions ...'

!     loop over orbitals

!      do iorb=1,norbt
      do iorb=1,norb
         if (inhyd(iorb).eq.0) goto 100
         i1beg=i1b(iorb)
         ngorb=i1si(iorb)
         call zeroArray(ngorb,psi(i1beg))

         if     (co1(iorb).ne.zero) then
            n1=mgx(1,iorb)
            l1=mgx(2,iorb)
            m1=mgx(3,iorb)
            ez1=eza1(iorb)
            call hydrogenOrbA(wk0,n1,l1,m1,ez1)
            call axpy(ngorb,co1(iorb),wk0,ione,psi(i1beg),ione)
         endif

         if (co2(iorb).ne.zero) then
            n2=mgx(4,iorb)
            l2=mgx(5,iorb)
            m2=mgx(6,iorb)
            ez2=eza2(iorb)
            call hydrogenOrbB(wk0,n2,l2,m2,ez2)
            call axpy (ngorb,co2(iorb),wk0,ione,psi(i1beg),ione)
         endif
 100     continue
      enddo

!01114 format(/1x,'  orbital        norm       (# primitive bf)')
!01115 format(1x,i2,1x,a5,e20.12,' (',i4,')')

!  initialize Coulomb and exchange potentials

      call initPot(psi,pot,excp,f2,f4,wk0)

    end subroutine initHyd
