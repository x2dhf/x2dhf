! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### scmc ###

!     Calculates self-consistent multiplicative constant (see Eq. 3 in
!     Karasiev, Ludenia, eq.3 PR 64 (2002) 062510)

module scmc_m
  implicit none
contains
  subroutine scmc (psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)
    use params
    use commons8
    use util

    use blas_m
    use dftex_m
    use excont_m
    use exint_m
    use rfdexch_m

    implicit none
    integer ::         ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iex,iex1,iorb,iorb1,iorb2,ipe1,ipe2,&
         isiorb,isiorb1,isiorb2,isipot,isipot1,isipot2,isiex,isiex1,ngrid,nmut,nmut1,nmut2

    real (PREC) :: ehfex,eps,oc,oc1,oc2,ocx1,ocx2,w,wdcoul,wex1,wex2
    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,&
         wk7,wk8,wk9,wk10,wk11,wk12,wk13

    data eps /1.e-7_PREC/

    if (nel.lt.2) return

    ! FIXME
    !   contribution from coulomb interaction within the same shell

    wdcoul=0.0_PREC

    do iorb=1,norb
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)

       if (oc.gt.1.0_PREC) then
          call copy (isipot,pot(ibpot),ione,wk2,ione)
          call prod  (isiorb,psi(iborb),wk2)
          call prod  (isiorb,psi(iborb),wk2)

          w=dot(isiorb,wgt2,ione,wk2,ione)
          !          wdcoul=wdcoul+oc*(oc-1.0_PREC)/2.0_PREC*w
          wdcoul=wdcoul+oc/2.0_PREC*w
       endif
    enddo

    !   contribution from coulomb and exchange interaction between shells

    wex1 =0.0_PREC
    wex2 =0.0_PREC

    do iorb1=1,norb

       !      retrieve from a disk file exchange potentials involving iorb1

       if (iform.eq.0.or.iform.eq.2) then
          call  rfdexch(iorb1,excp)
       endif

       iborb1=i1b(iorb1)
       isiorb1=i1si(iorb1)
       nmut1=i1mu(iorb1)
       ibpot1=i2b(iorb1)
       isipot1=i2si(iorb1)
       oc1=occ(iorb1)

       ipe1=mgx(6,iorb1)
       iex1 =iorb1+iorb1*(iorb1-1)/2

       call copy (isiorb1,psi(iborb1),ione,wk0,ione)
       call prod  (isiorb1,psi(iborb1),wk0)

       !      calculate exchange interaction within pi, delta, etc. open shell
       ! FIXME
       if (ipe1.gt.0.and.abs(oc1-1.0_PREC).gt.eps) then
          call exint (iorb1,ocx1)
          ibex=i3b(iex1)
          isiex1=i3si(iex1)

          ngrid=min(isiorb1,isiex1)
          call prod2  (ngrid,wk0,excp(ibex),wk1)
          w=dot (ngrid,wgt2,ione,wk1,ione)
          wex1=wex1+ocx1*w

          if (iprint(82).ne.0) then
             write(*,6020) iorb1,w
6020         format(' ','orb v(x-within shell) ',i4,2x,4e15.5)
          endif

          if (iprint(83).ne.0) then
             write(*,7020) iorn(iorb1),bond(iorb1),gut(iorb1),ipe1,ocx1
7020         format(' ','v(x-within shell)  ',i4,1x,a8,a1,i4,f6.2)
          endif
       endif

       do iorb2=iorb1+1,norb
          iborb2=i1b(iorb2)
          isiorb2=i1si(iorb2)
          nmut2=i1mu(iorb2)
          ibpot2=i2b(iorb2)
          isipot2=i2si(iorb2)
          oc2=occ(iorb2)

          ipe2=mgx(6,iorb2)
          iex=iorb1+iorb2*(iorb2-1)/2
          ibex=i3b(iex)
          isiex=i3si(iex)

          !         exchange interaction between shells (same lambda)

          call excont (iorb1,iorb2,ocx1,ocx2)
          ngrid=min(isiorb1,isiorb2,isiex)
          call prod2 (ngrid,excp(ibex),psi(iborb1),wk1)
          call prod  (ngrid,psi(iborb2),wk1)
          w=dot(ngrid,wgt2,ione,wk1,ione)
          wex1=wex1+ocx1*w
          if (iprint(82).ne.0) then
             print *,'ibex,wex1,ocx1,w',ibex,wex1,ocx1,w
          endif
          !         exchange interaction between shells (different lambda)

          if (ilc(iex).gt.1) then
             call prod2 (ngrid,excp(ibex+ngrid),psi(iborb1),wk1)
             call prod  (ngrid,psi(iborb2),wk1)
             w=dot (ngrid,wgt2,ione,wk1,ione)
             wex2=wex2+ocx2*w
             if (iprint(83).ne.0) then
                print *,'wex2,ocx2,w',wex2,ocx2,w
             endif
          endif
       enddo
    enddo

    ! FIXME
    ehfex=wdcoul+wex1+wex2

    alphaf=two/three
    !   exchange energy from DFT functionals
    call dftex (psi,pot,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)

    alphaf=abs(ehfex/edftex)*two/three

    if (iprint(84).ne.0) then
       write(*,'("scmc: ehfex,edftex,alphaf",3e15.6)') ehfex,edftex,alphaf
    endif
    if (iprint(85).ne.0) then
       write(*,'("scmc: wdcoul,wex1,wex2",3e15.6)') wdcoul,wex1,wex2
    endif

  end subroutine scmc
end module scmc_m
