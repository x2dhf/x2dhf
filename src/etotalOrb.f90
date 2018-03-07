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
! ### etotalOrb ###

!   Calculates total energy in the case of exchange potentials being
!   kept on a disk during SCF cycles, i.e. by summing up contributions
!   due to a given orbital

module etotalOrb_m
  implicit none
contains
  subroutine etotalOrb (iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)
    use params
    use discret
    use commons8
    use util
    
    use blas_m
    use diffmu_m
    use diffnu_m
    use excont_m
    use exint_m
    use putin_m
    use putout_m

    implicit none
    integer :: ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iex,iex1,iorb,iorb1,iorb2,&
         isiex,isiex1,isiorb,isiorb1,isiorb2,isipot,isipot1,isipot2,isym,ngrid,nmut,nmut1,nmut2,norb2
    real (PREC) :: epscharge,oc,oc1,oc2,ocx1,ocx2,w,wdcoul,wex1,wex2,wndc,woneel

    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13


    !   an orbital is considerd occupied if its charge is greater than epscharge
    data epscharge /1.e-7_PREC/

    if (iorb.eq.0) then
       ictot=0
       vkt  =0.0_PREC
       vnt  =0.0_PREC
       evt  =0.0_PREC
       epott=0.0_PREC
       return
    endif

    engt(1)=eng(1)

    if (nel.lt.2) return

    ictot=ictot+1

    woneel=0.0_PREC
    wdcoul=0.0_PREC

    engt(1)=0.0_PREC

    !   calculate first contributions from one particle operators and
    !   coulomb potential contributions within the same shell

    if (iprint(80).ne.0) then
       write(*,*) 'orb        t        v(n)      t+v(n)     v(c-within shell)'
    endif

    norb2=norb*norb

    iborb=i1b(iorb)
    isiorb=i1si(iorb)
    nmut=i1mu(iorb)
    ibpot=i2b(iorb)
    isipot=i2si(iorb)
    oc=occ(iorb)

    isym=isymOrb(iorb)

    !   calculate derivatives over mu and ni variables by means of matrix
    !   multiplication

    call putin (nni,nmut,isym,psi(iborb),wk3)
    call diffnu (nmut,wk3,wk0,wk1,wk2)
    call putout (nni,nmut,wk1,wk0)
    call diffmu (nmut,wk3,wk2)
    call putout (nni,nmut,wk0,wk2)

    !   add derivatives over mu and ni

    call add (isiorb,wk0,wk1)

    !   add contribution from phi part of laplasian
    !   e enters the expression with minus sign which is already
    !   incorporated in e

    if (mm(iorb).ne.0) then
       w=dble(mm(iorb)*mm(iorb))
       call copy (isiorb,e,ione,wk0,ione)
       if (mm(iorb).ne.1) then
          call scal (isiorb,w,wk0,ione)
       endif
       call prod (isiorb,psi(iborb),wk0)
       call add (isiorb,wk0,wk1)
    endif

    call copy (isiorb,wk1,ione,wk2,ione)
    call prod  (isiorb,psi(iborb),wk2)
    w=dot(isiorb,wgt1,ione,wk2,ione)
    vk(iorb)=w

    call copy (isiorb,f0,ione,wk0,ione)
    call prod (isiorb,psi(iborb),wk0)

    call prod2 (isiorb,psi(iborb),wk0,wk2)
    w=dot(isiorb,wgt1,ione,wk2,ione)
    vn(iorb)=w

    call add (isiorb,wk0,wk1)
    call prod (isiorb,psi(iborb),wk1)

    w =dot(isiorb,wgt1,ione,wk1,ione)
    woneel=oc*w

    !   contribution from coulomb interaction within the same shell

    if (oc.gt.1.0_PREC) then
       if (isipot.ne.isiorb) then
          write(iout6,*) 'Coulomb potentials and orbitals have to defined on the same number of subgrids'
          stop 'etotalOrb'
       endif
       call copy (isipot,pot(ibpot),ione,wk2,ione)
       call prod  (isiorb,psi(iborb),wk2)
       call prod  (isiorb,psi(iborb),wk2)

       w=dot(isiorb,wgt2,ione,wk2,ione)
       wdcoul= oc*(oc-1.0_PREC)/2.0_PREC*w
    endif

    vkt=vkt+oc*vk(iorb)
    vnt=vnt+oc*vn(iorb)

    if (iprint(80).ne.0) then
       write(*,6000) iorb, vk(iorb),vn(iorb),vk(iorb)+vn(iorb),w
    endif

6000 format(' ',i4,2x,4e20.9)
    !6010 format(' ',i4,2x,4e20.9)
    !6012 format(' ',2i4,2x,4e20.9)

    !   contribution from coulomb and exchange interaction between shells

    wndc =0.0_PREC
    wex1 =0.0_PREC
    wex2 =0.0_PREC

    iorb1=iorb
    iborb1=i1b(iorb1)
    isiorb1=i1si(iorb1)
    nmut1=i1mu(iorb1)
    ibpot1=i2b(iorb1)
    isipot1=i2si(iorb1)
    oc1=occ(iorb1)

    iex1 =iorb1+iorb1*(iorb1-1)/2

    call copy (isiorb1,psi(iborb1),ione,wk0,ione)
    call prod  (isiorb1,psi(iborb1),wk0)

    !   calculate exchange interaction within pi, delta, etc. open shell
    if (mm(iorb1).gt.0.and.abs(oc1-one).gt.epscharge) then
       call exint (iorb1,ocx1)
       ibex=i3b(iex1)
       isiex1=i3si(iex1)

       ngrid=min(isiorb1,isiex1)
       call prod2  (ngrid,wk0,excp(ibex),wk1)
       w=dot (ngrid,wgt2,ione,wk1,ione)
       wex1=wex1+ocx1*w

       if (iprint(80).ne.0) then
          write(*,6020) iorb1,w
6020      format(' ','orb v(x-within shell) ',i4,2x,4e15.5)
       endif
    endif

    do iorb2=iorb1+1,norb
       iborb2=i1b(iorb2)
       isiorb2=i1si(iorb2)
       nmut2=i1mu(iorb2)
       ibpot2=i2b(iorb2)
       isipot2=i2si(iorb2)
       oc2=occ(iorb2)

       !      coulomb interaction between shells

       ngrid=min(isipot2,isiorb1)
       call prod2 (ngrid,pot(ibpot2),wk0,wk1)

       w=dot(ngrid,wgt2,ione,wk1,ione)
       wndc=wndc+oc1*oc2*w

       if (iprint(80).ne.0) then
          write(*,6022) iorb1,iorb2,w
6022      format(' ','orb1 orb2 v(c-inter shell) ',2i4,2x,4e15.5)
       endif

       iex=iorb1+iorb2*(iorb2-1)/2
       ibex=i3b(iex)
       isiex=i3si(iex)

       !      exchange interaction between shells (same lambda)

       call excont (iorb1,iorb2,ocx1,ocx2)

       ngrid=min(isiorb1,isiorb2,isiex)
       call prod2 (ngrid,excp(ibex),psi(iborb1),wk1)
       call prod  (ngrid,psi(iborb2),wk1)

       w=dot(ngrid,wgt2,ione,wk1,ione)
       wex1=wex1+ocx1*w

       if (iprint(80).ne.0) then
          write(*,6024) iorb1,iorb2,w
6024      format(' ','orb1 orb2 v(x-inter shell)a',2i4,2x,4e15.5)
       endif

       if (ilc(iex).gt.1) then
          call prod2 (ngrid,excp(ibex+ngrid),psi(iborb1),wk1)
          call prod  (ngrid,psi(iborb2),wk1)
          w=dot (ngrid,wgt2,ione,wk1,ione)
          wex2=wex2+ocx2*w

          if (iprint(80).ne.0) then
             write(*,6026) iorb1,iorb2,w
6026         format(' ','orb1 orb2 v(x-inter shell)b',2i4,2x,4e15.5)
          endif
       endif
    enddo

    evt=evt+woneel+wdcoul+wndc-wex1-wex2
    epott=epott+oc*vn(iorb)+wdcoul+wndc-wex1-wex2

    virrat=(epott+z1*z2/r)/vkt

    if (ictot.eq.norb) then
       iready=ictot
       etot=evt+z1*z2/r
       engt(1)=evt

       ictot=0
       vkt  =0.0_PREC
       vnt  =0.0_PREC
       evt  =0.0_PREC
       epott=0.0_PREC
    endif

    if (iprint(81).ne.0) then
       write(*,*) 'total energy: ', evt
       write(*,*) 'contributions: one-electron     coulomb    exchange  '
       write(*,*)  woneel,wdcoul+wndc,wex1+wex2
    endif

  end subroutine etotalOrb
end module etotalOrb_m
