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
! ### Ea ###

!     Calculates the eigenvalue of h operator as <orb_a|h|orb_a>.
!     Orbitals are assumed to be normalized.

module Ea_m
  implicit none
contains
  subroutine Ea(iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
    use params
    use discret
    use scf
    use commons8
    use util

    use blas_m
    use diffmu_m
    use diffnu_m
    use pmtx_m
    use putin_m
    use putout_m

    implicit none
    integer :: i,iorb,i1beg,i1beg1,i2beg,i2beg1,i3beg,ihc,iorb1,ipc,isym,kex, &
         nmut,ngex,ngorb,ngorb1,ngpot,ngpot1,ngrid,ngrid2
    real (PREC) :: eshift,oc,w,wkin,wnucl,woneel,wtwoel
    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3

    wtwoel=zero

    i1beg=i1b(iorb)
    i2beg=i2b(iorb)
    nmut=i1mu(iorb)
    ngorb=i1si(iorb)
    ngpot=i2si(iorb)

    if (ngorb.ne.ngpot) then
       write(*,*) 'Orbitals and corresponding Coulomb potentials have to be defined on the same number of grid points'
       stop 'Ea'
    endif

    ipc=iorb+(iorb-1)*norb
    engo(ipc)=zero
    isym=isymOrb(iorb)

    !     calculate derivatives over mu and ni

    call putin (nni,nmut,isym,psi(i1beg),wk3)
    call diffnu (nmut,wk3,wk0,wk1,wk2)
    call putout (nni,nmut,wk1,wk0)

    call diffmu (nmut,wk3,wk2)
    call putout (nni,nmut,wk0,wk2)


    !     add contribution from derivatives over mu and ni

    call add (ngorb,wk0,wk1)

    if (iprint(60).ne.0) then
       call copy (ngorb,wk1,ione,wk2,ione)
       call prod  (ngorb,psi(i1beg),wk2)
       wkin=dot(ngorb,wgt1,ione,wk2,ione)
    endif

    call copy (ngorb,f0,ione,wk0,ione)

    if (mm(iorb).ne.0) then
       !        nuclear energy for non-sigma orbitals contains contribution
       !        from e term (in toten this term is correctly added to the
       !        kinetic energy contribution); e enters the expression with
       !        minus sign which is already incorporated in e

       w=dble(mm(iorb)*mm(iorb))
       call axpy (ngorb,w,e,ione,wk0,ione)
    endif

    if (iprint(60).ne.0) then
       call prod2 (ngorb,psi(i1beg),wk0,wk2)
       call prod  (ngorb,psi(i1beg),wk2)
       wnucl=dot(ngorb,wgt1,ione,wk2,ione)
    endif

    call proda (ngorb,psi(i1beg),wk0,wk1)
    call prod (ngorb,psi(i1beg),wk1)

    woneel=dot(ngorb,wgt1,ione,wk1,ione)

    if(iprint(60).ne.0) then
       write(*,'("Ea: ",i4,1x,a8,a1,"energy contributions")') iorn(iorb),bond(iorb),gut(iorb)
       write(*,'("    occupation #   ",e26.16)') occ(iorb)
       write(*,'("    kinetic        ",e26.16)') wkin
       write(*,'("    nuclear        ",e26.16)') wnucl
       write(*,'("    kinetic+nuclear",e26.16)') woneel
       !         write(*,'("    kinetic+nuclear",3e26.16)') wkin+wnucl,woneel
    endif

    if(iprint(62).ne.0) then
       write(*,*) 'Ea: wk0'
       call pmtx (nni,nmu(1),wk0,ione,ione,incrni,incrmu)
       write(*,*) 'Ea: wk1'
       call pmtx (nni,nmu(1),wk1,ione,ione,incrni,incrmu)
    endif

    if (nel.gt.1) then
       do i=1,mxsize
          wk2(i)=zero
       enddo

       if (norb.eq.1) then

          !           contribution from coulomb potential of one sigma orbital
          !           14/12/02
          !	    call copy (ngpot,pot(i2beg),ione,wk2,ione)
          oc=occ(norb)-one
          call axpy (ngpot,oc,pot(i2beg),ione,wk2,ione)
          call prod  (ngorb,psi(i1beg),wk2)
       else

          !           add contributions from coulomb and exchange potentials

          do iorb1=1,norb

             !              it is asumed that exchange potentials involving iorb1
             !              has already been retrieved from disk

             i1beg1=i1b(iorb1)
             i2beg1=i2b(iorb1)
             ngorb1=i1si(iorb1)
             ngpot1=i2si(iorb1)
             oc=occ(iorb1)
             kex=iorb+norb*(iorb1-1)

             if (iorb.le.iorb1) then
                ihc=iorb+iorb1*(iorb1-1)/2
             else
                ihc=iorb1+iorb*(iorb-1)/2
             endif

             i3beg=i3b(ihc)
             ngex =i3si(ihc)
             ngrid=min(ngorb,ngpot1)
             do i=1,mxsize
                wk0(i)=zero
             enddo

             call copy (ngpot1,pot(i2beg1),ione,wk0,ione)
             call prod  (ngrid,psi(i1beg),wk0)

             if (iorb.eq.iorb1) oc=oc-1.0_PREC
             if (oc.ne.1.0_PREC) then
                call scal (ngrid,oc,wk0,ione)
             endif

             if (iorb1.ne.iorb)  then
                ngrid2=min(ngorb1,ngex,ngorb)
                oc=gec(kex)
                call prodas (ngrid2,-oc,psi(i1beg1),excp(i3beg),wk0)
                if (ilc(ihc).gt.1) then
                   oc=gec(kex+norb*norb)
                   call prodas (ngrid2,-oc,psi(i1beg1),excp(i3beg+ngex),wk0)
                endif
             else
                if ((mm(iorb).gt.0).and.(ilc(ihc).gt.0)) then
                   ngrid2=min(ngorb1,ngex,ngorb)
                   oc=gec(kex)
                   call prodas(ngrid2,-oc,psi(i1beg1),excp(i3beg),wk0)
                endif
             endif

             if (iorb1.eq.1) then
                call copy (ngrid,wk0,ione,wk2,ione)
             else
                call add   (ngrid,wk0,wk2)
             endif
          enddo
       endif

       !        to complete the integrand wk2 has to be multiplied by psi(i1beg)

       call prod (ngorb,psi(i1beg),wk2)
       wtwoel=dot(ngorb,wgt2,ione,wk2,ione)
    endif

    eshift=zero
    engo(ipc) = woneel+wtwoel
    engo(ipc)=engo(ipc)+eshift
    eng(iorb) = engo(ipc)

    if(iprint(64).ne.0) then
       write(*,'("iorb,woneel,wtwoel,woneel+wtwoel",i6,3e22.14)') iorb,woneel,wtwoel,woneel+wtwoel
    endif

    if(iprint(65).ne.0) then
       write(*,*) 'wtwoel', wtwoel
       write(*,*) 'Ea: wk0'
       call pmtx (nni,nmu(1),wk0,ione,ione,incrni,incrmu)
       write(*,*) 'Ea: wk2'
       call pmtx (nni,nmu(1),wk2,ione,ione,incrni,incrmu)
    endif

    if (iprint(66).ne.0) then
       write(*,'("Ea: woneel, wtwoel, eng ",3e26.16)') woneel,wtwoel,eng(iorb)
    endif

  end subroutine Ea
end module Ea_m
