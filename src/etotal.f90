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
! ### etotal ###

!    Calculates total HF energy

subroutine etotal (psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)
  use params
  use discret
  use commons8

  implicit none
  integer :: ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iex,iex1,iorb,iorb1,iorb2,isiex,isiex1,isiorb,&
       isiorb1,isiorb2,isipot,isipot1,isipot2,isym,ngrid,nmut,nmut1,nmut2,norb2

  real (PREC) :: epscharge,etsum,oc,oc1,oc2,ocx2,ocx1,w,wcouldft,wdcoul,wex1,wex2,wndc,woneel
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6, &
       wk7,wk8,wk9,wk10,wk11,wk12,wk13
  real (PREC), external :: dot

  data epscharge /1.e-7_PREC/

  engt(1)=eng(1)

  call prepwexch

  engt(1)=0.0_PREC

  !     calculate first contributions from one particle operators and
  !     coulomb potential contributions within the same shell

  norb2=norb*norb
  woneel=0.0_PREC
  wdcoul=0.0_PREC

  vkt=0.0_PREC
  vnt=0.0_PREC

  etsum=zero
  wcouldft=0.0_PREC
  do iorb=1,norb
     iborb=i1b(iorb)
     isiorb=i1si(iorb)
     nmut=i1mu(iorb)
     ibpot=i2b(iorb)
     isipot=i2si(iorb)
     oc=occ(iorb)

     isym=isymOrb(iorb)

     !        calculate derivatives over mu and ni variables by means of matrix
     !        multiplication

     call putin (nni,nmut,isym,psi(iborb),wk3)
     call diffnu (nmut,wk3,wk0,wk1,wk2)
     call putout (nni,nmut,wk1,wk0)
     call diffmu (nmut,wk3,wk2)
     call putout (nni,nmut,wk0,wk2)

     !        add derivatives over mu and ni

     call add (isiorb,wk0,wk1)

     !        add contribution from phi part of laplasian
     !        e enters the expression with minus sign which is already
     !        incorporated in e

     if (mm(iorb).ne.0) then
        w=dble(mm(iorb)*mm(iorb))
        call copy (isiorb,e,ione,wk0,ione)
        if (mm(iorb).ne.1) then
           call scal (isiorb,w,wk0,ione)
        endif
        call prod  (isiorb,psi(iborb),wk0)
        call add   (isiorb,wk0,wk1)
     endif

     call copy (isiorb,wk1,ione,wk2,ione)
     call prod  (isiorb,psi(iborb),wk2)
     w=dot(isiorb,wgt1,ione,wk2,ione)
     vk(iorb)=w

     call copy (isiorb,f0,ione,wk0,ione)
     call prod  (isiorb,psi(iborb),wk0)

     call prod2 (isiorb,psi(iborb),wk0,wk2)
     w=dot(isiorb,wgt1,ione,wk2,ione)
     vn(iorb)=w

     call add (isiorb,wk0,wk1)
     call prod (isiorb,psi(iborb),wk1)

     w =dot(isiorb,wgt1,ione,wk1,ione)
     woneel=woneel+oc*w

     vkt=vkt+oc*vk(iorb)
     vnt=vnt+oc*vn(iorb)

     if (iprint(78).ne.0) then
        etsum=etsum+oc*vk(iorb)
        write(*,7028) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb),vk(iorb),oc,etsum
7028    format('<',i4,1x,a5,a1,'| T |',i4,1x,a5,a1,' >',26x,d25.16, f8.2, d25.16)

        etsum=etsum+oc*vn(iorb)
        write(*,7030) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb),vn(iorb),oc,etsum
7030    format('<',i4,1x,a5,a1,'| V |',i4,1x,a5,a1,' >',26x,d25.16, f8.2, d25.16)
     endif

     !          if (nel.lt.2) then
     !             vkt=vkt+oc*vk(iorb)
     !             vnt=vnt+oc*vn(iorb)
     !             evt=woneel
     !             epott=vnt
     !             virrat=(epott+z1*z2/r)/vkt
     !             etot=evt+z1*z2/r
     !          endif

     !        contribution from Coulomb interaction within the same shell

     if (oc.gt.one) then
        if (isipot.ne.isiorb) then
           write(iout6,*) 'coulomb potentials and orbitals have to be defined on the same number of subgrids'
           stop 'etotal'
        endif
        call copy (isipot,pot(ibpot),ione,wk2,ione)
        call prod  (isiorb,psi(iborb),wk2)
        call prod  (isiorb,psi(iborb),wk2)

        w=dot(isiorb,wgt2,ione,wk2,ione)
        wdcoul=wdcoul+oc*(oc-one)/two*w

        !  To calculate DFT exchange energy the Coulomb energy must include
        !  J_{ii}=K_{ii} terms

        wcouldft=wcouldft+oc*oc/2.0_PREC*w


        if (iprint(78).ne.0) then
           etsum=etsum+oc*(oc-one)/two*w
           write(*,7031) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb), &
                iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb),w,oc*(oc-one)/two,etsum
7031       format('<',i4,1x,a5,a1,'| J1(',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',d25.16, f8.2, d25.16)
        endif

     endif
  enddo

! 6000 format(' ',i4,2x,4e22.14)
! 6010 format(' ',i4,2x,4e22.14)
! 6012 format(' ',2i4,2x,4e22.14)

!     contribution from coulomb and exchange interaction between shells

  wndc =0.0_PREC
  wex1 =0.0_PREC
  wex2 =0.0_PREC

  do iorb1=1,norb

     !        retrieve from a disk file exchange potentials involving iorb1

     if (iform.eq.0.or.iform.eq.2) then
        call  rfdexch(iorb1,excp)
     endif

     iborb1=i1b(iorb1)
     isiorb1=i1si(iorb1)
     nmut1=i1mu(iorb1)
     ibpot1=i2b(iorb1)
     isipot1=i2si(iorb1)
     oc1=occ(iorb1)

     iex1 =iorb1+iorb1*(iorb1-1)/2

     call copy (isiorb1,psi(iborb1),ione,wk0,ione)
     call prod  (isiorb1,psi(iborb1),wk0)

     !        calculate exchange interaction within pi, delta, etc. shell

     if (mm(iorb1).gt.0.and.abs(oc1-one).gt.epscharge) then
        call exint (iorb1,ocx1)
        ibex=i3b(iex1)
        isiex1=i3si(iex1)

        ngrid=min(isiorb1,isiex1)
        call prod2  (ngrid,wk0,excp(ibex),wk1)
        w=dot (ngrid,wgt2,ione,wk1,ione)
        wex1=wex1+ocx1*w

        if (iprint(78).ne.0) then
           etsum=etsum-ocx1*w
           write(*,7032) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1), &
                iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1),w,-ocx1,etsum
7032       format('<',i4,1x,a5,a1,'| K (',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',d25.16, f8.2, d25.16)
        endif
     endif

     do iorb2=iorb1+1,norb
        iborb2=i1b(iorb2)
        isiorb2=i1si(iorb2)
        nmut2=i1mu(iorb2)
        ibpot2=i2b(iorb2)
        isipot2=i2si(iorb2)
        oc2=occ(iorb2)

        !           Coulomb interaction between shells

        ngrid=min(isipot2,isiorb1)
        call prod2 (ngrid,pot(ibpot2),wk0,wk1)

        w=dot(ngrid,wgt2,ione,wk1,ione)
        wndc=wndc+oc1*oc2*w
        wcouldft=wcouldft+oc1*oc2*w

        if (iprint(78).ne.0) then
           etsum=etsum+oc1*oc2*w
           write(*,7034) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb2),bond(iorb2),gut(iorb2), &
                iorn(iorb2),bond(iorb2),gut(iorb2),iorn(iorb1),bond(iorb1),gut(iorb1),w,oc1*oc2,etsum
7034       format('<',i4,1x,a5,a1,'| J (',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',d25.16, f8.2, d25.16)
        endif


        iex=iorb1+iorb2*(iorb2-1)/2
        ibex=i3b(iex)
        isiex=i3si(iex)

        !           exchange interaction between shells (same lambda)

        call excont (iorb1,iorb2,ocx1,ocx2)

        ngrid=min(isiorb1,isiorb2,isiex)
        call prod2 (ngrid,excp(ibex),psi(iborb1),wk1)
        call prod  (ngrid,psi(iborb2),wk1)

        w=dot(ngrid,wgt2,ione,wk1,ione)
        wex1=wex1+ocx1*w

        if (iprint(78).ne.0) then
           etsum=etsum-ocx1*w
           write(*,7036) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1), &
                iorn(iorb2),bond(iorb2),gut(iorb2),iorn(iorb2),bond(iorb2),gut(iorb2),w,-ocx1,etsum
7036       format('<',i4,1x,a5,a1,'| K (',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',d25.16, f8.2, d25.16)
        endif

        !           exchange interaction between shells (different lambda)

        if (ilc(iex).gt.1) then
           call prod2 (ngrid,excp(ibex+ngrid),psi(iborb1),wk1)
           call prod  (ngrid,psi(iborb2),wk1)
           w=dot (ngrid,wgt2,ione,wk1,ione)
           wex2=wex2+ocx2*w

           if (iprint(78).ne.0) then
              etsum=etsum-ocx2*w
              write(*,7038) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1), &
                   iorn(iorb2),bond(iorb2),gut(iorb2),iorn(iorb2),bond(iorb2),gut(iorb2),w,-ocx2,etsum
7038          format('<',i4,1x,a5,a1,'| K1(',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',d25.16, f8.2, d25.16)
           endif
        endif
     enddo
  enddo

  evt=woneel+wdcoul+wndc-wex1-wex2
  epott=vnt+wdcoul+wndc-wex1-wex2
  virrat=(epott+z1*z2/r)/vkt
  etot=evt+z1*z2/r

  enkin=vkt
  ennucel=vnt
  encoul=wdcoul+wndc
  enexch=-wex1-wex2
  entot=enkin+ennucel+encoul+enexch+z1*z2/r

  !     to be able to compare HF and DFT exchange energies these
  !     contributions have to be calculated as
  encouldft=wcouldft
  enexchdft=evt-woneel-wcouldft

  engt(1)=evt

  etsum=zero
  if (iprint(79).ne.0) then
     write(*,*)
     write(*,'(" total energy contributions and their sum: ")')
     etsum=etsum+woneel
     write(*,'("  one-electron:             ",2d25.16)') woneel,etsum
     etsum=etsum+wdcoul+wndc
     write(*,'("  Coulomb:                  ",2d25.16)') wdcoul+wndc,etsum
     etsum=etsum-wex1-wex2
     write(*,'("  exchange:                 ",2d25.16)') -wex1-wex2,etsum
     write(*,'("  two-electron:             ",2d25.16)') wdcoul+wndc-wex1-wex2
     write(*,'("  nuclear repulsion:        ",2d25.16)') z1*z2/r,etsum+z1*z2/r

     write(*,*)
     write(*,'("  Coulomb (DFT):            ",2d25.16)') encouldft
     write(*,'("  exchange (DFT):           ",2d25.16)') enexchdft
     write(*,'("  two-electron (DFT)        ",2d25.16)') encouldft+enexchdft
     write(*,*)

  endif

end subroutine etotal

