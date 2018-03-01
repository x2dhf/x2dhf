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
! c ### contrib ###
! c
! c     Calculates various contributions to total energy

subroutine contrib (psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
  use params
  use discret
  use commons8

  implicit none
  integer :: ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iex,iex1,iorb,iorb1,iorb2,ipe,ipe1,ipe2,&
        isiex,isiex1,isiorb,isiorb1,isiorb2,isipot,isipot1,isipot2,ngrid,nmut,nmut1,nmut2,norb2
  real (PREC) :: eps,oc,oc1,oc2,ocx1,ocx2,w,wdcoul,wex1,wex2,wndc,woneel
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3
  real (PREC),external :: dot

  data eps /1.e-7_PREC/

  engt(1)=eng(1)

  call prepwexch

  engt(1)=zero

  !   calculate first contributions from one particle operators and
  !   Coulomb potential contributions within the same shell

  norb2=norb*norb
  woneel=zero
  wdcoul=zero

  vkt=zero
  vnt=zero

  do iorb=1,norb
     iborb=i1b(iorb)
     isiorb=i1si(iorb)
     nmut=i1mu(iorb)
     ibpot=i2b(iorb)
     isipot=i2si(iorb)
     oc=occ(iorb)

     !      calculate derivatives over mu and ni variables by means of matrix
     !      multiplication

     call putin (nni,nmut,isymOrb(iorb),psi(iborb),wk3)
     call diffnu (nmut,wk3,wk0,wk1,wk2)
     call putout (nni,nmut,wk1,wk0)
     call diffmu (nmut,wk3,wk2)
     call putout (nni,nmut,wk0,wk2)

     !      add derivatives over mu and ni

     call add (isiorb,wk0,wk1)

     !      add contribution from phi part of laplacian
     !      e enters the expression with minus sign which is already
     !      incorporated in e

     if (ipe.ne.0) then
        w=dble(ipe*ipe)
        call copy (isiorb,e,ione,wk0,ione)
        if (ipe.ne.1) then
           call scal (isiorb,w,wk0,ione)
        endif
        call prod  (isiorb,psi(iborb),wk0)
        call add   (isiorb,wk0,wk1)
     endif

     call copy (isiorb,wk1,ione,wk2,ione)
     call prod  (isiorb,psi(iborb),wk2)
     w=dot(isiorb,wgt1,ione,wk2,ione)
     !       write(*,'(3e22.14)') w,wk1(isiorb/2),wk2(isiorb/2)
     !       write(*,'(3e22.14)') psi(iborb),psi(iborb+3*isiorb/4)
     !       write(*,*) isiorb,iborb
     vk(iorb)=w

     call copy (isiorb,f0,ione,wk0,ione)
     call prod  (isiorb,psi(iborb),wk0)

     call prod2 (isiorb,psi(iborb),wk0,wk2)
     w=dot(isiorb,wgt1,ione,wk2,ione)
     vn(iorb)=w

     call add (isiorb,wk0,wk1)
     call prod (isiorb,psi(iborb),wk1)

     w = dot(isiorb,wgt1,ione,wk1,ione)
     woneel=woneel+oc*w

     if (nel.lt.2) then
        vkt=vkt+oc*vk(iorb)
        vnt=vnt+oc*vn(iorb)
        evt=woneel
        epott=vnt
        virrat=(epott+z1*z2/r)/vkt
        etot=evt+z1*z2/r
     endif

     !      contribution from coulomb interaction within the same shell

     if (oc.gt.one) then
        if (isipot.ne.isiorb) then
           write(iout6,*) 'coulomb potentials and orbitals have to defined on the same number of subgrids'
           stop 'contrib'
        endif
        call copy (isipot,pot(ibpot),ione,wk2,ione)
        call prod  (isiorb,psi(iborb),wk2)
        call prod  (isiorb,psi(iborb),wk2)

        w=dot(isiorb,wgt2,ione,wk2,ione)
        wdcoul=wdcoul+oc*(oc-one)/2.0_PREC*w
     endif

     vkt=vkt+oc*vk(iorb)
     vnt=vnt+oc*vn(iorb)
     write(*,6000) iorn(iorb),bond(iorb),gut(iorb),occ(iorb),vk(iorb),vn(iorb),oc*(vk(iorb)+vn(iorb))
     write(*,6002) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb),occ(iorb),w,oc*(oc-one)/2.0_PREC*w
  enddo
6000 format(' ','Vk,Vn',i4,1x,a8,a1,f6.1,2x,5e22.14)
6002 format(' ','V(c) ',i4,1x,a8,a1,i4,1x,a8,a1,f6.1,5e22.14)

  !   contribution from coulomb and exchange interaction between shells

  wndc =zero
  wex1 =zero
  wex2 =zero

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

     if (ipe1.gt.0.and.abs(oc1-one).gt.eps) then
        call exint (iorb1,ocx1)
        ibex=i3b(iex1)
        isiex1=i3si(iex1)

        ngrid=min(isiorb1,isiex1)
        call prod2  (ngrid,wk0,excp(ibex),wk1)
        w=dot (ngrid,wgt2,ione,wk1,ione)
        wex1=wex1+ocx1*w
        write(*,7020) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb1),bond(iorb1),gut(iorb1),ocx1,w,wex1
7020    format(' ','v(x) ',i4,1x,a8,a1,i4,1x,a8,a1,f6.1,2e22.14)
     endif

     do iorb2=iorb1+1,norb
        iborb2=i1b(iorb2)
        isiorb2=i1si(iorb2)
        nmut2=i1mu(iorb2)
        ibpot2=i2b(iorb2)
        isipot2=i2si(iorb2)
        oc2=occ(iorb2)

        !         coulomb interaction between shells

        ngrid=min(isipot2,isiorb1)
        call prod2 (ngrid,pot(ibpot2),wk0,wk1)

        w=dot(ngrid,wgt2,ione,wk1,ione)
        wndc=wndc+oc1*oc2*w

        write(*,7022) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb2),bond(iorb2),gut(iorb2),oc1*oc2,w,oc1*oc2*w,wndc
7022    format(' ','v(c) ',i4,1x,a8,a1,i4,1x,a8,a1,f6.1,4e22.14)

        ipe2=mgx(6,iorb2)
        iex=iorb1+iorb2*(iorb2-1)/2
        ibex=i3b(iex)
        isiex=i3si(iex)

        !         exchange interaction between shells

        call excont (iorb1,iorb2,ocx1,ocx2)

        ngrid=min(isiorb1,isiorb2,isiex)
        call prod2 (ngrid,excp(ibex),psi(iborb1),wk1)
        call prod  (ngrid,psi(iborb2),wk1)

        w=dot(ngrid,wgt2,ione,wk1,ione)
        wex1=wex1+ocx1*w

        write(*,7024) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb2),bond(iorb2),gut(iorb2),ocx1,w,ocx1*w,wex1
07024   format(' ','v(x)a',i4,1x,a8,a1,i4,1x,a8,a1,f6.1,4e22.14)

        !         exchange interaction between shells (same non-zero lambdas )

        if (ilc(iex).gt.1) then
           call prod2 (ngrid,excp(ibex+ngrid),psi(iborb1),wk1)
           call prod  (ngrid,psi(iborb2),wk1)
           w=dot (ngrid,wgt2,ione,wk1,ione)
           wex2=wex2+ocx2*w

           write(*,7026) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb2),bond(iorb2),gut(iorb2),ocx2,w,ocx2*w,wex2
7026       format(' ','v(x)b',i4,1x,a8,a1,i4,1x,a8,a1,f6.1,4e22.14)
        endif
     enddo
  enddo

  evt=woneel+wdcoul+wndc-wex1-wex2
  epott=vnt+wdcoul+wndc-wex1-wex2
  virrat=(epott+z1*z2/r)/vkt
  etot=evt+z1*z2/r

  engt(1)=evt

  write(*,8000) woneel,wdcoul+wndc,wex1+wex2,evt
8000 format(' ','E1, E2C, E2x, ET ',4e22.14)
  write(*,*)
end subroutine contrib
