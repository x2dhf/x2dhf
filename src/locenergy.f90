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
! ### locenergy ###
!
!     Calculates the locale energy for a given orbital.

subroutine locenergy (iorb,psi,pot,excp,e,f0,f4,wgt1,wgt2,wk0,wk1,wk2,wk3)
  use params
  use discret
  use scf
  use commons8

  use blas_m
  use diffmu_m
  use diffnu_m

  implicit none
  integer :: i,i1beg1,i1beg,i2beg,i2beg1,i3beg,ibeg,ihc,im,imax,in,iorb,iorb1,ioutmat, &
       ipc,ipe,isym,kex,ngex,ngorb,ngorb1,ngpot,ngpot1,ngrid,ngrid2,nmut
  real (PREC) :: oc,w,w1,w2,wk2max,wtwoel
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,f4

  character*13 :: fn


  ioutmat=30

  open(ioutmat,file='vni',status='replace',form='formatted')
  write(ioutmat,*) nni
  write(ioutmat,1000) (vni(in),in=1,nni)
  close(ioutmat)

  open(ioutmat,file='vmu',status='replace',form='formatted')
  write(ioutmat,*) i1mu(1)
  write(ioutmat,1000) (vmu(im),im=1,i1mu(1))
  close(ioutmat)
01000 format(F25.15)

  wtwoel=zero

  i1beg=i1b(iorb)
  i2beg=i2b(iorb)
  nmut=i1mu(iorb)
  ngorb=i1si(iorb)
  ngpot=i2si(iorb)

  ipc=iorb+(iorb-1)*norb
  engo(ipc)=zero
  isym=isymOrb(iorb)

  !    calculate derivatives over mu and ni

  call putin (nni,nmut,isym,psi(i1beg),wk3)
  call diffnu (nmut,wk3,wk0,wk1,wk2)
  call putout (nni,nmut,wk1,wk0)

  call diffmu (nmut,wk3,wk2)
  call putout (nni,nmut,wk0,wk2)

  !    add contribution from derivatives over mu and ni

  call add (ngorb,wk0,wk1)

  if (ipe.eq.0) then
     call copy (ngorb,f0,ione,wk0,ione)
  else

     !        e enter the expression with minus sign which is already
     !        incorporated in e

     w=dble(ipe*ipe)
     call copy (ngorb,f0,ione,wk0,ione)
     call axpy (ngorb,w,e,ione,wk0,ione)
  endif

  call proda (ngorb,psi(i1beg),wk0,wk1)
  call prod (ngorb,wgt1,wk1)

  !  now wk1 contains (T+V)|psi>

  !	  call prod  (ngorb,psi(i1beg),wk1)
  !	  w=dot(ngorb,wgt1,ione,wk1,ione)
  !	  write(*,*) '<T+V> ',w
  !        stop

  do i=1,ngorb
     wk2(i)=zero
  enddo

  if (nel.gt.1) then
     if (norb.eq.1) then

        !           contribution from coulomb potential of one sigma orbital

        call copy (ngpot,pot(i2beg),ione,wk2,ione)
        call prod  (ngorb,psi(i1beg),wk2)
     else
        do i=1,mxsize
           wk2(i)=zero
        enddo

        !    add contributions from coulomb and exchange potentials

        do iorb1=1,norb
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
              if ((ipe.gt.0).and.(ilc(ihc).gt.0)) then
                 ngrid2=min(ngorb1,ngex,ngorb)
                 oc=gec(kex)
                 call prodas (ngrid2,-oc,psi(i1beg1),excp(i3beg),wk0)
              endif
           endif

           if (iorb1.eq.1) then
              call copy (ngrid,wk0,ione,wk2,ione)
           else
              call add   (ngrid,wk0,wk2)
           endif
        enddo
     endif

     !         to complete the integrand wk2 has to be multiplied by psi(i1beg)

     call prod (ngorb,wgt2,wk2)
  endif

  !     wk2 contains V(fock)|psi>

  i1beg=i1b(iorb)
  ngorb=i1si(iorb)

  call add   (ngorb,wk1,wk2)

  call copy (ngorb,wk2,ione,wk3,ione)

  !      call prod  (ngorb,psi(i1beg),wk2)
  !      w=dot(ngorb,wgt1,ione,wk2,ione)
  !      write(*,*) '<T+V> ',w

  w=-eng(iorb)

  call copy(ngorb,f4,ione,wk1,ione)
  call prod(ngorb,psi(i1beg),wk1)
  call prod(ngorb,wgt2,wk1)
  call copy (ngorb,wk1,ione,wk0,ione)
  call axpy(ngorb,w,wk1,ione,wk2,ione)

  w=zero
  w1=zero
  w2=zero
  wk2max=zero
  do i=1,ngorb

     !         exclude grid points where nuclei are located

     if (abs(wk2(i)).gt.wk2max) then
        wk2max=abs(wk2(i))
        imax=i
     endif
     if (.not.(i.eq.1.or.i.eq.(ngrid-mxnmu))) then
        w=w+wk2(i)*wk2(i)
        w1=w1+wk2(i)*wk2(i)*f4(i)*wgt2(i)
     endif
  enddo

  write(*,*) 'imax, wk2max',imax,wk2max

  ! FIXME
  if (idbg(496).ne.0) then
     go to ( 11, 12, 13, 14, 15, 16), iorb
00011 fn='bf-1p.lenergy'
     goto 100
00012 fn='bf-5s.lenergy'
     goto 100
00013 fn='bf-4s.lenergy'
     goto 100
00014 fn='bf-3s.lenergy'
     goto 100
00015 fn='bf-2s.lenergy'
     goto 100
00016 fn='bf-1s.lenergy'
00100 continue
     open(ioutmat,file=fn,status='replace',form='formatted')
     ibeg=i1b(iorb)
     call prtmat (nni,i1mu(1),wk2,ioutmat)
     close(ioutmat)
  endif

  !       write(*,1010) iorn(iorb),bond(iorb),w,sqrt(w),w1,w2,sqrt(w2)
  !       write(*,1010) iorn(iorb),bond(iorb),w,sqrt(w),w1,sqrt(w1)
  write(*,1010) iorn(iorb),bond(iorb),w1,sqrt(w1),sqrt(w1)/(-eng(iorb))
01010 format(10x,i2,1x,a5,'  ',4e15.5)

  !        call enrradius(wk2)

  ! FIXME
  if (idbg(497).ne.0) then
     call leexact1(wk3,wk0)
     !          call divide(wk0,psi(i1beg))
  else
     !         call leexact(wk1,wk2)
  endif

end subroutine locenergy


