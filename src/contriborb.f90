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
! c ### contriborb ###
! c
! c     Calculates various contributions to total energy due to a given
! c     orbital

module contriborb_m
  implicit none
contains
  subroutine contriborb(iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
    use params
    use discret
    use scf
    use commons8
    use util

    use blas_m
    use diffmu_m
    use diffnu_m
    use putin_m
    use putout_m

    implicit none
    integer :: i,i1beg,i1beg1,i2beg,i2beg1,i3beg,ihc,iorb,iorb1,ipc,isym,kex,ngex,&
         ngorb,ngorb1,ngpot,ngpot1,ngrid,ngrid2,nmut
    real (PREC) :: oc,w,wcoul,wcoul0,wexc0,wkin,woneel,woneel0,wexc1,wexc2,wpot,wtwoel
    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3
    wtwoel=zero

    i1beg=i1b(iorb)
    i2beg=i2b(iorb)
    nmut=i1mu(iorb)
    ngorb=i1si(iorb)
    ngpot=i2si(iorb)

    if (ngorb.ne.ngpot) then
       write(*,*) 'Orbitals and corresponding Coulomb potentials have to be defined on the same number of grid points'
       stop 'contriborb'
    endif

    ipc=iorb+(iorb-1)*norb
    engo(ipc)=zero
    isym=isymOrb(iorb)

    !   calculate derivatives over mu and ni

    call putin (nni,nmut,isym,psi(i1beg),wk3)
    call diffnu (nmut,wk3,wk0,wk1,wk2)
    call putout (nni,nmut,wk1,wk0)

    call diffmu (nmut,wk3,wk2)
    call putout (nni,nmut,wk0,wk2)

    !   add contribution from derivatives over mu and ni

    call add (ngorb,wk0,wk1)

    call copy (ngorb,wk1,ione,wk2,ione)
    call prod  (ngorb,psi(i1beg),wk2)
    wkin=dot(ngorb,wgt1,ione,wk2,ione)

    call copy (ngorb,f0,ione,wk0,ione)
    if (mm(iorb).ne.0) then
       !      nuclear energy for non-sigma orbitals contains contribution
       !      from e term (in toten this term is correctly added to the
       !      kinetic energy contribution); e enters the expression with
       !      minus sign which is already incorporated in e

       w=dble(mm(iorb)*mm(iorb))
       call axpy (ngorb,w,e,ione,wk0,ione)
    endif

    call prod2 (ngorb,psi(i1beg),wk0,wk2)

    call copy (ngorb,wk2,ione,wk3,ione)

    call prod  (ngorb,psi(i1beg),wk2)
    wpot=dot(ngorb,wgt1,ione,wk2,ione)

    call proda (ngorb,psi(i1beg),wk0,wk1)
    call prod (ngorb,psi(i1beg),wk1)

    woneel=dot(ngorb,wgt1,ione,wk1,ione)
    write(*,8020) iorn(iorb), bond(iorb), gut(iorb),wkin,wpot,woneel
08020 format(/' ','Ek  En  Ek+En      ',i4,1x,a8,a1,4e22.14)
    woneel0=woneel

    ! FIXME
    do iorb1=1,norb
       call copy (ngorb,wk3,ione,wk1,ione)
       i1beg1=i1b(iorb1)
       call prod (ngorb,psi(i1beg1),wk1)
       w=dot(ngorb,wgt1,ione,wk1,ione)
       write(*,8022) iorn(iorb), bond(iorb), gut(iorb),iorn(iorb1),bond(iorb1),gut(iorb1),w
08022  format(' ','En   ',i4,1x,a8,a1,i4,1x,a8,a1,4e22.14)
    enddo

    if (nel.eq.1) return

    wcoul0=zero
    wexc0 =zero

    do i=1,mxsize
       wk2(i)=zero
    enddo

    if (norb.eq.1) then

       !      contribution from coulomb potential of one sigma orbital

       !      14/12/02
       !      call copy (ngpot,pot(i2beg),ione,wk2,ione)
       oc=occ(norb)-one
       call axpy (ngpot,oc,pot(i2beg),ione,wk2,ione)
       call prod  (ngorb,psi(i1beg),wk2)

    else

       !      add contributions from coulomb and exchange potentials

       do iorb1=1,norb

          !         it is asumed that exchange potentials involving iorb1
          !         has already been retrieved from disk

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
             wk3(i)=zero
          enddo

          call copy (ngpot1,pot(i2beg1),ione,wk3,ione)
          call prod  (ngrid,psi(i1beg),wk3)

          if (iorb.eq.iorb1) oc=oc-1.0_PREC
          if (abs(oc-1.0_PREC) .lt. precis) then
             call scal (ngrid,oc,wk3,ione)
          endif

          call prod (ngorb,psi(i1beg),wk3)
          wcoul=dot(ngorb,wgt2,ione,wk3,ione)

          woneel=woneel+wcoul
          wcoul0=wcoul0+wcoul
          write(*,8024) iorn(iorb), bond(iorb), gut(iorb),iorn(iorb1),bond(iorb1),gut(iorb1),oc,wcoul/occ(iorb1),woneel
08024     format(' ','v(c) ',i4,1x,a8,a1,i4,1x,a8,a1,f5.1,4e22.14)

          if (iorb1.ne.iorb)  then

             do i=1,mxsize
                wk3(i)=zero
             enddo

             ngrid2=min(ngorb1,ngex,ngorb)
             oc=gec(kex)
             call prodas (ngrid2,-oc,psi(i1beg1),excp(i3beg),wk3)

             call prod (ngorb,psi(i1beg),wk3)
             wexc1=dot(ngorb,wgt2,ione,wk3,ione)

             woneel=woneel+wexc1
             wexc0=wexc0+wexc1

             write(*,8026) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb1),bond(iorb1),gut(iorb1),oc,wexc1/oc,woneel
08026        format(' ','v(x)1',i4,1x,a8,a1,i4,1x,a8,a1,f5.1,4e22.14)


             if (ilc(ihc).gt.1) then
                do i=1,mxsize
                   wk3(i)=0.0_PREC
                enddo

                oc=gec(kex+norb*norb)
                call prodas (ngrid2,-oc,psi(i1beg1),excp(i3beg+ngex),wk3)

                call prod (ngorb,psi(i1beg),wk3)
                wexc1=dot(ngorb,wgt2,ione,wk3,ione)

                woneel=woneel+wexc1
                wexc0=wexc0+wexc1

                write(*,8028) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb1),bond(iorb1),gut(iorb1),oc,wexc1/oc,woneel
08028           format(' ','v(x)2',i4,1x,a8,a1,i4,1x,a8,a1,f5.1,4e22.14)
             endif
          else
             do i=1,mxsize
                wk3(i)=zero
             enddo

             if ((mm(iorb).gt.0).and.(ilc(ihc).gt.0)) then
                ngrid2=min(ngorb1,ngex,ngorb)
                oc=gec(kex)
                call prodas(ngrid2,-oc,psi(i1beg1),excp(i3beg),wk3)
                call prod (ngorb,psi(i1beg),wk3)
                wexc2=dot(ngorb,wgt2,ione,wk3,ione)
                woneel=woneel+wexc2
                wexc0=wexc0+wexc2
                write(*,8030) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb1),bond(iorb1),gut(iorb1),oc,wexc2/oc,woneel
08030           format(' ','v(x)3',i4,1x,a8,a1,i4,1x,a8,a1,f6.1,4e22.14)
             endif
          endif
       enddo
    endif

    write(*,8032) iorn(iorb),bond(iorb),gut(iorb),woneel0,wcoul0,wexc0
    write(*,*)
08032 format(' ','Ek+En v(c) v(x)123 ',i4,1x,a8,a1,4e22.14)

  end subroutine contriborb
end module contriborb_m

