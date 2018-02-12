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
! ### doSCF ####

!     Controls SCF process.

subroutine doSCF (cw_sor,cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,iasympt,ibeg,ic,iend,imomen,inde,indn,iorb,istat,modv,nei,next,noenergydec,nonormdec,nthren
  real (PREC) :: ddmax,denmax,ddmaxprev,dnmaxprev,dnomax,thren,thrno,time1,time2,tscf

  integer, dimension(*) :: cw_sor

  real (PREC), dimension(maxorb) :: demaxt
  real (PREC), dimension(*) ::  cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch

  logical :: stop_x2dhf

  istat=1
  stop_x2dhf=.false.
  ic=0
  thren=10.0_PREC**dble(-ienterm)
  thrno=10.0_PREC**dble(-inoterm)

  !    the SCF procedure starts here: set some parameters and go

  ddmax=0.0_PREC
  do i=1,norb
     itouch(i)=0
     demax(i)=100.0_PREC
     denmax=100.0_PREC
     demaxt(i)=0.0_PREC
     engi(i)=eng(i)
  enddo
  dnomax=100.0_PREC

  tscf=0.0_PREC
  trelax=0.0_PREC
  tortho=0.0_PREC
  trayl =0.0_PREC
  tlagra=0.0_PREC
  tmomen=0.0_PREC
  ttoten=0.0_PREC

  !   SCF/SOR process is stopped if orbital energy differences are
  !   either too large or too small

  ddmaxprev=1.e35_PREC
  dnmaxprev=1.e35_PREC
  noenergydec=0
  nonormdec=0

  call getCpuTime (time1)
  tscf=time1

  iscf=0

  time1=0.0_PREC
  time2=0.0_PREC

  !   inforce C_i symmetry or orbitals
  if (ihomon.eq.2) then
     call setci(izero,cw_orb)
  endif

  !   calculate total energy in the case of exchange potentials being
  !   kept on disk during SCF cycles

  if (iform.eq.0.or.iform.eq.2) then
     iorb=0
     call etotalOrb (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
          cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),    &
          cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),    &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)), cw_sctch(i5b(13)),cw_sctch(i5b(14)))
  endif

  !   prepear an array needed for writing exchange potentials involving a given
  !   orbital

  call prepwexch

  !   ----- begin of scf loop -----

  ioo=1
  imomen=0

  if (maxscf.le.0) write(iout6,*) '... skipping scf iterations ...'

  nthren=0
  do iscf=1,maxscf

     if (iscf.eq.1) then
        if (iprint16.eq.0) then
           write(iout6,15000)
        else
           write(iout6,25000)
        endif
     endif

     !      ----- loop over orbitals (in reverse order) -----

     !      relaxation begins from the innermost orbital, i.e. the last one
     !      on the input data list unless iscforder specifies alternative
     !      ordering

     ibeg=1
     iend=norb

     ! FIXME
     !if (iorder(i).eq.5.and.mod(iscf,itwo).eq.0) then
     !   ioo=-1
     !else
     !   ioo= 1
     !endif
     ioo=1

     iasympt=norb

     do i=ibeg,iend
        if (iscforder(i).ne.0) then
           iorb=iscforder(i)
        else
           iorb=norb-i+1
        endif

        itouch(iorb)=1
        !          if (ifix(iorb).ne.0) goto 100

        call getCpuTime (time1)

        if (iform.eq.0.or.iform.eq.2) then
           !            retrieve from disk exchange potentials involving iorb
           call rfdexch(iorb,cw_exch)
        endif

        !         determines Coulomb and exchange potentials in the asymptotic region
        !         from the multipole expansion

        !          if (iasympt.ne.0) then
        if (iasympt.ne.0.and.ifix(iorb).eq.0) then
           call potAsympt(iorb,cw_coul,cw_exch)
           iasympt=iasympt-1
        endif

        !         in case of near-degenerate orbitals, i.e. when performing
        !         calculations for homonuclear molecules without inforced
        !         symmetry (no homo label), especially when external electric
        !         field is applied, one can improve convergence by changing
        !         the direction of the SOR sweeps (forward/backward sweeps for
        !         even/odd SCF)

        isstart= 1
        isstop = ngrd1
        isstep = 1

        if (ialtsweeps.eq.1) then
           if (mod(iorb,itwo).eq.0) then
              isstart= 1
              isstop = ngrd1
              isstep = 1
           else
              isstart= ngrd1
              isstop = 1
              isstep =-1
           endif
        endif

        !   perform relaxation of Coulomb and exchange potentials
        !   (maxsorpot(iorb) sweeps of each function) and orbitals
        !   (maxsororb(iorb) sweeps)

        call getCpuTime(time1)
        call doSOR (iorb,cw_sor,cw_orb,cw_coul,cw_exch,             &
             cw_suppl(i4b( 1)),cw_suppl(i4b( 2)),cw_suppl(i4b( 3)), &
             cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b( 6)), &
             cw_suppl(i4b( 7)),cw_suppl(i4b( 8)),cw_suppl(i4b( 9)), &
             cw_suppl(i4b(10)),cw_suppl(i4b(11)),cw_suppl(i4b(12)), &
             cw_suppl(i4b(14)),                                     &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)), &
             cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)), &
             cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)), &
             cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)), &
             cw_sctch(i5b(13)),cw_sctch(i5b(14)))
        call getCpuTime (time2)
        trelax=trelax+(time2-time1)

        !         normalize the current orbital
        call getCpuTime (time1)
        call norm (iorb,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
        call getCpuTime (time2)

        tortho=tortho+(time2-time1)

        !         enforce C_i symmetry of iorb orbital
        if (ihomon.eq.2) then
           call setci(iorb,cw_orb)
        endif

        !         orthogonalize the current orbital
        call getCpuTime (time1)
        if (norb.gt.1) then
           call ortho (iorb,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
        endif
        call getCpuTime(time2)
        tortho=tortho+(time2-time1)

        !         call rotate(iorb,cw_orb,cw_sctch(i5b(1)))

        !         calculate diagonal energy value

        call getCpuTime(time1)
        if (islat.eq.0) then
           call Ea (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
                cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
        else
           call EaDFT (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
                cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
        endif

        call getCpuTime(time2)
        trayl =trayl + (time2-time1)

        demaxt(iorb)=eng(iorb)-engi(iorb)

        !         sometimes demaxt becomes zero although it should not normally happen
        !         that is why it is increased (esspecially important for OED calculations
        !         when converged orbitals get 'fixed' status)

        if (abs(demaxt(iorb)).eq.00_PREC) demaxt(iorb)=two*thren

        engi(iorb)=eng(iorb)

        !         calculate off-diagonal Lagrange multipliers

        call getCpuTime(time1)
        if (islat.eq.0) then
           call Eab (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
                cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
        else
           call EabDFT (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
                cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
        endif
        call getCpuTime(time2)
        tlagra =tlagra + (time2-time1)

        if (iprint(50).ne.0) then
           if (islat.eq.0) then
              call contriborb (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
                   cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
           else
              call contriborbDFT (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
                   cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
           endif
        endif

        if (iform.eq.0.or.iform.eq.1) then
           call wtdexch(iorb,cw_exch)
        endif

        !        if (iprtlev.eq.1.and.ifix(iorb).eq.0) then
        if (iprtlev.eq.1) then
           if (iprint16.eq.0) then
              write(iout6,15100) iscf,iorn(iorb),bond(iorb),gut(iorb),eng(iorb),demaxt(iorb),1.0-area(iorb),wstorthog(iorb)
              !   &                 nodes(iorb,cw_orb)
           else
              write(iout6,25100) iscf,iorn(iorb),bond(iorb),gut(iorb),eng(iorb),demaxt(iorb),1.0-area(iorb),wstorthog(iorb)
              !   &                 nodes(iorb,cw_orb)
           endif
        endif

        if (nobckup.gt.0) then
           modv=mod(iscf,abs(nobckup))
           if (modv.eq.0) then
              if (iform.eq.0.or.iform.eq.2) then

                 !                  calculate 'intermidiate' total energy

                 call getCpuTime(time1)
                 call etotalOrb (iorb,cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),&
                      cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),   &
                      cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),   &
                      cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),cw_sctch(i5b(13)),cw_sctch(i5b(14)))
                 call getCpuTime(time2)
                 ttoten=ttoten+(time2-time1)
              endif
           endif
        endif

        !     Sometimes a lengthy calculation needs to be stoped earlier.
        !     kill <pid> will do the job but you risk that the data being
        !     written will be lost. To avoid this and force the program to end
        !     gracefully make
        !        touch stop_x2dhf
        !     or simply
        !        ./xhf stop

        !     If this file is detected the program ends.

        inquire (file ='stop_x2dhf',exist=stop_x2dhf)
        if (stop_x2dhf) goto 110

     enddo

     !   ----- end of orbital loop -----


     if (idbg(77).ne.0)  then
        !         assign maxsor for each orbital every iepoch iterations
        call schedSOR
     endif

00110 continue
     inde=1
     denmax=0.0_PREC

     indn=1
     dnomax=0.0_PREC

     if (norb.gt.0) then
        do nei=1,norb
!           if(ifix(nei).eq.0) then
              if (abs(demaxt(nei)).ge.abs(denmax)) then
                 denmax=demaxt(nei)
                 inde=nei
              endif
!           endif
        enddo

        do nei=1,norb
           if(ifix(nei).eq.0) then
              if (abs(1.0_PREC-area(nei)).ge.dnomax) then
                 dnomax=abs(1.0_PREC-area(nei))
                 if (dnomax.eq.00_PREC) dnomax=precis
                 indn=nei
              endif
           endif
        enddo
     endif

     !      inde and indn point to the worst converged orbital as far as
     !      orbital energies and norms are concerned, respectively

     if (ddmax.eq.0.0_PREC) ddmax=denmax

     if (iprtlev.eq.2) then
        if (iprint16.eq.0) then
           write(iout6,15110) iscf,iorn(inde),bond(inde),gut(inde),eng(inde),denmax,1.0-area(inde)
           write(iout6,15110) iscf,iorn(indn),bond(indn),gut(indn),eng(indn),demaxt(indn),1.0-area(indn)
        else
           write(iout6,25100) iscf,iorn(inde),bond(inde),gut(inde),eng(inde),denmax,1.0-area(inde)
           write(iout6,25100) iscf,iorn(indn),bond(indn),gut(indn),eng(indn),demaxt(indn),1.0-area(indn)
        endif

     endif


     modv=mod(iscf,abs(nobckup))
     if (modv.eq.0) then
        if (iprtlev.eq.3) then
           if (iprint16.eq.0) then
              write(iout6,15110) iscf,iorn(inde),bond(inde),gut(inde),eng(inde),denmax,1.0-area(inde)
              write(iout6,15110) iscf,iorn(indn),bond(indn),gut(indn),eng(indn),demaxt(indn),1.0-area(indn)
           else
              write(iout6,25100) iscf,iorn(inde),bond(inde),gut(inde),eng(inde),denmax,1.0-area(inde)
              write(iout6,25100) iscf,iorn(indn),bond(indn),gut(indn),eng(indn),demaxt(indn),1.0-area(indn)
           endif
        endif
     endif

     if (iscmc.eq.1) then
        call scmc (cw_orb,cw_coul,cw_exch,                                             &
             cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),cw_suppl(i4b(14)),  &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),  &
             cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),  &
             cw_sctch(i5b( 9)),cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
             cw_sctch(i5b(13)),cw_sctch(i5b(14)))
     endif

     !   calculate and print various contributions to total energy
     if (iprint(50).ne.0) then
        call contrib (cw_orb,cw_coul,cw_exch,                                         &
             cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),cw_suppl(i4b(14)), &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))
     endif


     !   Sometimes the convergence threshold (for energy and/or
     !   normalization) is set too low and cannot be satisfied on a given
     !   grid and as a result the SCF/SOR process continues in vain. To
     !   avoid these excessive iterations the convergence process is
     !   stopped if either denmax or dnomax has not decreased over a
     !   given number (nlast) most recent iterations.

     !   skip nscf2skip iterations before starting the search for
     !   saturation; examine the last nenlast and nnolast iterations in case
     !   of orbital energy and normalization, respectively

     if (iscf.gt.nscf2skip) then
        if (abs(denmax).ge.abs(ddmaxprev)) then
           noenergydec=noenergydec+1
        else
           ddmaxprev=denmax
           noenergydec=0
        endif

        if (abs(dnomax).ge.abs(dnmaxprev)) then
           nonormdec=nonormdec+1
        else
           dnmaxprev=dnomax
           nonormdec=0
        endif
        !         print *,noenergydec,nonormdec,ddmaxprev,dnmaxprev
        if (noenergydec.gt.nenlast.or.nonormdec.gt.nnolast) then
           write(*,*) '          '
           write(*,*) '... solution cannot be further improved ...'
           goto 9999
        endif
     else
        ddmaxprev=demax(1)
        dnmaxprev=dnomax
     endif

     !   check if the scf process is to be continued

     if (abs(demax(1)).gt.diver) then
        write(*,*) '           '
        write(*,*) '... maximal orbital energy error is too large ...'
        goto 9999
     endif

     if (iscf+1.gt.maxscf) then
        write(*,*) '          '
        write(*,*) '... scf iteration limit reached ... '
        goto 9999
     endif

     if (abs(denmax).lt.thren) then

        if (nthren.eq.0) then
           nthren=iscf
           !            recalculate multipole moments
           call mpoleMom(cw_orb,cw_suppl,cw_sctch)
           imomen=1

           ddmax=denmax
           iasympt=norb

           do iorb=1,norb
              itouch(iorb)=0
           enddo
           goto 8888
        else
           if (nthren+1.ne.iscf) then
              nthren=0
              goto 8888
           else
              write(*,*) '          '
              write(*,*) '... orbital energy threshold reached ...'
              goto 9999
           endif
        endif
     endif

     if ( thrno.ne.0.0_PREC ) then
        if (exlorb.eq.0.and.dnomax.lt.thrno) then
           write(*,*) '          '
           write(*,*) '... orbital normalization threshold reached ...'
           goto 9999
        endif
     endif

8888 continue

     !      in case of OED calculations skip relaxing orbitals for which
     !      either relaxation threshold has been achieved

     if (imethod.eq.2) then
        i=0
        do iorb=1,norb
           if (abs(demaxt(iorb)).lt.thren) ifix(iorb)=1
           if (ifix(iorb).eq.0) i=i+1
        enddo

        if (i.eq.0) then
           write(*,*) '          '
           write(*,*) '... orbital energy threshold reached reached ...'
           goto 9999
        endif
     endif


     if (stop_x2dhf) then
        write(*,*) '          '
        write(*,*) '... stop_x2dhf detected ... '
        write(*,*) '... program is terminating ...'
        write(*,*)
        next=0
        goto 9999
     endif

     !   recalculate multipole moments expansion coefficients

     if (abs(ddmax/denmax)==facmul.or. abs(denmax/ddmax)==facmul) then
        call mpoleMom(cw_orb,cw_suppl,cw_sctch)
        imomen=1

        ddmax=denmax
        iasympt=norb

        do iorb=1,norb
           itouch(iorb)=0
        enddo
     endif

     !      recalculate multipole coefficients and write functions to a disk
     !      file if nobckup > 0

     if (nobckup.gt.0) then
        modv=mod(iscf,abs(nobckup))
        if (modv.eq.0) then
           if (imomen.eq.0.and.facmul.gt.00_PREC) then
              call mpoleMom(cw_orb,cw_suppl,cw_sctch)
           endif
           imomen=0

           call writeDisk(cw_orb,cw_coul,cw_exch)

           !            calculate 'intermediate' total energy
           if (iform.eq.1.or.iform.eq.3) then
              if (islat.eq.0) then
                 call Etotal (cw_orb,cw_coul,cw_exch,                                           &
                      cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),cw_suppl(i4b(14)),  &
                      cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),  &
                      cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),  &
                      cw_sctch(i5b( 9)),cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
                      cw_sctch(i5b(13)),cw_sctch(i5b(14)))
              else
                 call EtotalDFT (cw_orb,cw_coul,cw_exch,                                        &
                      cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),cw_suppl(i4b(14)),  &
                      cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),  &
                      cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),  &
                      cw_sctch(i5b( 9)),cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
                      cw_sctch(i5b(13)),cw_sctch(i5b(14)))
              endif
              call getCpuTime(time2)
              ttoten=ttoten+(time2-time1)
           endif
           if (iprint16.eq.0) then
              write(iout6,6100)  etot
           else
              write(iout6,6101)  etot
           endif
        endif
     endif

     !      to run a job in a bach mode (or in the background) with the
     !      output redirected into a file and to be able to monitor every
     !      single iteration one has to flush the output buffer to disk

     call flush(iout6)
  enddo
  !   ----- end of scf loop	-----

09999 continue

  !   write functions to disk file
  if (nobckup.ge.0) then
     call writeDisk(cw_orb,cw_coul,cw_exch)
  endif

  if (iout4dd==1) call writeDisk4dd(cw_orb,cw_coul,cw_suppl(i4b(5)),cw_suppl(i4b(9)),&
       cw_sctch(i5b( 1)),cw_sctch(i5b(2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)))


  return

6100  format(1x,'total energy: ',e26.16)
6101  format(1x,'total energy: ',e38.28)
15000 format(/'   scf  orbital',11x,'   energy  ',12x,'energy diff ',7x,'1-norm',9x,'overlap')
15100 format(i5,i4,1x,a8,a1,2x,e22.16,3e16.2)
15110 format(i5,i4,1x,a8,a1,3x,e22.16,2e16.2,i5)
25000 format(/'  iter.   orbital',24x,'   energy  ',19x,'energy diff.',7x,'1-norm')
25100 format(i5,i4,1x,a8,a1,2x,e44.32,3e16.2)
end subroutine doSCF
