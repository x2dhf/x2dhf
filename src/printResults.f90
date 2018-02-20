! *x**************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                        *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### printResults ####
!
!     Prints total energy and its components, orbital energies, multipole
!     moments, etc upon completion of the SCF process

subroutine printResults (cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,ibeg,iorb,inuclder,jorb,k,ngrid

  real (PREC) :: etot_final,tcpu,xnorm,zarea
  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch
  real (PREC), dimension(10) :: totalref
!  real (PREC), external :: total

  !     calculate final total energy
  !     check orhogonality of orbitals
  if (iprint(30).ne.0) then
     print *,''
     print *,'orthogonality of orbitals:'
     iprint(20)=1
     do iorb=1,norb
        jorb=norb+1-iorb
        call checkOrtho (jorb,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
     enddo
     iprint(20)=0
  endif


  write (*,*)
  if (islat.eq.0) then
     call etotal (cw_orb,cw_coul,cw_exch,                          &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),   &
          cw_suppl(i4b(14)),                                       &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),   &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),   &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),   &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),   &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))
  else
     call etotalDFT (cw_orb,cw_coul,cw_exch,                       &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),   &
          cw_suppl(i4b(14)),                                       &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),   &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),   &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),   &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),   &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))
  endif

  if (ifermi.eq.2.and.iprint(590).ne.0) call etotalGauss
  call prttoten

  etot_final=etot

  if (iprint16.eq.0) then
     write(*,*)
     write(*,'(" total energy contributions: ")')
     write(*,'("     nuclear-electron attraction energy = ",f20.12)') ennucel
     write(*,'("     kinetic energy                     = ",f20.12)') enkin
     write(*,'("     one electron energy                = ",f20.12)') enkin+ennucel

     if (islat.eq.0) then
        ! FIXME print DFT functional used
        write(*,'("     coulomb energy                     = ",f20.12)') encoul
        write(*,'("     exchange energy                    = ",f20.12)') enexch
     end if

     ! FIXME print DFT functional used
     write(*,'("     nuclear repulsion energy           = ",f20.12)') z1*z2/r
     write(*,'("     coulomb energy (DFT)               = ",f20.12)') encouldft
     write(*,'("     exchange energy (DFT)              = ",f20.12)') enexchdft
     write(*,'("     correlation energy (DFT)           = ",f20.12)') edftcorr
     endif
  else

     write(*,*)
     write(*,'(" total energy contributions: ")')
     write(*,'("     nuclear-electron attraction energy = ",f28.22)') ennucel
     write(*,'("     kinetic energy                     = ",f28.22)') enkin
     write(*,'("     one electron energy                = ",f28.22)') enkin+ennucel
     if (islat.eq.0) then
        ! FIXME print DFT functional used
        write(*,'("     coulomb energy                     = ",f28.22)') encoul
        write(*,'("     exchange energy                    = ",f28.22)') enexch
     end if

     ! FIXME print DFT functional used
     write(*,'("     nuclear repulsion energy           = ",f28.22)') z1*z2/r
     write(*,'("     coulomb energy (DFT)               = ",f28.22)') encouldft
     write(*,'("     exchange energy (DFT)              = ",f28.22)') enexchdft
     write(*,'("     correlation energy (DFT)           = ",f28.22)') edftcorr
  endif


  if (idft.eq.1) then
     islat=1
     imethod=3
     idftex=1
     call etotalDFT (cw_orb,cw_coul,cw_exch,                      &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),  &
          cw_suppl(i4b(14)),                                      &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),  &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),  &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))
     write(*,*)
     if (iprint16.eq.0) then
        write(*,'("     exchange energy: LDA               = ",f20.12)') enexchdft
     else
        write(*,'("     exchange energy: LDA               = ",f28.22)') enexchdft
     endif

     idftex=2
     call etotalDFT (cw_orb,cw_coul,cw_exch,                      &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),  &
          cw_suppl(i4b(14)),                                      &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),  &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),  &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))

     if (iprint16.eq.0) then
        write(*,'("     exchange energy: B88               = ",f20.12)') enexchdft
     else
        write(*,'("     exchange energy: B88               = ",f28.22)') enexchdft
     endif

     idftex=3
     call etotalDFT (cw_orb,cw_coul,cw_exch,                      &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),  &
          cw_suppl(i4b(14)),                                      &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),  &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),  &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))

     if (iprint16.eq.0) then
        write(*,'("     exchange energy: PW86              = ",f20.12)') enexchdft
     else
        write(*,'("     exchange energy: PW86              = ",f28.22)') enexchdft
     endif

     idftex=4
     call etotalDFT (cw_orb,cw_coul,cw_exch,                      &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),  &
          cw_suppl(i4b(14)),                                      &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),  &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),  &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))

     if (iprint16.eq.0) then
        write(*,'("     exchange energy: PW91              = ",f20.12)') enexchdft
     else
        write(*,'("     exchange energy: PW91              = ",f28.22)') enexchdft
     endif

     idftex=0
     idftcorr=1
     call etotalDFT (cw_orb,cw_coul,cw_exch,                         &
             cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),  &
             cw_suppl(i4b(14)),                                      &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
             cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),  &
             cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),  &
             cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
             cw_sctch(i5b(13)),cw_sctch(i5b(14)))

     if (iprint16.eq.0) then
        write(*,'("     correlation energy: LYP            = ",f20.12)') edftcorr
     else
        write(*,'("     correlation energy: LYP            = ",f28.22)') edftcorr
     endif

     idftcorr=2
     call etotalDFT (cw_orb,cw_coul,cw_exch,                      &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),  &
          cw_suppl(i4b(14)),                                      &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),  &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),  &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)),  &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))

     if (iprint16.eq.0) then
        write(*,'("     correlation energy: VWN            = ",f20.12)') edftcorr
     else
        write(*,'("     correlation energy: VWN            = ",f28.22)') edftcorr
     endif
  endif

  !     write final orbital energies and normalization factors

  write(*,*)
  !      iprint16=1
  if (iprint16.eq.0) then
     write(iout6,1000)
     do i=1,norb
        write(iout6,1002) iorn(i),bond(i),gut(i),eng(i),1.0-area(i)
     enddo
  else
     write(iout6,7000)
     do i=1,norb
        write(iout6,7002) iorn(i),bond(i),gut(i),eng(i),1.0-area(i)
     enddo
  endif
  write(*,*)

  !     check how an error in the orbital normalization factors influence
  !     the accuracy of the total energy

  !     denormalize orbitals (they leave SCF process properly normalized)
  do iorb=1,norb
     ibeg = i1b (iorb)
     ngrid= i1si(iorb)
     xnorm=area(iorb)
     zarea = sqrt(xnorm)
     call scal (ngrid,zarea,cw_orb(ibeg),ione)
  enddo

  !     calculate final total energy
  write (*,*)
  etot_final=etot

  if (islat.eq.0) then
     call etotal (cw_orb,cw_coul,cw_exch,                        &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
          cw_suppl(i4b(14)),                                     &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)), &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)), &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)), &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)), &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))
  else
     call etotalDFT (cw_orb,cw_coul,cw_exch,                     &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)), &
          cw_suppl(i4b(14)),                                     &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)), &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)), &
          cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)), &
          cw_sctch(i5b(10)),cw_sctch(i5b(11)),cw_sctch(i5b(12)), &
          cw_sctch(i5b(13)),cw_sctch(i5b(14)))
  endif

  write(iout6,6111) abs(etot_final-etot),abs(etot_final/etot-1.0_PREC)*100

  !     normalize orbitals again to bring them to the previous state
  do iorb=1,norb
     ibeg = i1b (iorb)
     ngrid= i1si(iorb)
     xnorm=area(iorb)
     zarea = 1.0_PREC/sqrt(xnorm)
     call scal (ngrid,zarea,cw_orb(ibeg),ione)
  enddo

  !     check orthohonality
  if (iprint(191).ne.0) then
     write(*,*) 'checking orthogonality:'
     do iorb=1,norb
        i=norb+1-iorb
        call checkOrtho (i,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b(1)))
     enddo
  endif

  !     check symmetry of orbitals in homonuclear case

  if (abs(z1-z2).lt.homolevl.and.ibreak.eq.0) call checkSym(cw_orb)

!      call orbtails (cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch)

  if (iprint(192).ne.0) then
     write(*,*)
     write(*,*) 'Euclidean norms of (T+V(n)+V-E)|i>:'
     do iorb=1,norb
        call locenergy (iorb,cw_orb,cw_coul,cw_exch,                 &
             cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b( 9)),  &
             cw_suppl(i4b(13)),cw_suppl(i4b(14)),                    &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
             cw_sctch(i5b( 4)))
     enddo
  endif

  !     check multipole expansion contributions to potentials
  if (iprint(193).ne.0) then
     call checkPot (cw_coul,cw_exch)
  endif

  !     calculate multipole moments and expectation values


  if (imethod.ne.2) then
     call propet2 (cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 9)),      &
          cw_suppl(i4b(14)),                                      &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
          cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)))


     if (iprint(152).ne.0) then
        !           save total multipole moments
        do k=1,maxmpole
           totalref(k)=total(k)
        enddo

        !     denormalize orbitals (they leave SCF process properly normalized)
        do iorb=1,norb
           ibeg = i1b (iorb)
           ngrid= i1si(iorb)
           xnorm=area(iorb)
           zarea = sqrt(xnorm)
           call scal (ngrid,zarea,cw_orb(ibeg),ione)
        enddo

        !           calculate total multipole moments

        call propet3 (cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 9)),      &
             cw_suppl(i4b(14)),                                      &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
             cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)))

        !           normalize orbitals again to bring them to the previous state
        do iorb=1,norb
           ibeg = i1b (iorb)
           ngrid= i1si(iorb)
           xnorm=area(iorb)
           zarea = 1.0_PREC/sqrt(xnorm)
           call scal (ngrid,zarea,cw_orb(ibeg),ione)
        enddo


        !           print relative errors of multipole moments

        write(*,2000)
        do k=1,4
           if (abs(totalref(k)).gt.1d-12) then
              write(*,2002) k,abs(one - total(k)/totalref(k))
           else
              write(*,2002) k,zero
           endif
        enddo

        !           restore total multipole moments
        do k=1,maxmpole
           total(k)=totalref(k)
        enddo
     endif
  endif

  if (iprint(140).ne.0) then
     call propet4 (cw_orb,cw_coul,cw_exch,cw_suppl(i4b( 9)),      &
             cw_suppl(i4b(14)),                                      &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),  &
             cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)))
  endif



2000 format(//5x,"relative errors of moments due to orbital norms not being equal 1:")
2002 format(10x,"k=",i1,1P,e12.2)

  !    calculate the derivative of the orbital density with respect to z
  !    at A

  if (iprint(150).ne.0) then
     inuclder=0
     if (inuclder.eq.0) then
        call nuclder(cw_orb,                                                           &
             cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),cw_suppl(i4b(14)),  &
             cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),  &
             cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)))
     endif
  endif

  ! FIXME
  !     calling scmc is not needed at the end of calculations
  if (iscmc.eq.1) then
     call scmc (cw_orb,cw_coul,cw_exch,                                                          &
          cw_suppl(i4b( 4)),cw_suppl(i4b( 5)),cw_suppl(i4b(13)),cw_suppl(i4b(14)),                   &
          cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),cw_sctch(i5b( 5)), &
          cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),cw_sctch(i5b(10)), &
          cw_sctch(i5b(11)),cw_sctch(i5b(12)),cw_sctch(i5b(13)),cw_sctch(i5b(14)))
  endif

  !     Calculates the total radial densities relative to centres A
  !     (z-R/2) and B (z+R/2)

  if (iprint(110).ne.0.or.iprint(111).ne.0) then
     call radialden(cw_orb,cw_sctch(i5b( 1)))
  endif

  call separator

  !     CPU time statistics

  tcpu=trelax+tortho+trayl+tmomen

  write(*,'(" cpu summary (sec):")')
  !      write(*,5003) trayl
  !      write(*,5007) tlagra
  write(*,5009) tlagra+trayl
  write(*,5002) tortho
  write(*,5004) tmomen
  write(*,5001) trelax
  write(*,5008) ttoten
  write(*,5005) tcpu

  if (iscf.ne.0) write(*,5006) tcpu/dble(iscf)

1000 format(7x,'orbital',15x,'energy',21x,'1-norm')
1002 format(1x,i3,1x,a8,1x,a1,e28.16,e20.2)
5001 format(4x,'relaxation of orb. and pot.        ', f12.2)
5002 format(4x,'normalization+orthogonalization    ', f12.2)
!5003 format(4x,'diagonal Lagrange multipliers      ', f12.2)
!5007 format(4x,'non-diagonal Lagrange multipliers  ', f12.2)
5009 format(4x,'Lagrange multipliers               ', f12.2)
5004 format(4x,'multipole moments                  ', f12.2)
5005 format(4x,'all SCF iterations                 ', f12.2)
5006 format(4x,'single SCF iteration               ', f12.2/)
5008 format(4x,'total energy                       ', f12.2)
!6100 format(1x,'total electronic energy: ',e28.16)
!6110 format(1x,'total energy:            ',e28.16)
6111 format(1x,'total energy uncertainty due to orbital norms not being equal 1:'/,5x, &
  & 'absolute = +/-',e8.2,',  relative = +/-',e8.2,'%')
!6120 format(1x,'virial ratio:            ',e28.16)
7000 format(7x,'orbital',22x,'energy',30x,'1-norm')
7002 format(1x,i3,1x,a8,1x,a1,e44.32,2e20.2)
!7100 format(1x,'total electronic energy: ',e44.32)
!7110 format(1x,'total energy:            ',e44.32)
!7111 format(1x,'total energy uncertaininty due to orbital norms not being equal 1:'/, &
!          & '  absolute   +/-',e8.2,/,'  relative   +/-',e8.2,'%')
!7120 format(1x,'virial ratio:            ',e44.32)

end subroutine printResults










