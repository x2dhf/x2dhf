! **************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### printCase ###
!
!     This routine can be used to check if input data have been
!     correctly transformed into values of numerous variable used by the
!     program.

subroutine printCase
  use params
  use discret
  use memory
  use scf

  use commons8

  implicit none

  integer :: ib,ie,ig,ior,iorb,ip,izz1,izz2,isizeint,isizereal
  integer :: lengtht,length0,maxorb2,mm3,omb

  real (PREC) ::  heta,hxibeg,hxiend,rinfig

  !     calculate the total length of working arrays (in bytes)
  !     system under consideration


  izz1=nint(z1)
  izz2=nint(z2)

  write(iout6,*)
  write(iout6,*) '  molecular system: '
  write(iout6,*)
  write(iout6,1000) element(izz1), z1, element(izz2),z2 ,r,r*0.5291772490_PREC

  if     (ipot.eq.1) then
     write(iout6,1401) z1atmass,element(izz1),z2atmass,element(izz2)
  elseif (ipot.eq.2) then
     write(iout6,1402) z1atmass,element(izz1),z2atmass,element(izz2)
  endif

  write(iout6,*)
  if (imethod.eq.1) write(iout6,*) '  method: HF'
  if (imethod.eq.2) then
     write(iout6,*) '  method: OED'
     write(iout6,*)
     if (ipot.eq.0) write(iout6,*) '  nuclear potential: Coulomb '
     if (ipot.eq.1) write(iout6,*) '  nuclear potential: Fermi '
     if (ipot.eq.2) write(iout6,*) '  nuclear potential: Gauss '
     if (ipot.eq.3) write(iout6,*) '  nuclear potential: Coulomb 3D '
     if (ipot.eq.4) write(iout6,*) '  nuclear potential: Kramers-Henneberger'
     if (ipot.eq.5) write(iout6,*) '  nuclear potential: Green-Sellin-Zachor'
  endif

  if (imethod.eq.3) write(iout6,1001) alphaf
  if (imethod.eq.4) then
     write(iout6,1002)
     if (idftex.eq.1) then
        write(iout6,1003) cdftex(idftex),alphaf
     else
        write(iout6,1004) cdftex(idftex)
     endif

     if (idftcorr.ne.0) then
        write(iout6,1005) cdftcorr(idftcorr)
     endif
  endif
  if (imethod.eq.5) then
     write(iout6,1006)
  endif



  if (ifefield.eq.1) write(iout6,1300) ffield
  if (iharm2xy.eq.1) write(iout6,1310) harm2xy

  !     its electronic configuration

  call printElConf


  if (iprint(175).eq.1) call fockform



  write(iout6,*)
  write(iout6,*) '  grid:'
  if (ngrids.eq.1) then
     write(iout6,1040) nni,hni,nmu(ngrids),hmu(ngrids),rinf
  else
     write(iout6,1042) nni,hni,mxnmu,rinf
     write(iout6,*)
     !         write(iout6,*) '         subgrid   nni   nmu    hni      hmu',
     !     &        '      rb     heta     hxib     hxie'
     write(iout6,"(10x,'subgrid',3x,'nni',3x,'nmu',4x,'hni',6x,'hmu',6x,'rb',5x,'heta',5x,'hxib',5x,'hxie')")

     rinf=r*vxi(iemu(ngrids))/2.0_PREC
     heta=abs(veta(nni)-veta(nni-1))
     ib=1
     do ig=1,ngrids
        rinfig=r*vxi(iemu(ig))/2.0_PREC
        ie=iemu(ig)
        hxibeg=vxi(ib+1)-vxi(ib)
        hxiend=vxi(ie)-vxi(ie-1)
        write(iout6,'(14x,i2,3x,i4,2x,i4,2f9.5,f7.3,3f9.5)') ig,nni,nmu(ig),hni,hmu(ig),rinfig,heta,hxibeg,hxiend
        ib=ie
     enddo
  endif

  write(iout6,*)
  write(iout6,*) '  scf: '
  write(iout6,1050) maxscf,ienterm,inoterm,facmul

  write(iout6,*)
  if       (exlorb.eq.0.0_PREC) then
     write (iout6,1051)
  else
     write (iout6,1052)
  endif
  if       (exlcoul.eq.0.0_PREC) then
     write (iout6,1053)
  else
     write (iout6,1054)
  endif
  if (imethod.eq.1) then
     if     (exlexp.eq.0.0_PREC) then
        write (iout6,1055)
     elseif (exlexp.eq.2.0_PREC) then
        write (iout6,1056)
     else
        write (iout6,1057)
     endif
  endif

  write(iout6,1058) mpole


  write(iout6,*)
  write(iout6,*) '  (mc)sor:'

  if     (ipoiss.eq.1) then
     write(iout6,*) '         sor method used for relaxing orbitals and potentials'
  elseif (ipoiss.eq.2) then
     write(iout6,*) '         mcsor method used for relaxing orbitals and potentials'
  elseif (ipoiss.eq.3) then
     write(iout6,*) '         sor method used for relaxing orbitals and mcsor for potentials'
  endif

  ior=iorder(ngrids)
  if (ialtsweeps.eq.0) then
     if (ior.eq.1) write(iout6,1071)
     if (ior.eq.2) write(iout6,1072)
     if (ior.eq.3) write(iout6,1073)
     if (ior.eq.4) write(iout6,1074)
  else
     if (ior.eq.1) write(iout6,1081)
     if (ior.eq.2) write(iout6,1082)
     if (ior.eq.3) write(iout6,1083)
     if (ior.eq.4) write(iout6,1084)
  endif


  write(iout6,*)
  write(iout6,1060)
  do iorb=1,norb
     write(*,1061) iorn(iorb),bond(iorb),gut(iorb),maxsororb(iorb),maxsorpot(iorb)
  enddo

  write(iout6,*)
  write(iout6,1090)
  do ig=1,ngrids
     !         write(iout6,1092) ig,ovforb(ig),ovfcoul(ig),ovfexch(ig)
     write(iout6,1091) ovforb(ig),ovfcoul(ig),ovfexch(ig)
  enddo




! FIXME
  if (iplot.eq.1) then
     write(iout6,*)
     write(iout6,1105)
  elseif (iplot.eq.2) then
     write(iout6,*)
     write(iout6,1106)
  endif

  if (ifermi.eq.1) then
     ip=iprint(230)
     iprint(230)=1
     call fermi
     iprint(230)=ip
  endif

  write(iout6,*)
  write(iout6,'("   machine accuracy      = ",e11.2)') precis
  write(iout6,*)

  !      write(iout6,'("   constants: ",/,
  !     &              "               pi        = ",e25.16,/,
  !     &              "               bohr      = ",e25.16," angstroms",/,
  !     &                                          )') pii,bohr2ang

  if (lengthfp.eq.8) then
     write(iout6,'("   constants: "/&
          & "               pi        = ",e25.16/&
          & "               bohr      = ",e25.16," angstroms"/&
          & )') pii,bohr2ang
  else
     write(iout6,'("   constants: "/ &
          & "               pi        = ",e45.34/ &
          & "               bohr      = ",e45.34," angstroms"/ &
          & )') pii,bohr2ang
  endif


  write(iout6,*)
  write(iout6,1110)

  !     Text and data segments of the program have about 700 and 210 KB,
  !     respectively.

  !     The amount of memory required by the program depends on a given
  !     case and the memory is allocated dynamically (if alloc routine
  !     is supported by the system).

  !     The amount of 'statically' allocated memeory (i.e. for a given
  !     choice of maxnu and maxmu values) is due to several integer and
  !     real arrays depending on the current values of maxnu and maxmu.

  maxorb2=(maxorb*maxorb)/2+maxorb

  !     real arrays
  isizereal= (35+9*maxorb+4*maxorb2)*maxnu+(41)*maxmu+(10+4*maxorb+maxorb/2)*maxorb &
       +18*maxorb+8*maxorb*maxorb+(maxbasis+maxorb+5)*maxbasis+44000+maxorb*maxorb

! FIXME see inihf 1830 arrays are counted as 1860
  isizeint=6*maxorb*maxorb+11*maxorb2 +(25)*maxorb+4*maxbasis

  length0=lengthfp*isizereal+lengthint*isizeint

  lengtht=8*(length1+length2+length3+length4+length5)+lengthint*(length6)

  if (imethod.eq.2.or.imethod.eq.3) then
     nexch=1
  endif

  !     1MB=2^20B
  omb=2.0_PREC**20
  write(iout6,1200) dble(length0)/omb,dble(lengtht)/omb

  !     write(*,1202) 2*maxfock*8/1.0e6
  mm3=2*norb*mxsize+nexch*mxsize
  if (idbg(550).ne.0) mm3=3*mxsize
  !      write(*,1210) norb,norb,nexch,i1e(norb),i1e(norb)*8/omb,
  !     &     i2e(norb),i2e(norb)*8/omb,
  !     &     nexch*mxsize,nexch*mxsize*8/omb,mm3,mm3*8/omb

  write(iout6,1210) dble(2*i1e(norb)*lengthfp/omb),dble(nexch*mxsize*lengthfp/omb)

  call separator


1000 format(10x,a2,'(',f5.2,')',3x,a2,' (',f5.2,')',3x,'R =',f9.5,' bohr =',f8.5,' angstroms')
1401 format(/3x,'finite nuclei (Fermi nuclear charge distribution):',/10x,'atomic masses:', d16.9,' ( ',a2,')',d16.9,' ( ',a2,')')
1402 format(/3x,'finite nuclei (Gauss nuclear charge distribution):',/10x,'atomic masses:', d16.9,' ( ',a2,')',d16.9,' ( ',a2,')')
1001 format(3x,'method: HFS','  (alpha = ',f10.5,')')
1002 format(3x,'method: DFT and ')
1003 format(3x,'                 ',a4,' functional',' (alpha = ',f10.5,')')
1004 format(3x,'                 ',a4,' functional')
1005 format(3x,'                 ',a4,' functional')
1006 format(3x,'method: SCMC')
1040 format(10x,'nu (h_nu)  = ',i4,'  (',f7.5,')'/10x,'mu (h_mu)  = ',i4,'  (',f7.5,')'/,10x,'R_infinity =  ',f6.2)
1042 format(10x,'nu (h_nu) ',i6,'    (',f7.5,')'/10x,'mu (h_mu) ',i6,'              '/,10x,'R_infinity   ',f6.2)
1050 format(10x,'thresholds',/14x,'scf iterations           =',i7/14x,'orbital energy           = 1.00E-',i2/ &
          14x,'orbital norm             = 1.00E-',i2/14x,'multipole moments recalc = ',1p,e8.2)
1051 format(10x,'orbitals are relaxed')
1052 format(10x,'orbitals are kept frozen')
1053 format(10x,'Coulomb potentials are relaxed')
1054 format(10x,'Coulomb potentials are kept frozen')
1055 format(10x,'exchange potential for each pair of orbitals is relaxed twice per single scf iteration')
1056 format(10x,'exchange potential for each pair of orbitals is relaxed once per single scf iteration')
1057 format(10x,'exchange potentials are kept frozen')
1058 format(/,10x,'multipole expansion coefficients = ',i2)

! 1060 format(10x,i2,' (mc)sor iterations on a (sub)grid without',
!     &     ' boundary values',/,10x,'being updated (default 10)',
!     &     //10x,i2,' iteration(s) of the group of (mc)sor iterations',
!     &     ' (as defined above);'/10x,'between these iterations boundary ',
!     &     'values are recalculated')
! 1060 format(10x,i2,' (mc)sor iterations for orbitals (default 10)',
!     &      /10x,i2,' (mc)sor iterations for potentials (default 10)')
 1060 format(28x,'(mc)sor iterations',/27x,' orbital  potentials')
 1061 format(10x,i4,1x,a8,a1,1x,i8,i10)



 1071 format(10x, 'ordering: natural column-wise')
 1072 format(10x, 'ordering: middle ')
 1073 format(10x, 'ordering: natural row-wise')
 1074 format(10x, 'ordering: reversed natural row-wise')
 1081 format(10x, 'ordering: natural column-wise + forward/backward sweeps')
 1082 format(10x, 'ordering: middle + forward/backward sweeps')
 1083 format(10x, 'ordering: natural row-wise + forward/backward sweeps')
 1084 format(10x, 'ordering: reversed natural row-wise + forward/backward sweeps')

! 1090 format(10x, 'overrelaxation parameters:'/13x,'subgrid      orbitals       potentials ')
 1090 format(10x, 'overrelaxation parameters:','   orbitals       potentials ')
 1091 format(38x,2x,f5.3,7x,f5.3,3x,f5.3)
! 1092 format(13x,2x,i2,10x,f5.3,6x,f5.3,6x,f5.3)
!c 1100 format(3x,'Xalpha approximation:'/10x,'alpha = ',f10.5)

! 1100 format(11x,'alpha = ',f10.5)


 1105 format(3x,'Plot:'/,10x,'orbital values exported in (x,z) ','cartesian coordinates')
 1106 format(3x,'Plot:'/,10x,'orbital values exported in (ni,mu) ','prolate spheroidal coordinates')
 1110 format(3x,'memory usage:')
 1200 format(10x,'text+data                           ',9x,'  0.9 MB '/,&
           10x,'bss (common block contributions)    ',6x,f8.1,' MB'/,&
           10x,'dynamical allocation                ',6x,f8.1,' MB ')
 1210 format(&
           14x,'orbitals + Coulomb potentials         ',f8.1,' MB',&
          /14x,'exchange potentials                   ',f8.1,' MB')

 1300 format(11x,'finite electric field: ',1Pe9.2,' au')
 1310 format(11x,'2D harmonic potential: ',1Pe9.2,' au')
  return
end subroutine printCase
