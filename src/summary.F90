! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module summary
  use params
  use discrete
  use commons
  use printUtils
  implicit none

contains
  ! ### locenergy ###
  !
  !     Calculates the locale energy for a given orbital.
  !
  subroutine locenergy 
    use params
    use discrete
    use scfshr
    use commons
    use utils
    use blas
    use inout
    use sharedMemory
    use utils 
    implicit none
    integer (KIND=IPREC) :: i,i1beg1,i1beg,i2beg,i2beg1,i3beg,ibeg,ihc,im,imax,in,iorb,iorb1,ioutmat, &
         ipc,isym,kex,nmut
    real (PREC) :: oc,w,w1,w2,wk2max,wtwoel
    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f4,wgt1,wgt2,wk0,wk1,wk2,wk3

    character*13 :: fn

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)


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
    mxsize=i1si(iorb)
    mxsize=i2si(iorb)

    !ipc=iorb+(iorb-1)*norb
    !engo(ipc)=zero
    ee(iorb,iorb)=zero
    isym=isymOrb(iorb)

    !    calculate derivatives over mu and ni

    call putin (nni,nmut,isym,psi(i1beg:),wk3)
    call diffnu (nmut,wk3,wk0,wk1,wk2)
    call putout (nni,nmut,wk1,wk0)

    call diffmu (nmut,wk3,wk2)
    call putout (nni,nmut,wk0,wk2)

    !    add contribution from derivatives over mu and ni

    call add (mxsize,wk0,wk1)

    if (mm(iorb).eq.0) then
       call dcopy (mxsize,f0,ione,wk0,ione)
    else
       w=dble(mm(iorb)*mm(iorb))
       call dcopy (mxsize,f0,ione,wk0,ione)
       call daxpy (mxsize,w,e,ione,wk0,ione)
    endif

    call proda (mxsize,psi(i1beg:),wk0,wk1)
    call prod (mxsize,wgt1,wk1)

    ! now wk1 contains (T+V)|psi>

    do i=1,mxsize
       wk2(i)=zero
    enddo

    if (nel.gt.1) then
       if (norb.eq.1) then
          ! contribution from coulomb potential of one sigma orbital
          call dcopy (mxsize,excp(i2beg:),ione,wk2,ione)
          call prod  (mxsize,psi(i1beg:),wk2)
       else
          do i=1,mxsize
             wk2(i)=zero
          enddo
          ! add contributions from coulomb and exchange potentials
          do iorb1=1,norb
             i1beg1=i1b(iorb1)
             i2beg1=i2b(iorb1)
             oc=occ(iorb1)
             kex=iorb+norb*(iorb1-1)

             if (iorb.le.iorb1) then
                ihc=iorb+iorb1*(iorb1-1)/2
             else
                ihc=iorb1+iorb*(iorb-1)/2
             endif

             i3beg=i3b(ihc)
             do i=1,mxsize
                wk0(i)=zero
             enddo

             call dcopy (mxsize,excp(i2beg1:),ione,wk0,ione)
             call prod  (mxsize,psi(i1beg:),wk0)

             if (iorb.eq.iorb1) oc=oc-1.0_PREC
             if (abs(oc-one)>epsilon(zero)) then
                call dscal (mxsize,oc,wk0,ione)
             endif

             if (iorb1.ne.iorb)  then
                oc=gec(kex)
                call prodas (mxsize,-oc,psi(i1beg1:),excp(i3beg:),wk0)
                if (ilc(ihc).gt.1) then
                   oc=gec(kex+norb*norb)
                   call prodas (mxsize,-oc,psi(i1beg1:),excp(i3beg+mxsize:),wk0)
                endif
             else
                if ((mm(iorb).gt.0).and.(ilc(ihc).gt.0)) then
                   oc=gec(kex)
                   call prodas (mxsize,-oc,psi(i1beg1:),excp(i3beg:),wk0)
                endif
             endif

             if (iorb1.eq.1) then
                call dcopy (mxsize,wk0,ione,wk2,ione)
             else
                call add   (mxsize,wk0,wk2)
             endif
          enddo
       endif

       ! to complete the integrand wk2 has to be multiplied by psi(i1beg)
       call prod (mxsize,wgt2,wk2)
    endif

    ! wk2 contains V(fock)|psi>

    i1beg=i1b(iorb)
    mxsize=i1si(iorb)

    call add   (mxsize,wk1,wk2)

    call dcopy (mxsize,wk2,ione,wk3,ione)

    w=-ee(iorb,iorb)
    
    call dcopy(mxsize,f4,ione,wk1,ione)
    call prod(mxsize,psi(i1beg:),wk1)
    call prod(mxsize,wgt2,wk1)
    call dcopy(mxsize,wk1,ione,wk0,ione)
    call daxpy(mxsize,w,wk1,ione,wk2,ione)

    w=zero
    w1=zero
    w2=zero
    wk2max=zero
    do i=1,mxsize
       ! exclude grid points where nuclei are located
       if (abs(wk2(i)).gt.wk2max) then
          wk2max=abs(wk2(i))
          imax=i
       endif
       if (.not.(i.eq.1.or.i.eq.(mxsize-mxnmu))) then
          w=w+wk2(i)*wk2(i)
          w1=w1+wk2(i)*wk2(i)*f4(i)*wgt2(i)
       endif
    enddo

    write(*,*) 'imax, wk2max',imax,wk2max

    ! FIXME
    if (idebug(496).ne.0) then
       go to ( 11, 12, 13, 14, 15, 16), iorb
00011  fn='bf-1p.lenergy'
       goto 100
00012  fn='bf-5s.lenergy'
       goto 100
00013  fn='bf-4s.lenergy'
       goto 100
00014  fn='bf-3s.lenergy'
       goto 100
00015  fn='bf-2s.lenergy'
       goto 100
00016  fn='bf-1s.lenergy'
00100  continue
       open(ioutmat,file=fn,status='replace',form='formatted')
       ibeg=i1b(iorb)
       call prtmat (nni,i1mu(1),wk2,ioutmat)
       close(ioutmat)
    endif

    write(*,1010) iorn(iorb),bond(iorb),w1,sqrt(w1),sqrt(w1)/(-ee(iorb,iorb))
01010 format(10x,i2,1x,a5,'  ',4e15.5)

    ! FIXME
    !call leexact1(wk3,wk0)

  end subroutine locenergy

  ! ### leexact1 ###
  !
  subroutine leexact1 (vt,orb)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: idev2,imu,ini,iorb,ipb,ipoints
    real (PREC) :: c1,d1,ddnom,dev1,dev2,dnum,e,r1,r1max,s1,s2,s3,t,w
    real (PREC), dimension(nni,mxnmu) :: orb,vt

    iorb=1
    !w=eng(iorb)
    w=ee(iorb,iorb)
    dev1=zero
    dev2=zero
    idev2=0
    r1max=0.20_PREC
    ipoints=0

    do ini=1,nni
       do imu=1,mxnmu
          r1=(r/2.0_PREC)*(vxi(imu)+veta(ini))
          dnum=zero
          ddnom=zero
          s1=zero
          s2=zero
          s3=zero
          if (abs(r1).gt.precis) then
             do ipb=npbasis,1,-1
                d1=primexp(ipb)
                c1=primcoef(iorb,ipb)
                e=(2.0_PREC*d1/pii)**(3.0_PREC/4.0_PREC)*exp(-d1*r1*r1)
                s1=s1+c1*(-2.0_PREC*d1*d1*r1*r1)*e
                s2=s2+c1*(3.0_PREC*d1)*e
                if (abs(r1).gt.precis) s3=s3+c1*(-z1/r1)*e
                ddnom=ddnom+c1*e
             enddo
             dnum=s1+s2+s3
             if (idebug(497).ne.0.and.abs(ddnom).gt.precis) then
                t=(vt(ini,imu)-dnum/ddnom*orb(ini,imu))**2
                dev1=dev1+(vt(ini,imu)-dnum/ddnom*orb(ini,imu))**2
             else
                t=0.0_PREC
             endif
             if (abs(r1).gt.r1max) then
                dev2=dev2+t
                idev2=idev2+1
             endif
          endif
       enddo
    enddo
    write(*,*) 'deviation from exact local energy',sqrt(dev1)
    write(*,*) 'deviation from exact local energy outside radius'
    write(*,*) idev2,r1max,sqrt(dev2)

  end subroutine leexact1

  
  ! ### tden ###
  !
  !    This routine initializes arrays used for differentiation needed in nuclder
  subroutine tden(iorb,ngorb1,psi,wk2)
    implicit none
    integer (KIND=IPREC) :: i,iborb1,iorb,iorb1,ngorb1
    real (PREC) :: coo
    real (PREC), dimension(*) :: psi,wk2

    do i=1,ngorb1
       wk2(i)=.0_PREC
    enddo

    call extinorg(psi)

    !      do iorb1=1,norb
    do iorb1=iorb,iorb
       iborb1=i1b(iorb1)
       ngorb1=i1si(iorb1)
       coo=occ(iorb1)
       do i=1,ngorb1
          wk2(i)=wk2(i)+coo*psi(iborb1+i-1)*psi(iborb1+i-1)
       enddo
    enddo

  end subroutine tden

  ! ### extinorg ###
  subroutine extinorg(psi)
    implicit none
    integer (KIND=IPREC) :: i,iprt
    real (PREC), dimension(nni,*) :: psi

    iprt=1
    if (iprt.eq.0) then
       write(*,*) 'psi(nni-5..,1) ',(psi(i,1),i=nni-5,nni)
       write(*,*) 'psi(nni,1..)   ',(psi(nni,i),i=1,6)
    endif

    psi(1,1)= exeven(1)*psi(2,1)+exeven(2)*psi(3,1)+exeven(3)*psi(4,1)+exeven(4)*psi(5,1)+exeven(5)*psi(6,1)

    psi(nni,1)= exeven(1)*psi(nni-1,1)+exeven(2)*psi(nni-2,1)+exeven(3)*psi(nni-3,1) &
         +exeven(4)*psi(nni-4,1)+exeven(5)*psi(nni-5,1)

    iprt=1
    if (iprt.eq.0) then
       write(*,*)
       write(*,*) 'psi(nni-5..,1) ',(psi(i,1),i=nni-5,nni)
       write(*,*) 'psi(nni,1..)   ',(psi(nni,i),i=1,6)
    endif

  end subroutine extinorg

 
  ! ### radialden ###
  !
  !     Calculates total radial densities relative to centres A (z-R/2)
  !     and B (z+R/2) along the internuclear axis (from A to -\infty or B
  !     to +\infty).
  !
  subroutine radialden
    use sharedMemory
    use utils

    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb,iunit,ngorb
    real (PREC) :: coo
    real (PREC), dimension(:), pointer :: psi,wk0

    psi=>orbptr
    wk0 =>scratchptr(1:mxsize8)

    ngorb=i1si(norb)
    call zeroArray(ngorb,wk0)

    do iorb=1,norb
       iborb=i1b(iorb)
       coo=occ(iorb)
       do i=1,ngorb
          wk0(i)=wk0(i)+coo*psi(iborb+i-1)*psi(iborb+i-1)
       enddo
    enddo

    iunit=99
#ifdef PRINT
! print=110: radialDen: radial density relative to centre A along the internuclear axis    
    ! -R_{\infty}<=z<=-R/2
    if (iprint(110)==1) then
       !        print radial density relative to centre A along the
       open(iunit,file='density-A',status='replace',form='formatted')
       call prtdenA(nni,i1mu(1),wk0,iunit)
       close(iunit)
    endif
#endif

#ifdef PRINT
! print=111: radialDen: radial density relative to centre B along the internuclear axis    
    ! R/2<=z<=R_{\infty}
    if (iprint(111)==1) then
       open(iunit,file='density-B',status='replace',form='formatted')
       call prtdenB(nni,i1mu(1),wk0,iunit)
       close(iunit)
    endif
#endif
  end subroutine radialden

  ! ### prtdenA ###
  !
  !     FIXME
  !
  subroutine prtdenA (m,n,a,ioutmat)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: in,imu,ioutmat,m,n
    real (PREC) :: r1t
    real (PREC), dimension(m,n) :: a

    write(ioutmat,'(10x,"r(au)",12x,"total electronic density")')
    in=m
    do imu=1,n
       r1t=(half*r)*(vxi(imu)+veta(in))
       write(ioutmat,'(2E25.16)') r1t,a(in,imu)
    enddo

  end subroutine prtdenA

  ! ### prtdenB ###
  !
  !     FIXME
  !
  subroutine prtdenB (m,n,a,ioutmat)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: in,imu,ioutmat,m,n
    real (PREC) :: r2t
    real (PREC), dimension(m,n) :: a

    write(ioutmat,'(10x,"r(au)",12x,"total electronic density")')
    in=1
    do imu=1,n
       r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))
       write(ioutmat,'(2E25.16)') r2t,a(in,imu)
    enddo

  end subroutine prtdenB

  subroutine showTail(psi)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,in,imu,iborb,iorb
    real (PREC) :: tailmax
    real (PREC), dimension(*) :: psi

    !write(*,'(/"  largest values of orbitals in their tail regions, i.e. max(orb(i,mxnmu),i=1,nni):")')
    write(*,'(/"  largest values of orbitals in their tail regions:")')
    imu=mxnmu-1
    do iorb=norb,1,-1
       tailmax=zero
       iborb=i1b(iorb)
       do in=1,nni
          if (abs(psi(iborb+(imu-1)*nni+in-1)) > tailmax) then
             tailmax=abs(psi(iborb+(imu-1)*nni+in-1))
          endif
       enddo
       write(*,'("    ",i4,1x,a8,a1,2x,e12.2)') iorn(iorb),bond(iorb),gusym(iorb),tailmax
    enddo
    write(*,'()') 
  end subroutine showTail

  
  ! ### printResults ####
  !
  !     Prints total energy and its components, orbital energies, multipole
  !     moments, etc upon completion of the SCF process
  !
  subroutine printResults 
    use params
    use discrete
    use commons
    use blas
    use totalEnergy
#ifdef LIBXC    
    use totalEnergyLXC
#endif    
    use normOrtho
    use printInitData    
    use scfshr
    use scfUtils
    use sharedMemory
    implicit none
    integer (KIND=IPREC) :: i,ibeg,iorb,inuclder,jorb,k,ngrid

    real (PREC) :: etot_final,tcpu,trelax,xnorm,zarea
    real (PREC), dimension(10) :: totalref
    real (PREC), dimension(:), pointer ::  cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch,cw_scratch4lxc
    
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    cw_orb=>orbptr
    cw_coul=>exchptr
    cw_exch=>exchptr
    cw_suppl=>supplptr
    cw_sctch=>scratchptr
    cw_scratch4lxc=>scratch4lxcptr
    !     calculate final total energy
    !     check orhogonality of orbitals

#ifdef PRINT
! print= 30: printResults: checking orthogonality of orbitals    
    if (iprint(30).ne.0) then
       write(*,'(/1x,"orthogonality of orbitals:")')
       do iorb=norb,1,-1
          call checkOrtho (.true.,iorb)
       enddo
    endif
    write (*,*)
#endif

    ! etotal and etotalg give the same results for He but not for Be 
    if (HF) call etotalHF

  

    if (DFT.or.HFS.or.SCMC) call etotalDFT
    
#ifdef LIBXC
    if (LXC) then
       enexchdft=zero
       call etotalLXC
    endif
#endif          

#ifdef PRINT
! print=590: printResults: calculating total energy with Gaussian charge distributin
    if (lpotGauss.and.iprint(590).ne.0) then
       call etotalGauss
       call printTotalEnergy
    endif
#endif
    
    etot_final=etot

    if (iprint16.eq.0) then
       write(*,*)
            write(*,'(5x,"nuclear attraction energy:        ",f19.12)') ennucel
            write(*,'(5x,"kinetic energy:                   ",f19.12)') enkin
            write(*,'(5x,"one-electron energy:              ",f19.12)') enkin+ennucel
            write(*,'(5x,"Coulomb energy:                   ",f19.12)') encoul
            write(*,'(5x,"exchange energy:                  ",f19.12)') enexch
            write(*,'(5x,"nuclear repulsion energy:         ",f19.12)') z1*z2/r
       if (HF) &
            write(*,'(5x,"Coulomb energy (DFT/LXC):         ",f19.12)') encouldft
       if (HF) &
            write(*,'(5x,"exchange energy (DFT/LXC):        ",f19.12)') enexchdft
       if (DFT) &
            write(*,'(5x,"exchange-correlation energy (DFT):",f19.12)')  enexchdft
       if (LXC) &
            write(*,'(5x,"exchange-correlation energy (LXC):",f19.12)') wtwoelLXC
       if (lxcHyb) &
            write(*,'(5x,"exact-exchange energy (HYB):      ",f19.12)')  enexchdft
    else
            write(*,*)
            write(*,'("     nuclear attraction energy       ",f29.22)') ennucel
            write(*,'("     kinetic energy                  ",f29.22)') enkin
            write(*,'("     one-electron energy             ",f29.22)') enkin+ennucel
            write(*,'("     Coulomb energy                  ",f29.22)') encoul
            write(*,'("     exchange energy                 ",f29.22)') enexch
            write(*,'("     nuclear repulsion energy        ",f29.22)') z1*z2/r
            write(*,'("     Coulomb energy (DFT/LXC):       ",f29.22)') encouldft
       if (HF) &
            write(*,'(5x,"Coulomb energy (DFT/LXC):         ",f29.22)') encouldft
       if (HF) &
            write(*,'(5x,"exchange energy (DFT/LXC):        ",f29.22)') enexchdft
       if (DFT) &
            write(*,'(5x,"exchange-correlation energy (DFT):",f29.22)') enexchdft
       if (LXC) &
            write(*,'(5x,"exchange-correlation energy (LXC):",f29.22)') wtwoelLXC
       if (lxcHyb) &
            write(*,'(5x,"exact-exchange energy (HYB):      ",f29.22)') enexchdft
    end if

    ! if (iprint16.eq.0) then
    !    write(*,*)
    !    write(*,'("     total energy contributions: ")')
    !    write(*,'("     nuclear attraction energy          = ",f20.12)') ennucel
    !    write(*,'("     kinetic energy                     = ",f20.12)') enkin
    !    write(*,'("     one-electron energy                = ",f20.12)') enkin+ennucel
    !    write(*,'("     Coulomb energy                     = ",f20.12)') encoul
    !    write(*,'("     exchange energy                    = ",f20.12)') enexch
    !    write(*,'("     nuclear repulsion energy           = ",f20.12)') z1*z2/r
    !    if (HF) write(*,'("     Coulomb energy (DFT/LXC)           = ",f20.12)') encouldft
    !    if (HF) write(*,'("     exchange energy (DFT/LXC)          = ",f20.12)') enexchdft
    !    !write(*,'("     correlation energy (DFT)           = ",f20.12)') edftcorr
    !    if (DFT) write(*,'("     exchange-correlation energy (DFT)  = ",f20.12)')  enexchdft
    !    if (LXC) write(*,'("     exchange-correlation energy (LXC)  = ",f20.12)') wtwoelLXC
    !    !write(*,'("     exact-exchange energy (HYB)        = ",f20.12)') alphaflxc*enexchdft       
    !    if (lxcHyb) write(*,'("     exact-exchange energy (HYB)        = ",f20.12)')  enexchdft
    ! else
    !    write(*,*)
    !    write(*,'(" Total energy contributions: ")')
    !    write(*,'("     nuclear attraction energy          = ",f29.22)') ennucel
    !    write(*,'("     kinetic energy                     = ",f29.22)') enkin
    !    write(*,'("     one-electron energy                = ",f29.22)') enkin+ennucel
    !    write(*,'("     Coulomb energy                     = ",f29.22)') encoul
    !    write(*,'("     exchange energy                    = ",f29.22)') enexch
    !    write(*,'("     nuclear repulsion energy           = ",f29.22)') z1*z2/r
    !    write(*,'("     Coulomb energy (DFT/LXC)           = ",f29.22)') encouldft
    !    if (HF) write(*,'("     Coulomb energy (DFT/LXC)           = ",f29.22)') encouldft
    !    if (HF) write(*,'("     exchange energy (DFT/LXC)          = ",f29.22)') enexchdft
    !    if (DFT) write(*,'("     exchange-correlation energy (DFT)  = ",f29.22)') enexchdft
    !    if (LXC) write(*,'("     exchange-correlation energy (LXC)  = ",f29.22)') wtwoelLXC
    !    !write(*,'("     exact-exchange energy (HYB)        = ",f29.22)') alphaflxc*enexchdft
    !    if (lxcHyb) write(*,'("     exact-exchange energy (HYB)        = ",f29.22)') enexchdft
    ! end if

    ! FIXME print DFT functional used
    ! write(*,'("     Coulomb energy (DFT)               = ",f29.22)') encouldft
    ! write(*,'("     exchange energy (DFT)              = ",f29.22)') enexchdft
    ! write(*,'("     correlation energy (DFT)           = ",f29.22)') edftcorr

    ! write final orbital energies and normalization factors
    write(*,*)
    ! iprint16=1
    if (iprint16.eq.0) then
       write(iout6,1000)
       do i=1,norb
          write(iout6,1002) iorn(i),bond(i),gusym(i),ee(i,i),1.0_PREC-orbNorm(i)
       enddo
    else
       write(iout6,7000)
       do i=1,norb
          write(iout6,7002) iorn(i),bond(i),gusym(i),ee(i,i),1.0_PREC-orbNorm(i)
       enddo
    endif
    write(*,*)

    ! check how errors in the orbital normalization factors influence
    ! the accuracy of the total energy
    if (lcheckTotEnergy) then 
       ! denormalize orbitals (they leave SCF process properly normalized)
       do iorb=1,norb
          zarea = sqrt(orbNorm(iorb))
          call dscal (mxsize,zarea,cw_orb(i1b(iorb):),ione)
       enddo
       
       ! calculate final total energy
       write (*,*)
       etot_final=etot
       
       if (HF) call etotalHF

       if (DFT.or.HFS.or.SCMC) call etotalDFT 

#ifdef LIBXC
       if (LXC) then
          if (lxcHyb) then
             call etotallxcHyb
          else
             call etotalLXC
          endif
       endif
#endif          
       
       write(iout6,6111) abs(etot_final-etot),abs(etot_final/etot-1.0_PREC)*100
       
       ! normalize orbitals again to bring them to the previous state
       do iorb=1,norb
          zarea = one/sqrt(orbNorm(iorb))
          call dscal (mxsize,zarea,cw_orb(i1b(iorb):),ione)
       enddo
    endif

#ifdef PRINT
! print=191: printResults: checking orthogonality of orbitals
    ! check orthohonality
    if (iprint(191).ne.0) then
       write(*,*) 'Checking orthogonality:'
       do iorb=norb,1,iorb
          call checkOrtho (.false.,iorb)
       enddo
    endif
#endif
    
#ifdef PRINT
! print=192: printResults: checking Euclidean norms of (T+V(n)+V-E)|i>
    if (iprint(192).ne.0) then
       write(*,*)
       write(*,*) 'Euclidean norms of (T+V(n)+V-E)|i>:'
       do iorb=1,norb
          call locenergy 
       enddo
    endif
#endif

#ifdef PRINT    
! print=193: printResults: checking multipole expansion contributions to potentials
    if (iprint(193).ne.0) then
       call checkPot (cw_coul,cw_exch)
    endif
#endif

    ! check symmetry of orbitals in homonuclear case
    !if (abs(z1-z2).lt.homolevl.and.lbreakCi) call checkSym(cw_orb)
    if (abs(z1-z2).lt.homolevl) call checkSym(cw_orb)
    

    ! calculate multipole moments and expectation values
    !if (lmmoments .and. .not.OED) then
    if (lmmoments) then
       !call propet2
       call propet 
    endif
    
#ifdef PRINT
    ! print=152: printResults: checking dependence of multipole moments on normalisation
    if (iprint(152).ne.0) then
       ! save total multipole moments
       do k=1,maxmpole
          totalref(k)=total(k)
       enddo
       
       ! denormalize orbitals (they leave SCF process properly normalized)
       do iorb=1,norb
          ibeg = i1b (iorb)
          do i=1,mxsize
             cw_sctch(i)=cw_orb(ibeg+i-1)
          enddo
          ngrid= i1si(iorb)
          xnorm=orbNorm(iorb)
          zarea = sqrt(xnorm)
          call dscal (ngrid,zarea,cw_sctch,ione)
       enddo
       
       ! calculate total multipole moments
       call propet3 
       
       ! normalize orbitals again to bring them to the previous state
       do iorb=1,norb
          ibeg = i1b (iorb)
          do i=1,mxsize
             cw_sctch(i)=cw_orb(ibeg+i-1)
          enddo
          
          ngrid= i1si(iorb)
          xnorm=orbNorm(iorb)
          zarea = 1.0_PREC/sqrt(xnorm)
          call dscal (ngrid,zarea,cw_sctch(ibeg:),ione)
       enddo
       
       ! print relative errors of multipole moments
       write(*,2000)
       do k=1,4
          if (abs(totalref(k)).gt.1d-12) then
             write(*,2002) k,abs(one - total(k)/totalref(k))
          else
             write(*,2002) k,zero
          endif
2000      format(//5x,"relative errors of moments due to orbital norms not being equal 1:")
2002      format(10x,"k=",i1,1P,e12.2)
       enddo
       
       ! restore total multipole moments
       do k=1,maxmpole
          total(k)=totalref(k)
       enddo
    endif

#endif
    
    ! FIXME
    ! calling scmc is not needed at the end of calculations
    if (iscmc.eq.1) then
       call alphaSCMC 
    endif

#ifdef PRINT
! print=110: printResults: total radial density relative to centre A (z-R/2) 
    if (iprint(110)/=0) then
       call radialden
    endif
#endif

#ifdef PRINT
! print=111: printResults: total radial density relative to centre B (z+R/2)
    if (iprint(111)/=0) then
       call radialden
    endif
#endif

    if (ltail) call showTail(cw_orb)
    
    call separator

    ! CPU time statistics

    trelax=trelaxCpuOrb+trelaxCpuCoul+trelaxCpuExch
    tcpu=trelax+tortho+trayl+tmomen

    if (tcpu<86400) then
       write(*,'(" CPU summary (sec):")')
       write(*,5109) tlagra+trayl
       write(*,5102) tortho
       write(*,5104) tmomen
       write(*,5108) ttoten

       write(*,6102) trelaxCpuOrb
       if (lcoulexch) then
          write(*,6106) trelaxCpuExch
       else
          write(*,6103) trelaxCpuCoul
          write(*,6104) trelaxCpuExch
       endif
       write(*,5101) trelax
       !if (iscf.ne.0) write(*,5106) tcpu/dble(iscf)
       write(*,5105) tcpu
       
       write(*,'(/" System clock summary (sec):")')
       write(*,6102) trelaxRealOrb
       if (lcoulexch) then
          write(*,6106) trelaxRealExch
       else
          write(*,6103) trelaxRealCoul
          write(*,6104) trelaxRealExch
       endif
       write(*,6105) trelaxRealOrb+trelaxRealCoul+trelaxRealExch
       if (maxscf>0 .and. (openmp .or. mcsorpt) ) &
            write(*,'(/" Speedup:",38x,f9.2)') trelax/(trelaxRealOrb+trelaxRealCoul+trelaxRealExch)    
              
    else
       write(*,'(" CPU summary (sec):")')
       write(*,5009) tlagra+trayl
       write(*,5002) tortho
       write(*,5004) tmomen
       write(*,5008) ttoten
       write(*,6002) trelaxCpuOrb
       if (lcoulexch) then
          write(*,6006) trelaxCpuExch
       else
          write(*,6003) trelaxCpuExch
          write(*,6004) trelaxCpuExch
       endif
       write(*,5001) trelax
       !if (iscf.ne.0) write(*,5006) tcpu/dble(iscf)
       write(*,5005) tcpu       
       write(*,'(/" System clock summary (sec):")')
       write(*,6002) trelaxRealOrb
       if (lcoulexch) then
          write(*,6006) trelaxRealExch
       else
          write(*,6003) trelaxRealCoul
          write(*,6004) trelaxRealExch
       endif
       write(*,6005) trelaxRealOrb+trelaxRealCoul+trelaxRealExch
       if (maxscf>0 .and. (openmp .or. mcsorpt) ) &
            write(*,'(/" Speedup:",38x,f12.2)') trelax/(trelaxRealOrb+trelaxRealCoul+trelaxRealExch)    
              
    endif

    
    
1000 format(8x,'orbital',17x,'energy',13x,'1-norm')
1002 format(4x,i3,1x,a8,1x,a1,1Pe28.16,1Pe12.2)
5001 format(4x,'relaxation of orbitals & potentials .......', f12.2)
5002 format(4x,'normalization+orthogonalization ...........', f12.2)
5009 format(4x,'Lagrange multipliers ......................', f12.2)
5004 format(4x,'multipole moments .........................', f12.2)
5005 format(4x,'SCF iterations ............................', f12.2)
5006 format(4x,'single SCF iteration ......................', f12.2)
5008 format(4x,'total energy ..............................', f12.2)
5101 format(4x,'relaxation of orbitals & potentials .......', f9.2)
5102 format(4x,'normalization+orthogonalization ...........', f9.2)
5109 format(4x,'Lagrange multipliers ......................', f9.2)
5104 format(4x,'multipole moments .........................', f9.2)
5105 format(4x,'SCF iterations ............................', f9.2)
5106 format(4x,'single SCF iteration ......................', f9.2)
5108 format(4x,'total energy ..............................', f9.2)
    
6002 format(4x,'relaxation of orbitals ....................', f12.2)
6003 format(4x,'relaxation of Coulomb potentials  .........', f12.2)    
6004 format(4x,'relaxation of exchange potentials .........', f12.2)
6005 format(4x,'relaxation of orbitals & potentials .......', f12.2)        
6006 format(4x,'relaxation of Coulomb & exchange potentials', f12.2)    

6102 format(4x,'relaxation of orbitals ....................', f9.2)
6103 format(4x,'relaxation of Coulomb potentials  .........', f9.2)    
6104 format(4x,'relaxation of exchange potentials .........', f9.2)
6105 format(4x,'relaxation of orbitals & potentials .......', f9.2)        
6106 format(4x,'relaxation of Coulomb & exchange potentials', f9.2)    
       
6111 format(5x,'total energy uncertainty due to orbital norms not being equal 1:'/,5x, &
          & 'absolute = +/-',e8.2,',  relative = +/-',e8.2,'%')
    !6120 format(1x,'virial ratio:            ',e28.16)
7000 format(14x,'orbital',22x,'energy',24x,'1-norm')
7002 format(9x,i3,1x,a8,1x,a1,1Pe44.32,1Pe13.2)
    !7100 format(1x,'total electronic energy: ',e44.32)
    !7110 format(1x,'total energy:            ',e44.32)
    !7111 format(1x,'total energy uncertaininty due to orbital norms not being equal 1:'/, &
    !          & '  absolute   +/-',e8.2,/,'  relative   +/-',e8.2,'%')
    !7120 format(1x,'virial ratio:            ',e44.32)

  end subroutine printResults

  ! ### propet ####
  !
  !     Calculates dipole, quadrupole, octopole and hexadecapole moments
  !     relative to the geometrical centre and the centre of mass.
  !
  subroutine propet
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use sharedMemory
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i1beg,igp,imu,in,inioff,iorb,izz1,izz2,k,mu,ni,ngorb

    real (PREC) :: atw1,atw2,cm2zz,costh,qq,qtot,qktot,qktot1,qkz1,qkz2,r1,rr,&
         sum1,sum2,sum3,xxplusyy,w,z,zcm

    real (PREC), dimension(:), pointer :: cw_orb,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    cw_orb=>orbptr
    cw_exch=>exchptr
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk4 =>scratchptr( 4*mxsize8+1: 5*mxsize8)
    wk5 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    
    !     1 bohr electron = 2.541765 Debye -- Gaussian94 User's Reference
    !     data au2Debye /2.541765d0/
    !     data bohr2ang /0.529177249d0/


    do k=1,maxmpole
       elect(k)=0.0_PREC
       electDA(k)=0.0_PREC
       total(k)=0.0_PREC
       totalDA(k)=0.0_PREC
    enddo

    zcm=0.0_PREC
    qtot=0.0_PREC

    do iorb=1,norb
       qq=-occ(iorb)
       qtot=qtot+qq
       elect(1)=elect(1)+qq*cmulti(iorb)
       elect(2)=elect(2)+qq*cmulti(iorb+norb)
       if (mpole.ge.3) elect(3)=elect(3)+qq*cmulti(iorb+2*norb)
       if (mpole.ge.4) elect(4)=elect(4)+qq*cmulti(iorb+3*norb)
       !write(*,'(i5,4e12.4)') iorb,cmulti(iorb),cmulti(iorb+norb),cmulti(iorb+2*norb),cmulti(iorb+3*norb)
       !write(*,'(/i5,4e25.16)') iorb,qq*cmulti(iorb),elect(1)
    enddo

    do k=1,mpole
       qkz1=z1*abs(r/2.0_PREC)**k*plegendg(k,izero,-1.0_PREC)
       qkz2=z2*abs(r/2.0_PREC)**k*plegendg(k,izero, 1.0_PREC)
       qktot1=elect(k)+qkz1+qkz2
       total(k)=qktot1
       totalDA(k)=total(k)*au2Debye*bohr2ang**(k-1)
       electDA(k)=elect(k)*au2Debye*bohr2ang**(k-1)
    enddo

    if (iprint16.eq.0) then
       write(*,'(/," multipole moments relative to geometrical centre:"/)')

       write(*,'(36x,"electronic (au/Debye-Ang^k) ",2x,"total (au/Debye-Ang^k)")')

       k=1
       write(*,'(5x,"dipole (Mu_z, k=1)          ",2e28.16)') elect(k),total(k)
       write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

       k=2
       write(*,*)
       write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e28.16)') elect(k),total(k)
       write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

       if (mpole.ge.3) then
          k=3
          write(*,*)
          write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e28.16)') elect(k),total(k)
          write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
       endif

       if (mpole.ge.4) then
          k=4
          write(*,*)
          write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e28.16)') elect(k),total(k)
          write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
       endif
    else
       write(*,'(/," multipole moments relative to geometrical centre:"/)')

       write(*,'(44x,"electronic (au/Debye-Ang^k) ",17x,"total (au/Debye-Ang^k)")')

       k=1
       write(*,'(5x,"dipole (Mu_z, k=1)          ",2e41.32)') elect(k),total(k)
       write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)

       k=2
       write(*,*)
       write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e41.32)') elect(k),total(k)
       write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)

       if (mpole.ge.3) then
          k=3
          write(*,*)
          write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e41.32)') elect(k),total(k)
          write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)
       endif

       if (mpole.ge.4) then
          k=4
          write(*,*)
          write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e41.32)') elect(k),total(k)
          write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)
       endif
    endif

    write(*,'(//," multipole moments relative to centre of mass:")')

    izz1=nint(z1)
    izz2=nint(z2)
    atw1=atweight(izz1)
    atw2=atweight(izz2)
    zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC
    write(iout6,1000) element(izz1),atweight(izz1),-r/2.0_PREC,element(izz2), atweight(izz2),r/2.0_PREC,zcm

01000 format(/,6x,'centre     atomic weight       z',/8x,a2,4x,f14.8,3x,f10.4,/, &
         8x,a2,4x,f14.8,3x,f10.4,/,2x,'centre-of-mass',15x,f10.4/)

    ! multipole moments relative to the centre of mass

#ifdef PRINT
! print=214: printResults: round-off errors in orbital contributions to multipole moments 
    
    ! estimate round-off errors present in orbital contributions to
    ! multipole moments of a given order k
    if (iprint(214).ne.0) then
       write(iout6,'(4x,"k",5x,"orb",17x,"contribution",9x,"error" )')
    endif
#endif
    
    do k=1,mpole
       !  loop over grid points
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=inioff+in

             !              rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
             !              x=(r/2.0_PREC)*vxi1(imu)*veta1(in)

             z=(r/2.0_PREC)*vxi(imu)*veta(in)-zcm
             xxplusyy=((r/2.0_PREC)*vxi1(imu)*veta1(in))**2
             rr=sqrt(xxplusyy+z*z)

             if (abs(rr).lt.precis) then
                costh=0.0_PREC
             else
                costh=z/rr
             endif
             wk0(igp)=rr**k*plegendg(k,izero,costh)
          enddo
       enddo

       qtot=0.0_PREC
       cm2zz=0.0_PREC
       do iorb=1,norb
          i1beg=i1b(iorb)
          ngorb=i1si(iorb)
          qq=-occ(iorb)
          qtot=qtot+qq
          call prod2 (ngorb,cw_orb(i1beg:),cw_orb(i1beg:),wk3)

          nni2=nni/2
          sum1=0.0_PREC
          sum2=0.0_PREC
          do imu=mxnmu,1,-1
             inioff=(imu-1)*nni
             do in=1,nni2
                igp=inioff+in
                sum1=sum1+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
             enddo

             do in=nni2+1,nni
                igp=inioff+in
                sum2=sum2+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
             enddo
          enddo

          sum3=0.0_PREC
          do imu=1,mxnmu
             inioff=(imu-1)*nni
             do in=1,nni
                igp=inioff+in
                wk5(igp)=wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
                sum3=sum3+wk5(igp)
             enddo
          enddo

#ifdef PRINT
! print=214: printResults: round-off errors in orbital contributions to multipole moments 
          if (iprint(214).ne.0) then
             write(iout6,'(i5,i4,1x,a8,a1,e28.16,e10.2)') k, iorn(iorb),bond(iorb),gusym(iorb),sum3,sum1+sum2-sum3
          endif
#endif
          cm2zz=cm2zz+qq*(sum1+sum2)
       enddo

       elect(k)=cm2zz
       electDA(k)=elect(k)*au2Debye*bohr2ang**(k-1)

       qkz1=z1*abs(r/2.0_PREC+zcm)**k*plegendg(k,izero,-1.0_PREC)
       qkz2=z2*abs(r/2.0_PREC-zcm)**k*plegendg(k,izero,1.0_PREC)
       qktot=cm2zz+qkz1+qkz2
       total(k)=qktot
       totalDA(k)=total(k)*au2Debye*bohr2ang**(k-1)
    enddo

    if (iprint16.eq.0) then
       write(*,'(/,36x,"electronic (au/Debye-Ang^k) ",2x,"total (au/Debye-Ang^k)")')

       k=1
       write(*,'(5x,"dipole (Mu_z, k=1)          ",2e28.16)') elect(k),total(k)
       write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

       k=2
       write(*,*)
       write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e28.16)') elect(k),total(k)
       write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)

       if (mpole.ge.3) then
          k=3
          write(*,*)
          write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e28.16)') elect(k),total(k)
          write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
       endif

       if (mpole.ge.4) then
          k=4
          write(*,*)
          write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e28.16)') elect(k),total(k)
          write(*,'(5x,"                            ",2e28.16)') electDA(k),totalDA(k)
       endif
    else
       write(*,'(/,41x,"electronic (au/Debye-Ang^k) ",17x,"total (au/Debye-Ang^k)")')

       k=1
       write(*,'(5x,"dipole (Mu_z, k=1)          ",2e41.32)') elect(k),total(k)
       write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)

       k=2
       write(*,*)
       write(*,'(5x,"quadrupole (Theta_zz, k=2)  ",2e41.32)') elect(k),total(k)
       write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)

       if (mpole.ge.3) then
          k=3
          write(*,*)
          write(*,'(5x,"octopole (Omega_zzz, k=3)   ",2e41.32)') elect(k),total(k)
          write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)
       endif

       if (mpole.ge.4) then
          k=4
          write(*,*)
          write(*,'(5x,"hexadecapole (Phi_zzzz, k=4)",2e41.32)') elect(k),total(k)
          write(*,'(5x,"                            ",2e41.32)') electDA(k),totalDA(k)
       endif
    endif

#ifdef PRINT
! print=216: printResults: expectation values of r_1^k r and r^2
    if (iprint(216).ne.0) then
       write(*,*)
       write(*,*)
       write(*,*) 'Expectation values of r_1^k'

       do iorb=1,norb
          i1beg=i1b(iorb)
          ngorb=i1si(iorb)
          write(*,*)
          !write(*,*) 'orbital #',iorb
          write(*,'(27x,i4,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)
          do mu=1,mxnmu
             imu=(mu-1)*nni
             do ni=1,nni
                r1=(r/2.0_PREC)*(vxi(mu)+veta(ni))
                if (r1.lt.precis) then
                   wk0(imu+ni)=0.0_PREC
                else
                   wk0(imu+ni)=1.0_PREC/r1
                endif
                wk1(imu+ni)=r1
                wk2(imu+ni)=r1*r1
             enddo
          enddo
          call prod2 (ngorb,cw_orb(i1beg:),wk0,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<1/r> =",1Pe23.16)') w

          call prod2 (ngorb,cw_orb(i1beg:),wk1,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"< r > =",1Pe23.16)') w          

          call prod2 (ngorb,cw_orb(i1beg:),wk2,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<r^2> =",1Pe23.16)') w                    
       enddo

       write(*,*)
       write(*,*)
       write(*,*) 'Expectation values of r^k'

       do iorb=1,norb
          i1beg=i1b(iorb)
          ngorb=i1si(iorb)
          write(*,*)
          write(*,'(27x,i4,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)          
          do mu=1,mxnmu
             imu=(mu-1)*nni
             do ni=1,nni
                rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)*r/2.0_PREC
                if (rr.lt.precis) then
                   wk0(imu+ni)=0.0_PREC
                else
                   wk0(imu+ni)=1.0_PREC/rr
                endif
                wk1(imu+ni)=rr
                wk2(imu+ni)=rr*rr
             enddo
          enddo
          call prod2 (ngorb,cw_orb(i1beg:),wk0,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<1/r> =",1Pe23.16)') w                    

          call prod2 (ngorb,cw_orb(i1beg:),wk1,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"< r > =",1Pe23.16)') w                              

          call prod2 (ngorb,cw_orb(i1beg:),wk2,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<r^2> =",1Pe23.16)') w                                        

          call prod2 (ngorb,cw_orb(i1beg:),cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<   > =",1Pe23.16)') w                                                  
       enddo

       write(*,*)
       write(*,*)
       write(*,*) 'Expectation values of r_2^k'

       do iorb=1,norb
          i1beg=i1b(iorb)
          ngorb=i1si(iorb)
          write(*,*)
          write(*,'(27x,i4,1x,a8,a1)') iorn(iorb),bond(iorb),gusym(iorb)                    
          do mu=1,mxnmu
             imu=(mu-1)*nni
             do ni=1,nni
                r1=(r/2.0_PREC)*(vxi(mu)-veta(ni))
                if (r1.lt.precis) then
                   wk0(imu+ni)=0.0_PREC
                else
                   wk0(imu+ni)=1.0_PREC/r1
                endif
                wk1(imu+ni)=r1
                wk2(imu+ni)=r1*r1
             enddo
          enddo
          call prod2 (ngorb,cw_orb(i1beg:),wk0,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<1/r> =",1Pe23.16)') w                                                  

          call prod2 (ngorb,cw_orb(i1beg:),wk1,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"< r > =",1Pe23.16)') w                                                            

          call prod2 (ngorb,cw_orb(i1beg:),wk2,wk3)
          call prod  (ngorb,cw_orb(i1beg:),wk3)
          call prod  (ngorb,f4,wk3)
          w=ddot(ngorb,wgt2,ione,wk3,ione)
          write(*,'(15x,"<r^2> =",1Pe23.16)') w                                                                      
       enddo
    endif
#endif
    
    ! calculating total electronic density at nuclei

    elect(1)=0.0_PREC
    elect(2)=0.0_PREC
    do iorb=1,norb
       qq=-occ(iorb)
       i1beg=i1b(iorb)
       ngorb=i1si(iorb)
       elect(1)=elect(1)+qq*cw_orb(i1beg)*cw_orb(i1beg)
       elect(2)=elect(2)+qq*cw_orb(i1beg+nni-1)*cw_orb(i1beg+nni-1)
    enddo

    if (iprint16.eq.0) then
       write(*,'(/,5x,"total charge density at (0,0,-R/2) and (0,0,+R/2)",e25.16 )') 
       write(*,'(36x,e25.16 )') elect(2)
       write(*,'(36x,e25.16 )') elect(1)       
    else
       write(*,'(/,5x,"total charge density at (0,0,-R/2) and (0,0,+R/2)",e25.16 )') 
       write(*,'(30x,e44.32 )') elect(2)
       write(*,'(30x,e44.32 )') elect(1)       
    endif

  end subroutine propet


  ! ### propet3 ####
  !
  !     Calculates dipole, quadrupole, octopole and hexadecapole moments
  !     relative to the geometrical centre and the centre of mass.
  !
  subroutine propet3 
    use params
    use blas
    use discrete
    use scfshr
    use commons
    use sharedMemory
    use utils    
    implicit none

    integer (KIND=IPREC) :: i1beg,igp,imu,in,inioff,iorb,izz1,izz2,k,ngorb
    real (PREC) :: atw1,atw2,cm2zz,costh,qq,qtot,qktot,qkz1,qkz2,rr,&
         sum1,sum2,sum3,xxplusyy,z,zcm

    real (PREC), dimension(:), pointer :: cw_orb,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5

    !     1 bohr electron = 2.541765 Debye -- Gaussian94 User's Reference
    !     data au2Debye /2.5417650_PREC/
    !     data bohr2ang /0.5291772490_PREC/
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    cw_orb=>orbptr
    cw_exch=>exchptr
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk4 =>scratchptr( 4*mxsize8+1: 5*mxsize8)
    wk5 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    
    do k=1,maxmpole
       total(k)=0.0_PREC
    enddo

    zcm=0.0_PREC
    qtot=0.0_PREC


    izz1=nint(z1)
    izz2=nint(z2)
    atw1=atweight(izz1)
    atw2=atweight(izz2)
    zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC

    ! multipole moments relative to the centre of mass


    izz1=nint(z1)
    izz2=nint(z2)
    atw1=atweight(izz1)
    atw2=atweight(izz2)
    zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC

    ! multipole moments relative to the centre of mass

    do k=1,mpole

       !  loop over grid points

       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=inioff+in

             !              rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
             !              x=(r/2.0_PREC)*vxi1(imu)*veta1(in)

             z=(r/2.0_PREC)*vxi(imu)*veta(in)-zcm
             xxplusyy=((r/2.0_PREC)*vxi1(imu)*veta1(in))**2
             rr=sqrt(xxplusyy+z*z)

             if (abs(rr).lt.precis) then
                costh=0.0_PREC
             else
                costh=z/rr
             endif
             wk0(igp)=rr**k*plegendg(k,izero,costh)
          enddo
       enddo

       qtot=0.0_PREC
       cm2zz=0.0_PREC
       do iorb=1,norb
          i1beg=i1b(iorb)
          ngorb=i1si(iorb)
          qq=-occ(iorb)
          qtot=qtot+qq
          call prod2 (ngorb,cw_orb(i1beg:),cw_orb(i1beg:),wk3)

          nni2=nni/2
          sum1=0.0_PREC
          sum2=0.0_PREC
          do imu=mxnmu,1,-1
             inioff=(imu-1)*nni
             do in=1,nni2
                igp=inioff+in
                sum1=sum1+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
             enddo

             do in=nni2+1,nni
                igp=inioff+in
                sum2=sum2+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
             enddo
          enddo

          sum3=0.0_PREC
          do imu=1,mxnmu
             inioff=(imu-1)*nni
             do in=1,nni
                igp=inioff+in
                wk5(igp)=wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
                sum3=sum3+wk5(igp)
             enddo
          enddo

#ifdef PRINT
! print=214: propet3: k iorb sum3 sum1+sum2-sum3           
          if (iprint(214).ne.0) then
             write(iout6,'(i5,i4,1x,a8,a1,e28.16,e10.2)') k, iorn(iorb),bond(iorb),gusym(iorb),sum3,sum1+sum2-sum3
          endif
#endif
          cm2zz=cm2zz+qq*(sum1+sum2)
       enddo

       qkz1=z1*abs(r/2.0_PREC+zcm)**k*plegendg(k,izero,-1.0_PREC)
       qkz2=z2*abs(r/2.0_PREC-zcm)**k*plegendg(k,izero,1.0_PREC)
       qktot=cm2zz+qkz1+qkz2
       total(k)=qktot
    enddo

  end subroutine propet3
  
  subroutine prtdenAB (m,n,a,ioutmat)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,np,ioutmat,m,n
    real (PREC) :: r2t
    real (PREC), dimension(*) :: a

    write(ioutmat,'(10x,"#total electronic density along z-axis:R=",f5.2)') r

    np=0
    j=nni
    do i=mxnmu-4,2,-1
       np=np+1
       write(ioutmat,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),a((i-1)*nni+j)
    enddo

    i=1
    do j=nni,1,-1
       np=np+1
       write(ioutmat,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),a((i-1)*nni+j)
    enddo

    j=1
    do i=2,mxnmu-4
       np=np+1
       write(ioutmat,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),a((i-1)*nni+j)
    enddo
    

  end subroutine prtdenAB
  
  ! ### checkPot ###
  !
  !     Checks contributions of multipole moments to Coulomb and exchange
  !     potentials at the practical infinity.
  !
  !     Determines asymptotic (boundary) values of Coulomb and exchange
  !     potentials from the multipole expansion of a given order.
  !
  subroutine checkPot (pot,excp)
    use params
    use commons
    use scfUtils
    
    implicit none
    integer (KIND=IPREC) :: iax,ibeg,idel,ido,iorb,iorb1,iorb2,ipc

    real (PREC), dimension(*) :: pot,excp

    ! write(*,*) 'asympt: for two different nonsigma orbitals',&
    !     	  'the case (+-) is now also included'

    write(*,*)
    write(*,*) 'Checking multipole expansion for Coulomb potentials'
    do iorb=norb,1,-1
       ibeg = i1b (iorb)
       call coulAsymptCont(iorb,pot(ibeg))
    enddo

    if (DFT.or.HFS.or.SCMC) return

    write(*,*)
    write(*,*) 'Checking multipole expansion for exchange potentials'
    do iorb1=1,norb
       do iorb2=iorb1,norb
          if (iorb1.eq.iorb2.and.mgx(6,iorb1).eq.0 ) goto 10
          if ((iorb1.eq.iorb2).and.(ilc(iorb1*(iorb1+1)/2).lt.1)) goto 10

          !           orbitals in increasing order

          ipc=iorb1+iorb2*(iorb2-1)/2
          iax=i3b(ipc)

          idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
          if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)

          ido=0
1234      ido=ido+1
          if (ido.eq.2) then
             idel=mgx(6,iorb2)+mgx(6,iorb1)
             ipc=ipc+norb*(norb+1)/2
             iax=iax+i3si(ipc)
          endif

          write(*,1000) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb2),bond(iorb2),gusym(iorb2)
          write(*,*) '------ ilc(ipc) ',iorb1,iorb2,'...',idel,'...',ipc,ilc(ipc)
1000      format(/i4,1x,a8,a1,3x,i4,1x,a8,a1,3x)

          call exchAsymptCont (idel,ipc,excp(iax))

          if (ilc(ipc).eq.2.and.ido.eq.1) go to 1234
10        continue
       enddo
    enddo

  end subroutine checkPot

  ! ### coulAsymptCont ###
  !
  !     Evaluates multipole moment contributions to the boundary values of
  !     Coulomb potential for a given orbital at imu=mxnmu, inu=(nn1-1)/2 (Pi/2)
  !
  subroutine coulAsymptCont(iorb,pot)
    use params
    use discrete
    use scfshr
    use commons

    implicit none
    integer (KIND=IPREC) :: i,iorb,itt,j,kk,kxk,m,n
    real (PREC) :: costh,pe,rr,xr,xrr
    real (PREC), dimension(10) :: dome,pottmp
    real (PREC), dimension(*) :: pot

    ! potentials are calculated for ni=Pi/2

    j=mxnmu
    itt=(j-1)*nni
    i=(nni-1)/2
    kk=i+itt
    rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
    costh=veta(i)*vxi(j)/rr
    xr=1.0_PREC/(rr*r2)

    dome(1)=costh
    dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
    do n=2,mpole-1
       dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
    enddo

    pe=0.0_PREC
    xrr=xr
    do m=1,mpole
       xrr=xrr*xr
       kxk=iorb+(m-1)*norb
       pe=pe+cmulti(kxk)*dome(m)*xrr
       pottmp(m)=r2*vxi(j)*(pe+xr)
    enddo

    write(*,1000) iorn(iorb),bond(iorb),gusym(iorb),pottmp(1),(pottmp(m)-pottmp(m-1),m=2,mpole),pottmp(mpole)
1000 format(/i4,1x,a8,a1,3x,/4e13.5/5e13.5)
  end subroutine coulAsymptCont
  
  ! ### exchAsymptCont ###
  !
  !     Evaluates multipole moment contributions to the boundary values of
  !     exchange for at imu=mxnmu, inu=(nn1-1)/2 (Pi/2)
  !
  subroutine exchAsymptCont (idel,ipc,excp)
    use params
    use discrete
    use commons
    use scfUtils

    implicit none
    integer (KIND=IPREC) :: i,idel,ipc,j,m

    real (PREC), dimension(*) ::  excp
    real (PREC), dimension(maxmpole) :: excptmp,pe

    j=mxnmu
    i=(nni-1)/2
    call vexch(i,j,idel,ipc,pe)
    do m=1,mpole
       excptmp(m)= r2*vxi(j)*pe(m)
    enddo
    write(*,1000) excptmp(1),(excptmp(m)-excptmp(m-1),m=2,mpole),excptmp(mpole)
1000 format(4e16.8/5e16.8)

  end subroutine exchAsymptCont

  subroutine checkOrbSym (nmut,orb,ihsym)
    use params
    use discrete

    implicit none
    integer (KIND=IPREC) :: ihsym,imu,inu,nmut
    real (PREC), dimension(nni,nmut) :: orb
    ! imu=13 corresponds to the smallest grid that can be used
    imu=13
    inu=2
    if (orb(inu,imu)*orb(nni-inu+1,imu).gt.0.0_PREC) then
       ihsym= 1
    else
       ihsym=-1
    endif

  end subroutine checkOrbSym
  
  ! ### checkSym ###
  !
  !     Checks Ci symmetry of all orbitals
  !
  subroutine checkSym(psi)
    use params
    use commons

    implicit none
    integer (KIND=IPREC) :: iorb,ibeg,ihsym,nmut
    real (PREC), dimension(*) :: psi
    character*8 :: sigma,pi,delta,phi

    data sigma/'sigma'/,pi/'pi'/,delta/'delta'/,phi/'phi'/

    write(*,1000)
    do iorb=1,norb
       ibeg = i1b(iorb)
       nmut = i1mu(iorb)
       call checkOrbSym(nmut,psi(ibeg),ihsym)
       if (ihsym.eq.1) then
          if (bond(iorb).eq.sigma.or.bond(iorb).eq.delta) then
             write(*,1005) iorn(iorb),bond(iorb),gusym(iorb)
          else
             write(*,1010) iorn(iorb),bond(iorb),gusym(iorb)
          endif
       endif
       if (ihsym.eq.-1) then
          if(bond(iorb).eq.sigma.or.bond(iorb).eq.delta) then
             write(*,1010) iorn(iorb),bond(iorb),gusym(iorb)
          else
             write(*,1005) iorn(iorb),bond(iorb),gusym(iorb)
          endif
       endif
    enddo
01000 format(/,' checking symmetry of orbitals:'/,'          required    actual ')
01005 format(1x,i3,1x,a8,1x,a1,10x,'g')
01010 format(1x,i3,1x,a8,1x,a1,10x,'u')
  end subroutine checkSym
  
end module summary
