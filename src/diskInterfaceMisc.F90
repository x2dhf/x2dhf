! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

module diskInterfaceMisc
  implicit none

contains
  ! ### writeDisk4pair ### 
  !
  !     Writes orbital in formatted form for the x5dhf program to use.
  subroutine writeDisk4pair
    use blas
    use commons
    use data4II
    use discrete
    use params
    use printUtils
    use inout
    use scfshr
    use sharedMemory
    use solver
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i4,i1beg,ngrid,isym4nu,isym4mu,isymmetry
    real (PREC), dimension(:), pointer :: psi,excp,f0,f4,&
              wk0,wk1,wk2,wk3,r8mxsize,r8mxsize1
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    r8mxsize =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    r8mxsize1 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    
    open(9999,file='out4pair.dat', status='unknown',form='formatted')
    
    !  call getDateTime(datetime)
    write(9999,'(80a1)') title
    !  write(9999,'(a80)') datetime
    write(9999,'(" nnu nmu ")') 
    write(9999,formint) nni,nmu(1)
    write(9999,'(" r  rinf")') 
    write(9999,formfp64) r,rinf
    write(9999,'(" z1 z2 ")') 
    write(9999,formfp64) z1,z2
    write(9999,'(" norb nel  ")')
    write(9999,formint) norb,nel
    
    iorb=1
    write(9999,'(" ordering of orbitals ")')
    write(9999,'(1x,i3,1x,a8,1x,a1)') (iorn(iorb),bond(iorb),gusym(iorb),iorb=1,norb)
    
    write(9999,'(" orbital energies: ")')
    do iorb=1,norb
       write(9999,formfp64) ee(iorb,iorb)
    enddo
    
    do iorb=1,norb
       i1beg=i1b(iorb)
       ngrid=i1si(iorb)
       
       do i4=1,ngrid
          r8mxsize(i4)=excp(i1beg+i4-1)
       enddo
       
       if (OED.or.initFuncsOED) then
          do i4=1,ngrid
             r8mxsize(i4)=zero
          enddo
       endif
       
       ! multiply \tilde{V}_C by 2/(R\xi) to get V_C
       
       ! isymmetry=isymOrb(iorb) 
       ! orbital 1sigma is even, its derivatives are odd functions
       isymmetry=-1

       !write(*,'("Warning! no call to pop2pot!!!")')
       call pot2pot(r8mxsize,f4)
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," Coulomb potential" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," Coulomb potential" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! write(9999,'(/" orbital=",i3,1x,a8,1x,a1," Coulomb potential" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)
       
       isym4nu=1
       isym4mu=1
       call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize,wk3)
       call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
       call putout (nni,mxnmu,r8mxsize,wk0)
       r8mxsize1=r8mxsize
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)
       
       call diff1mu (mxnmu,wk3,wk2)
       call putout (nni,mxnmu,r8mxsize,wk2)
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)

       isym4nu=-1
       isym4mu=1
       call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize1,wk3)
       
       call diff1mu (mxnmu,wk3,wk2)
       call putout (nni,mxnmu,r8mxsize,wk2)
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu,mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu,mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)
       
       do i4=1,ngrid
          r8mxsize(i4)=psi(i1beg+i4-1)
       enddo
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," orbital function")') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," orbital function")') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)
       
       ! printout suitable for comparison with JM data
       ! call prtmatcw1(nni,mxnmu,r8mxsize,9999)
       
       isym4nu=1
       call putin1 (nni,mxnmu,isym4nu,isym4mu,psi(i1beg:),wk3)
       call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
       call putout (nni,mxnmu,r8mxsize,wk0)
       r8mxsize1=r8mxsize
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / nu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," / nu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)

       isym4mu=1
       call diff1mu (mxnmu,wk3,wk2)
       call putout (nni,mxnmu,r8mxsize,wk2)
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," / mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)
       
       isym4nu=-1
       isym4mu=1
       call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize1,wk3)
       call diff1mu (mxnmu,wk3,wk2)
       call putout (nni,mxnmu,r8mxsize,wk2)
       
       write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / nu,mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       call prtmatcw(nni,mxnmu,r8mxsize,9999)
       ! write(6,'(/" orbital=",i3,1x,a8,1x,a1," / nu,mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
       ! call prtmatrw(nni,mxnmu,r8mxsize,9999)
    enddo
    
    do i4=1,ngrid
       r8mxsize(i4)=f0(i4)
    enddo
    
    write(9999,'(/" f0=")') 
    call prtmatcw(nni,mxnmu,r8mxsize,9999)
    isym4nu=1
    isym4mu=1
    call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize,wk3)
    call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
    call putout (nni,mxnmu,r8mxsize,wk0)
    r8mxsize1=r8mxsize
    
    ! call checkd1nu(nni,mxnmu,r8mxsize)

    write(9999,'(/" f0=",i3,1x,a8,1x,a1," / nu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
    call prtmatcw(nni,mxnmu,r8mxsize,9999)
    
    call diff1mu (mxnmu,wk3,wk2)
    call putout (nni,mxnmu,r8mxsize,wk2)
    
    ! call checkd1mu(nni,mxnmu,r8mxsize)
  
    write(9999,'(/" f0=",i3,1x,a8,1x,a1," / mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
    call prtmatcw(nni,mxnmu,r8mxsize,9999)
    
    isym4nu=-1
    isym4mu=1
    call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize1,wk3)
    call diff1mu (mxnmu,wk3,wk2)
    call putout (nni,mxnmu,r8mxsize,wk2)
    
    write(9999,'(/" f0=",i3,1x,a8,1x,a1," / nu,mu derivative" )') iorn(iorb),bond(iorb),gusym(iorb)
    call prtmatcw(nni,mxnmu,r8mxsize,9999)
    
    close(9999)
    
  end subroutine writeDisk4pair

  ! ### writeDisk4dft ### 
  !
  !     Writes orbitals, potenials, Lagrange multipliers (diagonal 
  !     and off-diagonal) and multipole expansion coefficients to a disk 
  !     file in either formatted or unformatted form
  !
  subroutine writeDisk4dft 

    use blas
    use commons
    use data4II
    use dftvxc    
    use discrete
    use elocc
    use memory
    use params
    use printUtils
    use scfshr
    use sharedMemory
    use solver
    use utils

    implicit none
    
    integer (KIND=IPREC) :: iborb,iborb1,ibpot,ibpot1,iorb1,ipc,&
         ngorb,ngorb1,ngpot,ngpot1,isiorb1
    
    real (PREC) :: const13,ocdown,ocup,tmpf
    real (PREC),dimension(:), pointer :: psi,excp,e,f0,f1,f2,f3,f4,fock1,fock2,&
         wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13,work

    parameter (const13=1.0_PREC/3.0_PREC)

    e=>supplptr(i4b(4):)
    psi=>orbptr
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)

    work=>scratchptr    
    fock1   => work(           1: 1*mxsize8)
    fock2   => work( 1*mxsize8+1: 2*mxsize8)

    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk4 =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    wk5 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    wk6 =>scratchptr( 6*mxsize8+1: 7*mxsize8)
    wk7 =>scratchptr( 7*mxsize8+1: 8*mxsize8)            
    wk8 =>scratchptr( 8*mxsize8+1: 9*mxsize8)
    wk9 =>scratchptr( 9*mxsize8+1:10*mxsize8)
    wk10=>scratchptr(10*mxsize8+1:11*mxsize8)
    wk11=>scratchptr(11*mxsize8+1:12*mxsize8)
    wk12=>scratchptr(12*mxsize8+1:13*mxsize8)
    wk13=>scratchptr(13*mxsize8+1:14*mxsize8)
    
    !   contributions from one-electron terms
    iorb=1
    iborb=i1b(iorb)
    ngorb=i1si(iorb)
    ibpot=i2b(iorb)
    ngpot=i2si(iorb)
    
    call zeroArray(mxsize,fock1)
    call zeroArray(mxsize,fock2)
    
    call zeroArray(mxsize,excp(length3-mxsize:))
    call zeroArray(mxsize,wk2)
    
    ! contributions from the local exchange approximation
    
    if (nel.gt.1) then
       ! LDA
       if (idftex.eq.1) then
          do iorb1=1,norb
             if (inhyd(iorb1).eq.1) cycle
             iborb1 =i1b (iorb1)
             isiorb1=i1si(iorb1)
             call exocc (iorb1,ocup,ocdown)
             call prodas (isiorb1,ocup,  psi(iborb1:),psi(iborb1:),fock2)
             call prodas (isiorb1,ocdown,psi(iborb1:),psi(iborb1:),wk2)
          enddo
          
          do i=1,mxsize
             fock2(i)=(fock2(i))**const13
             wk2(i)=(wk2(i))**const13
          enddo
          
          call add (mxsize,fock2,wk2)
          
          ! multiply the local exchange potential by f4 to make it
          ! commensurate with the coulomb potential
          tmpf=fdftpot(alphaf)
          do i=1,mxsize
             wk2(i)=tmpf*wk2(i)*f4(i)
          enddo
          ! B88
       elseif (idftex.eq.2) then
          call fbe88(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
       elseif (idftex.gt.2) then
          write(*,'(/"Warning: unsupported exchange potential")')
          stop 'fockDFT'
       endif
       
       !      store the local exchange potential in a separate array
       call dcopy (mxsize,wk2,ione,fock1,ione)
       
       call zeroArray(mxsize,wk2)
       
       !      add contributions from correlations potentials
       if (idftcorr.eq.1) then
          ! LYP
          call flypcs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
          ! VWN
       elseif (idftcorr.eq.2) then
          call fvwncs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
          
       elseif (idftcorr.gt.2) then
          write(*,'(/"Warning: unsupported correlation potential")')
          stop 'fockDFT'
       endif
       
       ! add correlation contributions to the local exchange ones
       call add (mxsize,wk2,fock1)
    endif
    
    ! add local exchange and correlations contributions ones to one-electron contribution 
    !  call dcopy (mxsize,wk2,ione,fock1,ione) 
    
    !  call dcopy (mxsize,fock2,ione,wk0,ione) 
    
    call zeroArray(mxsize,wk2)
    call zeroArray(mxsize,fock2)
    
    if (nel.gt.1) then
       ! add contributions from Coulomb and off-diagonal Lagrange multipliers.
       do iorb1=1,norb
          if (inhyd(iorb1).eq.1) cycle
          iborb1=i1b(iorb1)
          ngorb1=i1si(iorb1)
          ibpot1=i2b(iorb1)
          ngpot1=i2si(iorb1)
          
          ipc=iorb1+norb*(iorb-1)
          
          !         in the local exchange approximation the Coulomb potential also
          !         includes the contribution from the orbital
          
          call daxpy (ngpot1,occ(iorb1),excp(ibpot1:),ione,wk2,ione)
          
          if (iorb.ne.iorb1.and.abs(ee(iorb1,iorb))>epsilon(zero)) then
             do i=1,ngorb1
                fock2(i)=fock2(i)+ee(iorb1,iorb)*f4(i)*psi(iborb1+i-1)
             enddo
          endif
       enddo
     
       ! store the local exchange potential in exch array as its not
       ! used in HFS/DFT (to be used by EaDFT and EabDFT)
       ! at the end of excp array (important when scmc is on)
       call dcopy (mxsize,fock1,ione,excp(length3-mxsize:),ione)
       
       ! add the coulomb potential to the local exchange one
       call add (mxsize,wk2,fock1)
       
       ! multiply coulomb/exchange potentials and off-diagonal Lagrange
       ! multipliers by f2
     
       call prod (mxsize,f2,fock1)
       call prod (mxsize,f2,fock2)
       ! nel.gt.1 
    endif
    
    do i=1,mxsize
       fock2(i)=fock2(i)+fock1(i)
    enddo
    
    open(9999,file='out4dft.dat', status='unknown',form='formatted')
    
    write(9999,'(80a1)') title
    write(9999,'(" nnu nmu ")') 
    write(9999,formint) nni,nmu(1)
    write(9999,'(" r  rinf")') 
    write(9999,formfp64) r,rinf
    write(9999,'(" z1 z2 ")') 
    write(9999,formfp64) z1,z2
    write(9999,'(" norb nel  ")')
    write(9999,formint) norb,nel
    
    write(9999,'(/" DFT potential" )')
    call prtmatcw(nni,mxnmu,fock2,9999)
    
    close(9999)
    
  end subroutine writeDisk4dft

  ! ### writeDisk4lxc ### 
  !
  !     Writes orbitals, potenials, Lagrange multipliers (diagonal 
  !     and off-diagonal) and multipole expansion coefficients to a disk 
  !     file in either formatted or unformatted form
  !
  subroutine writeDisk4lxc 

    use blas
    use commons
    use data4II
    use dftvxc    
    use discrete
    use elocc
    use memory
    use params
    use printUtils
    use scfshr
    use sharedMemory
    use solver
    use utils

    implicit none
    
    integer (KIND=IPREC) :: iborb,iborb1,ibpot,ibpot1,iorb1,ipc,&
         ngorb,ngorb1,ngpot,ngpot1,isiorb1
    
    real (PREC) :: tmpf,zshift
    real (PREC),dimension(:), pointer :: fock1,fock2,work

    work=>scratchptr    
    fock1   => work(           1: 1*mxsize8)
    fock2   => work( 1*mxsize8+1: 2*mxsize8)
    
    !   contributions from one-electron terms
    iorb=1
    iborb=i1b(iorb)
    ngorb=i1si(iorb)
    ibpot=i2b(iorb)
    ngpot=i2si(iorb)
    
    open(9999,file='out4dft-1.dat', status='unknown',form='formatted')
    
    write(9999,'("###",80a1)') title
    write(9999,'("### nnu nmu ")') 
    write(9999,formint) nni,nmu(1)
    write(9999,'("### r  rinf")') 
    write(9999,formfp64) r,rinf
    write(9999,'("### z1 z2 ")') 
    write(9999,formfp64) z1,z2
    write(9999,'("### norb nel  ")')
    write(9999,formint) norb,nel
    
    write(9999,'(/"### DFT potential" )')
    !call prtmatcw(nni,mxnmu,fock1,9999)
    zshift=r/two

    i=1
    do j=nni-1,1,-1
       write(9999,'(f10.4,f20.10)') (r/two*vxi(i)*veta(j)+zshift),-fock1((i-1)*nni+j)
    enddo
    
    j=1
    do i=2,mxnmu-4
       write(9999,'(f10.4,f20.10)') (r/two*vxi(i)*veta(j)+zshift),-fock1((i-1)*nni+j)
    enddo
    !endif
    close(9999)

    open(9999,file='out4dft-2.dat', status='unknown',form='formatted')
    
    write(9999,'(/"### -1/r" )')
    zshift=r/two

    ! j=nni
    ! do i=mxnmu-4,2,-1
    !    write(9999,'(f10.4,f20.10)') (r/two*vxi(i)*veta(j)+zshift),-one/(r/two*vxi(i)*veta(j)+zshift)                 
    ! enddo
       
    i=1
    do j=nni-1,1,-1
       write(9999,'(f10.4,f20.10)') (r/two*vxi(i)*veta(j)+zshift),-one/(r/two*vxi(i)*veta(j)+zshift)                 
    enddo
    
    j=1
    do i=2,mxnmu-4
       write(9999,'(f10.4,f20.10)') (r/two*vxi(i)*veta(j)+zshift),-one/(r/two*vxi(i)*veta(j)+zshift)                 
    enddo
    !endif
    close(9999)
    
  end subroutine writeDisk4lxc

  
  ! ### writeDisk4kinpot ### 
  !
  !     Calculates and writes to disk the Pauli and von Weizsaecker kinetic potentials.
  !     In case of HF, OC, HLi systems z-coordinate is shifted.
  subroutine writeDisk4kinpot
    use blas
    use commons
    use dftvxc    
    use discrete
    use elocc
    use memory
    use nabla    
    use params
    use scfshr
    use sharedMemory
    use solver
    use utils

    implicit none

    integer (KIND=IPREC),parameter :: iunit=9999
    integer (KIND=IPREC) :: i,ii,iborb1,ibpot1,iorb,iorb1,igp,imu,inu,j,&
         ngorb,ngorb1,ngpot,ngpot1,np

    real (PREC) :: xnorm,zshift
    real (PREC),dimension(:), pointer :: psi,excp,e,f0,f1,f2,f3,f4,wgt1,wgt2,&
         wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    e=>supplptr(i4b(4):)
    psi=>orbptr
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)
    
    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk4 =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    wk5 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    wk6 =>scratchptr( 6*mxsize8+1: 7*mxsize8)
    wk7 =>scratchptr( 7*mxsize8+1: 8*mxsize8)            
    wk8 =>scratchptr( 8*mxsize8+1: 9*mxsize8)
    wk9 =>scratchptr( 9*mxsize8+1:10*mxsize8)
    wk10=>scratchptr(10*mxsize8+1:11*mxsize8)
    wk11=>scratchptr(11*mxsize8+1:12*mxsize8)
    wk12=>scratchptr(12*mxsize8+1:13*mxsize8)
    wk13=>scratchptr(13*mxsize8+1:14*mxsize8)

    call zeroArray(mxsize,wk0)        
    call zeroArray(mxsize,wk2)    

    ! total density: wk0 = \rho = \sum_i q_i \rho_i
    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call prodas (mxsize,occ(iorb1),psi(iborb1:),psi(iborb1:),wk0)
    enddo

    ! tau: wk2 = 1/2 \sum_i q_i \nabla psi_i \nabla psi_i     
    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call nfng(psi(iborb1:),psi(iborb1:),wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
       call dscal (mxsize,half*occ(iorb1),wk10,ione)
       call add(mxsize,wk10,wk2)
    enddo

#ifdef PRINT
    ! print=140: printResults: 1/2\sum_i q_i <|\nabla \phi_i|^2>
    if (iprint(140).ne.0) then
       call dcopy (mxsize,wk2,ione,wk8,ione)
       call prod (mxsize,f4,wk8)
       xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
       write(*,'(/3x,"  1/2\sum_i q_i <|\nabla \phi_i|^2>")')
       write(*,'(3x,"  via nfng                    ",1pe25.16)') xnorm       

    ! tau: w13 
       call zeroArray(mxsize,wk10)
       call tau(psi,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)
       call dscal (mxsize,half,wk13,ione)
       
       call dcopy (mxsize,wk13,ione,wk8,ione)
       call prod (mxsize,f4,wk8)
       xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
       write(*,'(3x,"  via tau                     ",1pe25.16)') xnorm
    endif
#endif    

    ! w12 = \nabla^2 \rho
    !call n2f(wk0,wk5,wk6,wk7,wk12)
    !call n2rho(psi,wk4,wk5,wk6,wk7,wk12)
    call laplace(psi,e,wk5,wk6,wk7,wk8,wk12)

    ! w12 = \nabla^2 \rho + 4 \tau
    call dcopy (mxsize,wk13,ione,wk8,ione)
    call dscal(mxsize,four,wk8,ione)
    call add (mxsize,wk8,wk12)

#ifdef PRINT
    ! print=141: printResults: \sum_i (q_i \laplace psi_i ce) + 4 \tau  
    if (iprint(141).ne.0) then
       call dcopy (mxsize,wk12,ione,wk8,ione)
       call prod (mxsize,f4,wk8)
       xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
       write(*,'(5x,"\sum_i (q_i \laplace psi_i) + 4 \tau  ")')
       write(*,'(33x,1pe25.16)') xnorm       
    endif
#endif
    
    call zeroArray(mxsize,wk10)
    
    ! wk4 == Pauli kinetic potential
    ! wk5 == von Weizsaecker kinetic potential

    ! wk10 = |\nabla \rho |^2
    call nfng(wk0,wk0,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)

    ! wk3 = sum_i (E_{HOMO} -E_i) \rho_i
    call zeroArray(mxsize,wk3)
    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call prod2 (mxsize,psi(iborb1:),psi(iborb1:),wk4)
       call dscal (mxsize, occ(iorb1),wk4,ione)
       call dscal (mxsize, ee(ione,ione)-ee(iorb1,iorb1),wk4,ione)
       call add (mxsize,wk4,wk3)
    enddo

    !call n2rho(psi,wk4,wk5,wk6,wk7,wk8)
    call n2f(wk0,wk4,wk5,wk6,wk8)

    ! call dcopy (mxsize,wk8,ione,wk7,ione)
    ! ! correct the values of Laplacian along the internuclear axis
    ! call setBVgen(ione,ione,wk7)

    ! call dcopy (mxsize,wk12,ione,wk8,ione)
    ! ! correct the values of Laplacian along the internuclear axis
    ! call setBVgen(ione,ione,wk8)
    
    ! wk6 = tau^W = |\nabla \rho|^2/(8\rho)
    !                    wk10          wk0
    ! wk3 = sum_i (E_{HOMO} -E_i) \rho_i / \rho
    do inu=1,nni
       do imu=1,mxnmu
          igp=(imu-1)*nni + inu
          if (abs(wk0(igp))>precis) then
             wk6(igp)=wk10(igp)/(eight*wk0(igp))
             wk3(igp)=wk3(igp)/wk0(igp)
             wk4(igp)=( wk13(igp) - wk6(igp)) /wk0(igp) + wk3(igp)
             wk5(igp)=wk6(igp)/wk0(igp)-wk12(igp)/(four*wk0(igp))
          else
             wk3(igp)=zero
             wk4(igp)=zero
             wk5(igp)=zero
             wk6(igp)=zero             
          endif
       enddo
    enddo

    open(iunit,file='out4kinpot.dat', status='unknown',form='formatted')
    !zshift=zero
    write(iunit,1000)
    ! v_k=V^P + v^W
1000 format (&
          "      #z+R/2    ",&
          "       rho      ",&
          "  nabla^2 rho   ",&
          "       tau      ",&
          "      tau^W     ",&
          "       v^P      ",&
          "       v^W      ",&
          "     v^P+v^W    ",&
          " v^P (2nd term) ")

    call dcopy (mxsize,wk12,ione,wk8,ione)
    ! correct the values of Laplacian along the internuclear axis
    call setBVgen(ione,ione,wk8)

    if (abs(z1-one)<epsilon(zero).and.abs(z2-nine)<epsilon(zero)) then    
    !F (A) 0.04618936
    !H (A) -0.87071064

       zshift=-0.87071064_PREC/bohr2ang
       print *,"r(H) r(F): ", zshift, zshift+r
    endif
    
    if (abs(z1-eight)<epsilon(zero).and.abs(z2-six)<epsilon(zero)) then 
       ! O (A) -0.48351647
       ! C (A)  0.64448353 
       zshift=-0.48351647_PREC/bohr2ang
       print *,"r(O) r(C): ", zshift, zshift+r
    endif

    if (abs(z1-one)<epsilon(zero).and.abs(z2-three)<epsilon(zero)) then        
       ! H (A) -2.64
       ! Li (A) 
       zshift=-1.133_PREC
       print *,"r(H) r(Li): ", zshift, zshift+r
    endif
    
    np=0
    
    if (abs(z2)<epsilon(zero)) then
       !print *, "single nucleus at r=0"
       i=1
       do j=nni,1,-1
          np=np+1
          write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two),&       
               wk0((i-1)*nni+j),wk8((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
               wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
       enddo

    else
       j=nni
       do i=mxnmu-4,2,-1
          !do i=2,mxnmu-4
          np=np+1
          write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two+zshift),&       
               wk0((i-1)*nni+j),wk8((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
               wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
       enddo
       
       
       i=1
       do j=nni,1,-1
          np=np+1
          write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two+zshift),&       
               wk0((i-1)*nni+j),wk8((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
               wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
       enddo
       
       j=1
       do i=2,mxnmu-4
          np=np+1
          write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two+zshift),&       
               wk0((i-1)*nni+j),wk8((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
               wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
       enddo
    endif
    close(iunit)

  end subroutine writeDisk4kinpot

  ! ### writeDisk4kinpot ### 
  !
  !     Calculates and writes to disk the Pauli and von Weizsaecker kinetic potentials.
  !
  subroutine writeDisk4kinpot2 (psi,pot,excp,e,f0,f1,f2,f3,f4,wgt1,wgt2,&
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)

    use blas
    use commons
    use dftvxc    
    use discrete
    use elocc
    use memory
    use nabla    
    use params
    use scfshr
    use solver
    use utils

    implicit none

    integer (KIND=IPREC),parameter :: iunit=9999
    integer (KIND=IPREC) :: i,ii,iborb1,ibpot1,iorb,iorb1,igp,imu,inu,j,&
         ngorb,ngorb1,ngpot,ngpot1,np

    real (PREC) :: xnorm,zshift
    real (PREC),dimension(*) :: psi,pot,excp,e,f0,f1,f2,f3,f4,wgt1,wgt2,wk0,wk1,&
         wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    call zeroArray(mxsize,wk0)        
    call zeroArray(mxsize,wk2)    

    ! wk0 = \rho = \sum_i \rho_i
    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call prodas (mxsize,occ(iorb1),psi(iborb1),psi(iborb1),wk0)
    enddo

    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call nfng(psi(iborb1),psi(iborb1),wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
       call dscal (mxsize,half*occ(iorb1),wk10,ione)
       call add(mxsize,wk10,wk2)
    enddo
    
    call dcopy (mxsize,wk2,ione,wk8,ione)
    call prod (mxsize,f4,wk8)
    xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
    write(*,'(" 1/2\sum <|\nabla \phi|^2> = ",1pe16.8)') xnorm

    ! w13 = tau
    call zeroArray(mxsize,wk10)
    call tau(psi,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)
    call dscal (mxsize,half,wk13,ione)

    call dcopy (mxsize,wk13,ione,wk8,ione)
    call prod (mxsize,f4,wk8)
    xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
    write(*,'("    call tau:       <\tau> = ",1pe16.8)') xnorm

    ! w12 = \nabla^2 \rho
    !call n2f(wk0,wk5,wk6,wk7,wk12)
    !call n2rho(psi,wk4,wk5,wk6,wk7,wk12)
    call laplace(psi,e,wk5,wk6,wk7,wk8,wk12)

    call dcopy (mxsize,wk13,ione,wk8,ione)
    call dscal(mxsize,four,wk8,ione)
    call add (mxsize,wk8,wk12)
    !call dscal(mxsize,two,wk12,ione)

    ! 
    ! do i=1,mxnmu
    !    ii=(i-1)*nni
    !    do j=1,nni
    !       wk12(ii+j)=four*pii*(r/two*vxi(i)*veta(j)+r/two)*wk12(ii+j)
    !    enddo
    ! enddo

    call dcopy (mxsize,wk12,ione,wk8,ione)
    !call prod (mxsize,wk0,wk8)
    call prod (mxsize,f4,wk8)
    xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
    write(*,'("wgt2:  \nabla^2 \rho (tau) = ",1pe16.8)') xnorm


    !call n2f(wk0,wk5,wk6,wk7,wk12)
    ! call dcopy (mxsize,wk12,ione,wk8,ione)
    ! call prod (mxsize,wk0,wk8)
    ! call prod (mxsize,f4,wk8)
    ! xnorm=ddot (mxsize,wgt2,ione,wk8,ione)
    ! write(*,'("wgt2:  \nabla^2 \rho (n2f) = ",1pe16.8)') xnorm

    call zeroArray(mxsize,wk10)
    
    ! call prod (mxsize,f4,wk2)
    ! xnorm=ddot (mxsize,wgt2,ione,wk2,ione)
    ! print *,"tau=",xnorm
    
    ! wk4 == Pauli kinetic potential
    ! wk5 == von Weizsaecker kinetic potential

    ! wk10 = |\nabla \rho |^2
    call nfng(wk0,wk0,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)


    ! wk3 = sum_i (E_{HOMO} -E_i) \rho_i
    call zeroArray(mxsize,wk3)

    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call prod2 (mxsize,psi(iborb1),psi(iborb1),wk4)
       call dscal (mxsize, occ(iorb1),wk4,ione)
       call dscal (mxsize, ee(ione,ione)-ee(iorb1,iorb1),wk4,ione)
       call add (mxsize,wk4,wk3)
    enddo

    call n2rho(psi,wk4,wk5,wk6,wk7,wk8)
    
    ! wk6 = tau^W = |\nabla \rho|^2/(8\rho)
    !                    wk10          wk0
    do inu=1,nni
       do imu=1,mxnmu
          igp=(imu-1)*nni + inu
          if (abs(wk0(igp))>precis) then
             wk6(igp)=wk10(igp)/(eight*wk0(igp))
             wk3(igp)=wk3(igp)/wk0(igp)
             !wk4(igp)=( wk13(igp) - wk6(igp)) /wk0(igp) + wk3(igp)/wk0(igp)
             wk4(igp)=( wk13(igp) - wk6(igp)) /wk0(igp) + wk3(igp)

             wk5(igp)=wk6(igp)/wk0(igp)-wk12(igp)/(four*wk0(igp))
             !print *,igp,inu,imu,wk0(igp),wk4(igp),wk3(igp)
          else
             wk3(igp)=zero
             wk4(igp)=zero
             wk5(igp)=zero
             wk6(igp)=zero             
          endif
       enddo
    enddo

    open(iunit,file='out4kinpot.dat', status='unknown',form='formatted')
    
    write(iunit,1000)
1000 format("     #z+R/2",11x,"rho",6x,4x,"nabla^2 rho",8x,"tau",&
          13x,"tau^W",11x,"v^P",14x,"v^W",13x,"v_k",13x,"allee")
    

    call dcopy (mxsize,wk12,ione,wk8,ione)
    ! correct the values of Laplacian along the internuclear axis
    call setBVgen(ione,ione,wk8)

    ! j=nni
    ! do i=2,mxnmu-4
    !    np=np+1
    !    write(iunit,'(1p12e16.8)') -(r/two*vxi(i)*veta(j)+r/two),&       
    !         wk12((i-1)*nni+j),(wk12((i-1)*nni+j)-wk8((i-1)*nni+j))
    ! enddo

    ! write(iunit,'("tau")')    
    
    ! call dcopy (mxsize,wk13,ione,wk8,ione)
    ! ! correct the values of Laplacian along the internuclear axis
    ! call setBVgen(ione,ione,wk8)

    ! do i=2,mxnmu-4
    !    np=np+1
    !    write(iunit,'(1p12e16.8)') -(r/two*vxi(i)*veta(j)+r/two),&       
    !         wk13((i-1)*nni+j),(wk13((i-1)*nni+j)-wk8((i-1)*nni+j))
    ! enddo


    ! write(iunit,'("tau^W")')
    
    ! call dcopy (mxsize,wk6,ione,wk8,ione)
    ! ! correct the values of Laplacian along the internuclear axis
    ! call setBVgen(ione,ione,wk8)

    ! do i=2,mxnmu-4
    !    np=np+1
    !    write(iunit,'(1p12e16.8)') -(r/two*vxi(i)*veta(j)+r/two),&       
    !         wk6((i-1)*nni+j),(wk6((i-1)*nni+j)-wk8((i-1)*nni+j))
    ! enddo


    ! write(iunit,'("v^P")')
    
    ! call dcopy (mxsize,wk4,ione,wk8,ione)
    ! ! correct the values of Laplacian along the internuclear axis
    ! call setBVgen(ione,ione,wk8)

    ! do i=2,mxnmu-4
    !    np=np+1
    !    write(iunit,'(1p12e16.8)') -(r/two*vxi(i)*veta(j)+r/two),&       
    !         wk4((i-1)*nni+j),(wk4((i-1)*nni+j)-wk8((i-1)*nni+j))
    ! enddo


    ! write(iunit,'("v^W")')
    
    ! call dcopy (mxsize,wk5,ione,wk8,ione)
    ! ! correct the values of Laplacian along the internuclear axis
    ! call setBVgen(ione,ione,wk8)

    ! do i=2,mxnmu-4
    !    np=np+1
    !    write(iunit,'(1p12e16.8)') -(r/two*vxi(i)*veta(j)+r/two),&       
    !         wk5((i-1)*nni+j),(wk5((i-1)*nni+j)-wk8((i-1)*nni+j))
    ! enddo
    
    ! return

    
    !F (A) 0.04618936
    !H (A) -0.87071064

    !zshift=-0.87071064/bohr2ang
    !print *,"r(H) r(F): ", zshift, zshift+r

    ! O (A) -0.48351647
    ! C (A)  0.64448353 
    zshift=-0.48351647/bohr2ang
    print *,"r(O) r(C): ", zshift, zshift+r
    
    np=0

    j=nni
    do i=mxnmu-4,2,-1
       !do i=2,mxnmu-4
       np=np+1
       write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two+zshift),&       
            wk0((i-1)*nni+j),wk12((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
            wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
    enddo

    
    i=1
    do j=nni,1,-1
       np=np+1
       write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two+zshift),&       
            wk0((i-1)*nni+j),wk8((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
            wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
    enddo

    j=1
    do i=2,mxnmu-4
       np=np+1
       write(iunit,'(1p12e16.8)') (r/two*vxi(i)*veta(j)+r/two+zshift),&       
            wk0((i-1)*nni+j),wk8((i-1)*nni+j),wk13((i-1)*nni+j),&                                    
            wk6((i-1)*nni+j),wk4((i-1)*nni+j),wk5((i-1)*nni+j),wk4((i-1)*nni+j)+wk5((i-1)*nni+j),wk3((i-1)*nni+j)
    enddo

    close(iunit)

  end subroutine writeDisk4kinpot2
  
end module diskInterfaceMisc
  
