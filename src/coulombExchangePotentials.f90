! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module coulombExchangePotentials
  implicit none
contains
  ! ### initCoulomb ###
  !
  !     Initializes Coulomb potentials using the Thomas-Fermi model or SAP
  !     and exchange potentials via the local exchange or 1/r approximation
  !
  subroutine initCoulomb
    use params
    use sharedMemory
    use discrete
    use memory
    use commons
    use blas
    use utils
    
    implicit none

    integer (KIND=IPREC) :: i3beg,igp,imu,in,inioff,iorb,iorb1,iorb2,iorb2t,irec,ishift,k
    real (PREC), parameter :: slim=1.0_PREC
    real (PREC), parameter :: zc1=1.0_PREC
    real (PREC), parameter :: zc2=1.0_PREC
    real (PREC) :: ch1,ch2,crt1,crt2,ez1,ez2,ra1,ra2,vetat,vxit
    real (PREC), dimension(:), pointer :: excp,excp1,f4

    excp=>exchptr
    f4=>supplptr(i4b(9):)
    
    if (OED.or.initFuncsOED) return

    print *,'... initializing Coulomb potentials (pottf) ...'

    if (.not.linitFuncsNoexch) then
       ! Initialization of Coulomb potentials
       ! loop over orbitals
       do iorb=1,norb
          if (inhyd(iorb).eq.0) cycle
          ishift=i1b(iorb)-1
          ez1=eza1(iorb)
          ez2=eza2(iorb)
          ch1=-1.00_PREC
          ch2=-1.00_PREC

          if (ez1.lt.precis) ch1=-2.00_PREC
          if (ez2.lt.precis) ch2=-2.00_PREC

          crt1=abs(co1(iorb))
          crt2=abs(co2(iorb))
          do imu=1,mxnmu
             inioff=(imu-1)*nni
             vxit=vxi(imu)
             do in=1,nni
                igp=ishift+inioff+in
                vetat=veta(in)
                if (imu.eq.1.and.in.eq.1) then
                   excp(igp)=0.0_PREC
                else
                   excp(igp)=pottf(r, vetat,vxit,zc1,ch1,slim)*crt1+pottf(r,(-vetat),vxit,zc2,ch2,slim)*crt2
                   !if (mod(igp,2000)==0) write(*,'(2i5,7e12.4)') imu,in,r,crt1,crt2,ez1,ez2,excp(igp)
                endif
             enddo
          enddo
          excp1=>exchptr(i1b(iorb):)
          call prod(mxsize,f4,excp1)
       enddo
    endif

  end subroutine initCoulomb

  ! This routine initializes exchange potentials via the local exchange or 1/r approximation
  subroutine initExchange
    use params
    use sharedMemory
    use discrete
    use memory
    use commons
    use blas
    use utils

    implicit none

    integer (KIND=IPREC) :: i3beg,igp,imu,in,inioff,iorb,iorb1,iorb2,iorb2t,irec,ishift,k
    real (PREC) :: ch1,ch2,crt1,crt2,ez1,ez2,ra1,ra2,slim,vetat,vxit,zc1,zc2
    !real (PREC), dimension(*) :: psi,excp,f4

    real (PREC), dimension(:), pointer :: excp,excp1,excp2,f4,psi

    excp=>exchptr
    f4=>supplptr(i4b(9):)
    psi=>orbptr

    if (OED.or.initFuncsOED) return

    ! Initialization of exchange potentials

    if (LXC.or.DFT.or.HFS.or.SCMC) then
       print *,'... initializing Slater exchange potential ...'
       call slaterPot
       return
    else
       print *,'... initializing exchange potentials ...'
       
       ! Initialization of exchange potentials that are all kept in memory
       
       do iorb1=1,norb
         do iorb2=iorb1,norb
             k=iorb1+iorb2*(iorb2-1)/2
             if (iorb1.eq.iorb2.and.ll(iorb2).eq.0) cycle
             ishift=i3b(k)-1

             ! loop over grid points
             do imu=1,mxnmu
                if (imu.eq.1) cycle
                inioff=(imu-1)*nni
                vxit=vxi(imu)
                do in=1,nni
                   if (in.eq.1.or.in.eq.nni) cycle
                   igp=ishift+inioff+in
                   vetat=veta(in)
                   ra1=(r/2.0_PREC)*(vxit+vetat)
                   ra2=(r/2.0_PREC)*(vxit-vetat)
                   if (ra1.gt.1.d-8) then
                      excp(igp)=co1(iorb1)/ra1
                   else
                      excp(igp)=0.0_PREC
                   endif
                   if (ra2.gt.1.d-8) then
                      excp(igp)=excp(igp)+co2(iorb1)/ra2
                   endif
                enddo
             enddo

             excp1=>exchptr((i3b(k)):)
             excp2=>exchptr((i3b(k))+mxsize:)
             call prod(mxsize,f4,excp1)
             if (iorb1.ne.iorb2) then
                if (ll(iorb1).ne.0.and.ll(iorb2).ne.0) then
                   ishift=ishift+mxsize
                   call dcopy(mxsize,excp1,ione,excp2,ione)
                endif
             endif
          enddo
       enddo
    endif

  end subroutine initExchange
  

  subroutine initCoulombSAP 
    !USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE    
    use params
    use sharedMemory
    use discrete
    use memory
    use commons
    use blas
    use sapPotential
    use utils

    implicit none

    integer (KIND=IPREC) :: i3beg,igp,imu,in,inioff,iorb,iorb1,iorb2,iorb2t,irec,ishift,iz1,iz2,k
    real (PREC) :: ra1,ra2,r1t,r2t,vetat,vxit
    real (PREC) :: infinity
    real (PREC), dimension(:), pointer :: excp,excp1,f4

    excp=>exchptr
    f4=>supplptr(i4b(9):)
    
    infinity=huge(one)
    
    if (OED.or.initFuncsOED) return

    print *,'... initializing Coulomb potentials (effective_coulomb_charge) ...'
    
    if (.not.linitFuncsNoexch) then
       ! Initialization of Coulomb potentials
       iz1=nint(z1)
       iz2=nint(z2)

       do iorb=1,norb

          ishift=i1b(iorb)-1

          !  loop over grid points
          do imu=1,mxnmu
             inioff=(imu-1)*nni
             vxit=vxi(imu)
             do in=1,nni
                igp=ishift+inioff+in
                vetat=veta(in)
                r1t=(r/2.0_PREC)*(vxi(imu)+veta(in))
                r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))
                
                excp(igp)=z1-effective_coulomb_charge(iz1,r1t)*co1lda+&
                      z2-effective_coulomb_charge(iz2,r2t)*co2lda
                if (abs(excp(igp))>infinity) then
                   write(*,'("initCoulombSAP",i2,a8,3i5,2e12.4)') &
                        iorn(iorb),bond(iorb),imu,in,iz1,r1t,effective_coulomb_charge(iz1,r1t)
                   write(*,'("initCoulombSAP",i2,a8,3i5,2e12.4)') &
                        iorn(iorb),bond(iorb),imu,in,iz2,r2t,effective_coulomb_charge(iz2,r2t)
                   excp(igp)=excp(igp-1)
                endif
             enddo
          enddo
          excp1=>excp(i1b(iorb):)
          call prod(mxsize,f4,excp1)
       enddo
    endif

  end subroutine initCoulombSAP

  ! ### pottf ###
  !
  !     Calculates a Thomas-Fermi potential (a modified version of dentfa by Desclaux)
  !
  function pottf(r,ek,qk,dz,ch,slim)
    use params

    implicit none
    real (PREC) :: pottf
    real (PREC) :: ch,dz,ek,qk,r,slim,t,w

    pottf=0.0_PREC
    if ((dz+ch).lt.1.e-04_PREC)  return
    if (abs(qk+ek).lt.1.e-10_PREC) return

    w=r*(qk+ek)/2.0_PREC*(dz+ch)**(1.0_PREC/3.0_PREC)
    w=sqrt(w/0.885340_PREC)
    t=w*(0.6011200_PREC*w+1.8106100_PREC)+1.0_PREC
    w=w*(w*(w*(w*(0.0479300_PREC*w+0.2146500_PREC)+0.7711200_PREC)+1.3951500_PREC)+1.8106100_PREC)+1.0_PREC
    pottf=slim*(1.0_PREC-(t/w)*(t/w))*2.0_PREC/(r*(qk+ek))

  end function pottf

  ! ### slaterPot ###
  !
  !     Calculates exchange potentials according to the local Slater
  !     approximation
  !
  subroutine slaterPot
    use params
    use blas
    use discrete
    use commons
    use elocc
    use sharedMemory
    use utils

    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb

    real (PREC) :: const13,coo,xa,ocdown,ocup
    real (PREC), dimension(:), pointer :: f4,rhodown,rhoup,wk1,wk2
    
    parameter (const13=1.0_PREC/3.0_PREC)

    if (nel.eq.1) return

    f4=>supplptr(i4b(9):)
    rhodown=>scratchptr(           1: 1*mxsize8)
    rhoup  =>scratchptr( 1*mxsize8+1: 2*mxsize8)    
    wk1    =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2    =>scratchptr( 3*mxsize8+1: 4*mxsize8)

    call zeroArray(mxsize,exchptr(1:))
    call zeroArray(mxsize,rhodown)
    call zeroArray(mxsize,rhoup)

    do iorb=1,norb
       !if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)

       call exocc (iorb,ocup,ocdown)
       call prod2 (mxsize,orbptr(iborb:),orbptr(iborb:),wk1)
       call dscal (mxsize,ocup,wk1,ione)

       call prod2 (mxsize,orbptr(iborb:),orbptr(iborb:),wk2)
       call dscal (mxsize,ocdown,wk2,ione)

       ! store total spin densities
       call add(mxsize,wk1,rhoup)
       call add(mxsize,wk2,rhodown)
    enddo

    xa=-three/two*alphaf*(three/pii)**const13    
    call add (mxsize,rhodown,rhoup)
    
    ! multiply exchange potential by f4 to make it commensurate with
    ! Coulomb potential

    do i=1,mxsize
       exchptr(i)=xa*f4(i)*(rhoup(i))**const13
    enddo

  end subroutine slaterPot
  
  
end module coulombExchangePotentials
