! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2023  Jacek Kobus 

module integrals
  implicit none
contains

  ! ### coulij ###
  !
  function coulij (iorb1,iorb2,psi,pot,wgt2,wk0)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use utils    

    implicit none
    integer (KIND=IPREC) :: iorb1,iorb2
    real (PREC) :: coulij
    real (PREC), dimension(*) :: psi,pot,wgt2,wk0
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    call prod2 (mxsize,psi(i1b(iorb2)),pot(i2b(iorb1)),wk0)
    call prod  (mxsize,psi(i1b(iorb2)),wk0)
    coulij=ddot(mxsize,wgt2,ione,wk0,ione)
  end function coulij

  ! ### oneelii ###
  !
  !     This version of oneelij function is only used in etotal routine
  !     where vk(iorb1) and vn(iorb1) values are needed.
  !
  real (PREC) function oneelii(iorb1)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use inout
    use sharedMemory
    use utils
    
    implicit none
    integer (KIND=IPREC) :: iorb1,isym,nmut,i,ii,j,k,np
    real (PREC) :: w,wkin,wnucl
    real (PREC), dimension(:), pointer :: psi1
    real (PREC), dimension(:), pointer :: psi,e,f0,wgt1,&
              wk0,wk1,wk2,wk3
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    e=>supplptr(i4b(4):)
    f0=>supplptr(i4b(5):)
    psi=>orbptr
    wgt1=>supplptr(i4b(13):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    nmut=i1mu(iorb1)
    isym=isymOrb(iorb1)
    psi1=>psi(i1b(iorb1):)    
    call putin (nni,nmut,isym,psi1,wk3)
    call diffnu (nmut,wk3,wk0,wk1,wk2)
    call putout (nni,nmut,wk1,wk0)

    call diffmu (nmut,wk3,wk2)
    call putout (nni,nmut,wk0,wk2)

    call add (mxsize,wk0,wk1)

    if (mm(iorb1).ne.0) then
       ! nuclear energy for non-sigma orbitals contains contribution from e term (in toten
       ! this term is correctly added to the kinetic energy contribution); e enters the
       ! expression with minus sign which is already incorporated in e
       w=dble(mm(iorb1)*mm(iorb1))
       call prodas (mxsize,w,e,psi1,wk1)
    endif

    ! FIXME (debug pragma?)
    if (idebug(707)==1.or.idebug(807)==1) then
       j=1
       print *,"oneelii before interpolation: i=1,j=1", wk1(j)
       wk1(j)=zero
       do i=2,6
          ii=(i-1)*nni
          wk1(j)=wk1(j)+ wk1(ii+j)*exeven(i-1)
       end do
       wk1(j)=-wk1(j)
       print *,"         after interpolation:         ", wk1(j)
       
       i=1
       j=nni
       print *,"oneelii before interpolation: i=1,j=nni", wk1(j)
       wk1(j)=zero
       do k=1,5
          wk1(j)=wk1(j) + wk1(j-k)*exeven(k)
       end do
       print *,"         after interpolation:           ", wk1(j)
       
       open(999,file='t.dat', status='unknown',form='formatted')
       
       write(999,'("### oneelii: T(nu)+T(mu) ")')
       np=0
       j=nni
       do i=mxnmu-4,2,-1
          !do i=mxnmu,1,-1
          np=np+1
          write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk1((i-1)*nni+j)
       enddo
       
       i=1
       do j=nni,1,-1
          np=np+1
          write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk1((i-1)*nni+j)
       enddo
       
       j=1
       do i=2,mxnmu-4
          np=np+1
          write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk1((i-1)*nni+j)
       enddo
       
       close(999)
    endif
    
    call dcopy (mxsize,wk1,ione,wk2,ione)
    call prod  (mxsize,psi1,wk2)
    wkin=ddot(mxsize,wgt1,ione,wk2,ione)
    vk(iorb1)=wkin

    call dcopy (mxsize,f0,ione,wk0,ione)
    
    call prod2 (mxsize,psi1,wk0,wk2)
    call prod  (mxsize,psi1,wk2)
    wnucl=ddot(mxsize,wgt1,ione,wk2,ione)
    vn(iorb1)=wnucl

    call proda (mxsize,psi1,wk0,wk1)
    call prod  (mxsize,psi1,wk1)

    oneelii=ddot(mxsize,wgt1,ione,wk1,ione)

#ifdef PRINT
! print= 65: oneelii: 1-electron contrbutions    
    if(iprint(65).ne.0) then
       write(*,'(" oneelii: energy contributions")') 
       write(*,'(4x,"<",i2,1x,a5,a1,"| T |",i2,1x,a5,a1,"> =",1Pe23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),wkin
       write(*,'(4x,"<",i2,1x,a5,a1,"| V |",i2,1x,a5,a1,"> =",1Pe23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),wnucl
       write(*,'(4x,"<",i2,1x,a5,a1,"|T+V|",i2,1x,a5,a1,"> =",1Pe23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),wkin+wnucl
    endif
#endif
  end function oneelii


  real (PREC) function oneelii4plot(iorb1,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use inout
    use utils
    
    implicit none
    integer (KIND=IPREC) :: iorb1,isym,nmut,i,ii,j,k,np
    real (PREC) :: t1,w,wkin,wnucl,w1,w2,w4
    real (PREC), dimension(*) :: psi,e,f0,wgt1,wk0,wk1,wk2,wk3

    real (PREC), dimension(9) :: aa1,aa2,a1,a2
    real (PREC), dimension(nni,mxnmu) :: borb

    data aa1/ 3.0_PREC, -32.0_PREC, 168.0_PREC, -672.0_PREC, 0.0_PREC, 672.0_PREC, -168.0_PREC, 32.0_PREC, -3.0_PREC /
    data aa2/ -9.0_PREC,  128.0_PREC, -1008.0_PREC, 8064.0_PREC, -14350.0_PREC, 8064.0_PREC,-1008.0_PREC, 128.0_PREC, -9.0_PREC /

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    nmut=i1mu(iorb1)
    isym=isymOrb(iorb1)

    w1=1.0_PREC/(840.0_PREC)
    w2=1.0_PREC/(5040.0_PREC)
    
    do i=1,mxnmu
       w4=vxi2(i)
       do j=1,nni
          borb(j,i)=w4
       enddo
    enddo
  
    ! calculate derivatives over mu and ni
    call putin (nni,nmut,isym,psi(i1b(iorb1)),wk3)
    call diffnu (nmut,wk3,wk0,wk1,wk2)
    call putout (nni,nmut,wk1,wk0)

    call diffmu (nmut,wk3,wk2)
    call putout (nni,nmut,wk0,wk2)

    call add (mxsize,wk0,wk1)

    ! take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
          if (t1.lt.precis) then
             wk1(ii+j)=0.0_PREC
          else
             wk1(ii+j)=wk1(ii+j)*(four/t1)
          endif

       enddo
    enddo

    
    open(999,file='t.dat', status='unknown',form='formatted')

    write(999,'("### oneelii4plot: T")')
    np=0
    j=nni
    do i=mxnmu-4,2,-1
       np=np+1
       write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk1((i-1)*nni+j)
    enddo

    i=1
    do j=nni,1,-1
       np=np+1
       write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk1((i-1)*nni+j)
    enddo

    ! j=1
    ! do i=2,mxnmu-4
    !    np=np+1
    !    write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk0((i-1)*nni+j)
    ! enddo

    

    close(999)

     !derivatives over mu variable
     do k=1,9
        a2(k)=w2*aa2(k)/(hmu(1)*hmu(1))
        a1(k)=w1*aa1(k)/hmu(1)
     enddo
    
     do i=2,mxnmu
        do k=1,9
           dmu(k,i)=a2(k)-borb(1,i)*a1(k)
           d1mu(k,i)=-a1(k)
        enddo
     enddo

     call putin (nni,nmut,isym,psi(i1b(iorb1)),wk3)
     call diffnu (nmut,wk3,wk0,wk1,wk2)
     call putout (nni,nmut,wk1,wk0)

     call diffmu (nmut,wk3,wk2)
     call putout (nni,nmut,wk0,wk2)

     call add (mxsize,wk0,wk1)

     if (mm(iorb1).ne.0) then
        ! nuclear energy for non-sigma orbitals contains contribution from e term (in toten
        ! this term is correctly added to the kinetic energy contribution); e enters the
        ! expression with minus sign which is already incorporated in e
        w=dble(mm(iorb1)*mm(iorb1))
        call prodas (mxsize,w,e,psi(i1b(iorb1)),wk1)
     endif


    ! take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
          if (t1.lt.precis) then
             wk1(ii+j)=0.0_PREC
          else
             wk1(ii+j)=wk1(ii+j)*(four/t1)
          endif

       enddo
    enddo
     
    
     ! j=1
     ! print *,"oneelii4plot before interpolation: i=1,j=1", wk1(j)
     ! wk1(j)=zero
     ! do i=2,6
     !    ii=(i-1)*nni
     !    wk1(j)=wk1(j)+ wk1(ii+j)*exeven(i-1)
     ! end do
     ! wk1(j)=-wk1(j)
     ! print *,"              after interpolation:         ", wk1(j)
    
     ! i=1
     ! j=nni
     ! print *,"oneelii4plot before interpolation: i=1,j=nni", wk1(j)
     ! wk1(j)=zero
     ! do k=1,5
     !    wk1(j)=wk1(j) + wk1(j-k)*exeven(k)
     ! end do
     ! print *,"              after interpolation:           ", wk1(j)
    
    open(999,file='t.dat', status='old',access='append',form='formatted')

    j=1
    do i=2,mxnmu-4
       np=np+1
       write(999,'(3i5,4e16.8)') np,j,i,r/two*vxi(i)*veta(j),psi(i1b(iorb1)+(i-1)*nni+j-1),wk0((i-1)*nni+j)
    enddo
    close(999)

    call dcopy (mxsize,wk1,ione,wk2,ione)
    call prod  (mxsize,psi(i1b(iorb1)),wk2)
    wkin=ddot(mxsize,wgt1,ione,wk2,ione)
    vk(iorb1)=wkin

    call dcopy (mxsize,f0,ione,wk0,ione)
    
    call prod2 (mxsize,psi(i1b(iorb1)),wk0,wk2)
    call prod  (mxsize,psi(i1b(iorb1)),wk2)
    wnucl=ddot(mxsize,wgt1,ione,wk2,ione)
    vn(iorb1)=wnucl

    call proda (mxsize,psi(i1b(iorb1)),wk0,wk1)
    call prod  (mxsize,psi(i1b(iorb1)),wk1)

    oneelii4plot=ddot(mxsize,wgt1,ione,wk1,ione)

#ifdef PRINT
! print= 64: oneelii4plot: T+V 
    if (iprint(64).ne.0.and.iprint(66).eq.0) then
       write(*,'(4x,"<",i2,1x,a5,a1,"|T+V|",i2,1x,a5,a1,"> =",e23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),oneelii4plot
    endif
#endif

#ifdef PRINT
! print= 65: oneelii4plot: T, V, T+V 
    if(iprint(65).ne.0) then
       write(*,'(" oneelii4plot: energy contributions")') 
       write(*,'(18x,"<",i2,1x,a5,a1,"| T |",i2,1x,a5,a1,"> =",e23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),wkin
       write(*,'(18x,"<",i2,1x,a5,a1,"| V |",i2,1x,a5,a1,"> =",e23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),wnucl
       write(*,'(18x,"<",i2,1x,a5,a1,"|T+V|",i2,1x,a5,a1,"> =",e23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb1),bond(iorb1),gusym(iorb1),wkin+wnucl
    endif
#endif
  end function oneelii4plot
  
  ! ### exchij ###
  !
  function exchij (exchType,iorb1,iorb2,psi,excp,wgt2,wk0)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use utils    

    implicit none
    integer (KIND=IPREC) :: ihc,iorb1,iorb2,kex,exchType
    real (PREC) :: exchij
    real (PREC), dimension(*) :: psi,excp,wgt2,wk0
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    kex=iorb1+norb*(iorb2-1)
          
    if (iorb1.le.iorb2) then
       ihc=iorb1+iorb2*(iorb2-1)/2
    else
       ihc=iorb2+iorb1*(iorb1-1)/2
    endif
    
    if (exchType==0) then
       call prod2 (mxsize,psi(i1b(iorb1)),excp(i3b(ihc)),wk0)
    else
       call prod2 (mxsize,psi(i1b(iorb1)),excp(i3b(ihc)+mxsize),wk0)
    endif
       
    call prod (mxsize,psi(i1b(iorb2)),wk0)
    exchij=ddot(mxsize,wgt2,ione,wk0,ione)
  end function exchij


  ! ### oneelij ###
  !
  ! FIXME    <1|T+V_n|2> ?????? or <2|T+V_n|1> ??????
  !
  real (PREC) function oneelij(iorb1,iorb2,psi,e,f0,wgt1,wk0,wk1,wk2,wk3)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use inout
    use utils
    
    implicit none
    integer (KIND=IPREC) :: iorb1,iorb2,isym
    real (PREC) :: w,wkin,wnucl
    real (PREC), dimension(*) :: psi,e,f0,wgt1,wk0,wk1,wk2,wk3
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    isym=isymOrb(iorb2)

    ! FIXME iorb1->iorb2 below?
    ! calculate derivatives over mu and ni
    call putin (nni,i1mu(iorb1),isym,psi(i1b(iorb1)),wk3)
    call diffnu (i1mu(iorb1),wk3,wk0,wk1,wk2)
    call putout (nni,i1mu(iorb1),wk1,wk0)

    call diffmu (i1mu(iorb1),wk3,wk2)
    call putout (nni,i1mu(iorb1),wk0,wk2)

    ! add contribution from derivatives over mu and ni
    call add (mxsize,wk0,wk1)

    if (mm(iorb1).ne.0) then
       ! nuclear energy for non-sigma orbitals contains contribution from e term (in toten
       ! this term is correctly added to the kinetic energy contribution); e enters the
       ! expression with minus sign which is already incorporated in e
       w=dble(mm(iorb1)*mm(iorb1))
       call prodas (mxsize,w,e,psi(i1b(iorb1)),wk1)
    endif

#ifdef PRINT
! print= 62: oneelij: T 
    if (iprint(62).ne.0) then
       call dcopy (mxsize,wk1,ione,wk2,ione)
       call prod  (mxsize,psi(i1b(iorb2)),wk2)
       wkin=ddot(mxsize,wgt1,ione,wk2,ione)
       !vk(iorb2)=wkin
    endif
#endif
    
    call dcopy (mxsize,f0,ione,wk0,ione)    

#ifdef PRINT
! print= 62: oneelij: V
    if (iprint(62).ne.0) then
       call prod2 (mxsize,psi(i1b(iorb1)),wk0,wk2)
       call prod  (mxsize,psi(i1b(iorb2)),wk2)
       wnucl=ddot(mxsize,wgt1,ione,wk2,ione)
       !vn(iorb1)=wnucl
    endif
#endif
    call proda (mxsize,psi(i1b(iorb1)),wk0,wk1)
    call prod  (mxsize,psi(i1b(iorb2)),wk1)

    oneelij=ddot(mxsize,wgt1,ione,wk1,ione)

#ifdef PRINT
! print= 63: oneelij: T+V
    if(iprint(63).ne.0) then
       write(*,'(" oneelij: energy contributions")') 
       write(*,'(4x,"<",i2,1x,a5,a1,"|T+V|",i2,1x,a5,a1,"> =",1Pe23.16)') iorn(iorb1),bond(iorb1),gusym(iorb1),&
            iorn(iorb2),bond(iorb2),gusym(iorb2),oneelij
    endif
#endif
    
  end function oneelij


  ! ### twoelij ###
  !
  !  <2|V_C|1>
  function twoelij (iorb1,iorb2,psi,pot,excp,wgt2,wk0,wk1,wk2,wk3)
    use params
    use discrete
    use scfshr
    use commons
    use blas
    use inout
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ihc,iorb,iorb1,iorb2,kex

    real (PREC) :: coo,w,woneel,wtwoel,wtwoelCoul,wtwoelExch,twoelij
    real (PREC), dimension(*) :: psi,pot,excp,wgt2,wk0,wk1,wk2,wk3
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    ! add contributions from Coulomb (wk0) and exchange potentials (wk1)

    ! coo=occ(iorb)
    ! if (iorbRef.eq.iorb) coo=coo-one
    ! call daxpy (mxsize,occ(iorb),pot(i2b(iorb)),ione,wk0,ione)
    
    ! call prod  (mxsize,psi(i1b(iorbRef)),wk0)


    ! do iorb=1,norb
    !    coo=occ(iorb)
    !    if (iorbRef.eq.iorb) coo=coo-one
    !    call dcopy (mxsize,pot(i2b(iorb)),ione,wk0,ione)
    !    call prod  (mxsize,psi(i1b(iorbRef)),wk0)
    !    if (coo/=one) then
    !       call dscal (mxsize,coo,wk0,ione)
    !    endif
    ! enddo


    ! The code below could be streamlined a bit by first calculating the Coulomb potential
    ! contribution and then the exchange one. However, if this is done, some open-shell
    ! cases fail to converge, e.g. Li one.

    call zeroArray (mxsize,wk0)
    call zeroArray (mxsize,wk1)
    call zeroArray (mxsize,wk2)

    ! add contributions from Coulomb and exchange potentials
    do iorb=1,norb
       ! it is asumed that exchange potentials involving iorb1
       ! has already been retrieved from disk
       coo=occ(iorb)
       kex=iorb1+norb*(iorb-1)
          
       if (iorb1.le.iorb) then
          ihc=iorb1+iorb*(iorb-1)/2
       else
          ihc=iorb+iorb1*(iorb1-1)/2
       endif
       
       if (iorb1.eq.iorb) coo=coo-one
       call dcopy (mxsize,pot(i2b(iorb)),ione,wk1,ione)
       call dscal (mxsize,coo,wk1,ione)
       call prod (mxsize,psi(i1b(iorb1)),wk1)
       
       call add (mxsize,wk1,wk0)
       
       if (iorb.ne.iorb1)  then
          call prodas (mxsize,-gec(kex),psi(i1b(iorb)),excp(i3b(ihc)),wk2)
          if (ilc(ihc).gt.1) then
             call prodas (mxsize,-gec(kex+norb*norb),psi(i1b(iorb)),excp(i3b(ihc)+mxsize),wk2)
          endif
       else
          if ((mm(iorb).gt.0).and.(ilc(ihc).gt.0)) then
             call prodas(mxsize,-gec(kex),psi(i1b(iorb)),excp(i3b(ihc)),wk2)
          endif
       endif
    enddo
     
    ! to complete the integrand wk2 has to be multiplied by psi(i1beg)
    call prod (mxsize,psi(i1b(iorb2)),wk0)
    wtwoelCoul=ddot(mxsize,wgt2,ione,wk0,ione)
    
    call prod (mxsize,psi(i1b(iorb2)),wk2)
    wtwoelExch=ddot(mxsize,wgt2,ione,wk2,ione)
    
    wtwoel=wtwoelCoul+wtwoelExch
    !call prod (ngorb,psi(i1beg),wk2)
    !wtwoel=ddot(ngorb,wgt2,ione,wk2,ione)
    
    twoelij=wtwoelCoul+wtwoelExch

    ! if(iprint(64).ne.0) then
    !    write(*,'(4x,"<",i2,1x,a5,a1,"|J+K|",i2,1x,a5,a1"> =",1Pe23.16)') &
    !         iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoel
    ! endif

#ifdef PRINT
! print= 64: twoelij: 2-electron contributions
    if (iprint(64).ne.0) then
       write(*,'(4x,"twoelij: energy contributions for",i4,1x,a8,a1)') iorn(iorb1),bond(iorb1),gusym(iorb1)
       write(*,'("    two-electron-Coul         =",1Pe23.16)') wtwoelCoul
       write(*,'("    two-electron-Exch         =",1Pe23.16)') wtwoelExch
       write(*,'("    two-electron              =",1Pe23.16)') twoelij
    endif
#endif

#ifdef PRINT
! print= 65: twoelij: 2-electron contributions <i|J|j> and <i|K|j>
    if(iprint(65).ne.0) then
       !write(*,'(18x,"DFT: 2-electron energy contributions:")')
       write(*,'(4x,"<",i2,1x,a5,a1,"| J |",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoelCoul
       write(*,'(4x,"<",i2,1x,a5,a1,"| K |",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoelExch
       write(*,'(4x,"<",i2,1x,a5,a1,"|J+K|",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoel
    endif
#endif
    
  end function twoelij

  ! ### twoelijLXC ###
  !
  !     FIXME
  !
  real (PREC) function twoelijLXC (iorb1,iorb2,psi,pot,excp,wgt2,wk0,wk1,wk2)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use inout
    use sharedMemory
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ihc,iorb,iorb1,iorb2,kex
    integer (KIND=IPREC) :: ibpot

    real (PREC) :: coo,w,woneel,wtwoel,wtwoelCoul,wtwoelCoulExch,wtwoelExch,&
         wtwoelHyb1,wtwoelHyb2,twoijDFT
    real (PREC) :: wtwoelCoul1,wtwoelExch1
    real (PREC), dimension(*) :: psi,pot,excp,wgt2,wk0,wk1,wk2
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    call zeroArray (mxsize,wk0)
    call zeroArray (mxsize,wk1)
    call zeroArray(mxsize,wk2)
    
    wtwoel=zero
    wtwoelExch=zero
    wtwoelHyb1=zero
    wtwoelHyb2=zero

    ! Add contributions from Coulomb and local exchange potential.
    ! In the local exchange approximation the coulomb potential
    ! includes also the contribution from the orbital
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) cycle
       if (inDFT(iorb).eq.0) cycle
       call daxpy (mxsize,occ(iorb),pot(i2b(iorb)),ione,wk2,ione)
    enddo
    call prod  (mxsize,psi(i1b(iorb1)),wk2)
    call prod (mxsize,psi(i1b(iorb2)),wk2)
    wtwoelCoul=ddot(mxsize,wgt2,ione,wk2,ione)

    ! multiply the local exchange/correlation potential by psi(i1beg)psi(i2beg) and
    ! add the result to the Coulomb potential

    !call prod2 (mxsize,psi(i1b(iorb1)),pot(length2-2*mxsize+1),wk1)
    !call prod (mxsize,psi(i1b(iorb2)),wk1)

    call prod2 (mxsize,psi(i1b(iorb1)),coulombptr(i2b(iorb1):),wk1)
    call prod (mxsize,psi(i1b(iorb2)),wk1)
    wtwoelHyb1=ddot(mxsize,wgt2,ione,wk1,ione)

    !???? wtwoelHyb1 already contains contributions from -alpha times  <i|K(i)|i>=<i|J(i)|i> exchange terms

    ! Now, wtwoelHyb1 contains ONLY contribution from DFT functional
    wtwoelCoul=wtwoelCoul+wtwoelHyb1
    wtwoel=wtwoelCoul
    if (lxcHyb) then
       call dcopy (mxsize,pot(i2b(iorb1)),ione,wk2,ione)
       call prod  (mxsize,psi(i1b(iorb1)),wk2)
       call prod (mxsize,psi(i1b(iorb2)),wk2)
       wtwoelCoulExch=-alphaf*ddot(mxsize,wgt2,ione,wk2,ione)
       ! wtwoelCoulExch contains contribution from -alpha times  <i|K(i)|i>=<i|J(i)|i> exchange terms       
       !call dcopy (mxsize,pot(length2-mxsize+1),ione,wk1,ione)

       call dcopy (mxsize,exchangeptr(i2b(iorb2):),ione,wk1,ione)   
       call prod (mxsize,psi(i1b(iorb2)),wk1)
       ! wtwoelHyb2 contains ONLY contribution from -alpha <i|K(j)|i>, i/=j terms
       wtwoelHyb2=-ddot(mxsize,wgt2,ione,wk1,ione)

       print *,"wtwoelHyb2",iorb1,iorb2,wtwoelHyb2

       call dcopy (mxsize,exchangeptr(i2b(iorb2):),ione,wk1,ione)   
       call prod (mxsize,psi(i1b(iorb1)),wk1)
       ! wtwoelHyb2 contains ONLY contribution from -alpha <i|K(j)|i>, i/=j terms
       wtwoelHyb2=-ddot(mxsize,wgt2,ione,wk1,ione)
       print *,"wtwoelHyb2",iorb1,iorb2,wtwoelHyb2
       
       !wtwoel=wtwoel+wtwoelCoulExch+wtwoelHyb2
       !wtwoel=wtwoel+wtwoelHyb2
       ! exchange contribution is already included in
       ! pot(length2-2*mxsize+1),wk1)
       ! no need to include it again
       !wtwoel=wtwoel+wtwoelCoulExch
       ! add contributions from exchange potentials

       call zeroArray(mxsize,wk2)       
       do iorb=1,norb
          !coo=occ(iorb)
          kex=iorb1+norb*(iorb-1)
          
          if (iorb1.le.iorb) then
             ihc=iorb1+iorb*(iorb-1)/2
          else
             ihc=iorb+iorb1*(iorb1-1)/2
          endif
          
          !if (iorb1.eq.iorb) coo=coo-one
          !call dcopy (mxsize,pot(i2b(iorb)),ione,wk1,ione)
          !call dscal (mxsize,coo,wk1,ione)
          !call prod (mxsize,psi(i1b(iorb1)),wk1)
          
          !call add (mxsize,wk1,wk0)
          
          if (iorb.ne.iorb1)  then
             call prodas (mxsize,-gec(kex),psi(i1b(iorb)),excp(i3b(ihc)),wk2)
             if (ilc(ihc).gt.1) then
                call prodas (mxsize,-gec(kex+norb*norb),psi(i1b(iorb)),excp(i3b(ihc)+mxsize),wk2)
             endif
          else
             if ((mm(iorb).gt.0).and.(ilc(ihc).gt.0)) then
                call prodas(mxsize,-gec(kex),psi(i1b(iorb)),excp(i3b(ihc)),wk2)
             endif
          endif
       enddo
     
       ! to complete the integrand wk2 has to be multiplied by psi(i1beg)
    
       call prod (mxsize,psi(i1b(iorb2)),wk2)
       wtwoelExch=alphaf*ddot(mxsize,wgt2,ione,wk2,ione)
       !print *,"3",wtwoelExch,-alphaf*wtwoelCoul/occ(iorb1)
       ! contribution i the Coulomb part
       !wtwoel=wtwoel+wtwoelExch-alphaf*wtwoelCoul/occ(iorb1)
       !wtwoel=wtwoel+wtwoelExch-alphaf*wtwoelCoul/occ(iorb1)
       wtwoel=wtwoel+wtwoelExch
    endif
 
    twoelijLXC=wtwoel

#ifdef PRINT
! print= 64: twoelijDFT: 2-electron contributions 
    if (iprint(64).ne.0) then
       !write(*,'("twoelijDFT: energy contributions")') 
       write(*,'(4x,"twoelijDFT: energy contributions for",i4,1x,a8,a1)') iorn(iorb1),bond(iorb1),gusym(iorb1)
       write(*,'("    two-electron-Coul         =",1Pe23.16)') wtwoelCoul
       write(*,'("    two-electron-CoulExch     =",1Pe23.16)') wtwoelCoulExch
       write(*,'("    two-electron-Exch         =",1Pe23.16)') wtwoelExch
       write(*,'("    two-electron-hyb1         =",1Pe23.16)') wtwoelHyb1
       write(*,'("    two-electron-hyb2         =",1Pe23.16)') wtwoelHyb2
       write(*,'("    two-electron              =",1Pe23.16)') wtwoel
    endif
#endif

        
#ifdef PRINT
! print= 65: twoelijDFT: 2-electron contributions <i|J|j> and <i|K|j>
    if(iprint(65).ne.0) then
       !write(*,'(18x,"DFT: 2-electron energy contributions:")')
       write(*,'(4x,"<",i2,1x,a5,a1,"| J |",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoelCoul
       write(*,'(4x,"<",i2,1x,a5,a1,"| K |",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoelExch
       write(*,'(4x,"<",i2,1x,a5,a1,"|J+K|",i2,1x,a5,a1"> =",1Pe23.16)') &
            iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),wtwoel
    endif
#endif
    
    call dcopy (mxsize,coulombptr(i2b(iorb2):),ione,wk1,ione)
    call prod (mxsize,psi(i2b(iorb2)),wk1)
    call prod (mxsize,psi(i2b(iorb1)),wk1)
    wtwoelHyb1=ddot(mxsize,wgt2,ione,wk1,ione)


#ifdef PRINT
! print= 64: twoelijDFT: 2-electron hybrid contribution
    if (iprint(64).ne.0) then
       !write(*,'("twoelijDFT: energy contributions")') 
       write(*,'(4x,"twoelijDFT: energy contributions for",i4,1x,a8,a1)') iorn(iorb1),bond(iorb1),gusym(iorb1)
       write(*,'("    two-electron-hyb1         =",1Pe23.16)') wtwoelHyb1
    endif
#endif
    
    ! call dcopy (mxsize,exchangeptr(i2b(iorb2):),ione,wk1,ione)
    ! call prod (mxsize,psi(i1b(iorb1)),wk1)
    ! wtwoelExch1=ddot(mxsize,wgt2,ione,wk1,ione)

end function twoelijLXC


end module integrals
