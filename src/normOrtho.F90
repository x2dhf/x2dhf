! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module normOrtho
  implicit none
contains

  ! ### norm ###
  !
  !     Normalizes a given orbital.
  
  subroutine norm(iorb)
    use params
    use commons
    use discrete 
    use utils
    use blas
    use sharedMemory
    use discrete
    implicit none
    integer (KIND=IPREC) :: i,ibeg,iorb
    real (PREC) :: xnorm,zarea
    real (PREC), dimension(:), pointer :: f4,psi,psi1,psi3,wgt2,wk0

#ifdef BLAS
    real (PREC) ddot
    external ddot
#endif
    
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr

    !if (ifixorb(iorb).ne.0) return

    ibeg = i1b (iorb)

    do i=1,mxsize
       wk0(i)=psi(ibeg-1+i)*psi(ibeg-1+i)
    enddo
    call prod  (mxsize,f4,wk0)

    xnorm=ddot (mxsize,wgt2,ione,wk0,ione)
    orbNorm(iorb)=xnorm
    zarea = one/sqrt(xnorm)
    psi1=>orbptr(ibeg:)
    call dscal (mxsize,zarea,psi1,ione)

#ifdef PRINT
! print= 36: norm: checking norms of orbitals    
    if (iprint(36).ne.0) then
       write(*,'("1-norm: ",i4,1x,a8,a1,3x,e16.2)') iorn(iorb),bond(iorb),gusym(iorb),abs(one-orbNorm(iorb))
    endif
#endif    

  end subroutine norm
  
  ! ### norm94 ###
  !
  subroutine norm94 (iorb,xnorm)
    use params
    use commons
    use discrete
    use sharedmemory
    use utils
    use blas

    implicit none
    integer (KIND=IPREC) :: i,ibeg,iorb
    real (PREC) :: xnorm
    real (PREC), dimension(:), pointer :: f4,psi,psi1,psi3,wgt2,wk0

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr


    ibeg = i1b (iorb)
    do i=1,mxsize
       wk0(i)=psi(ibeg-1+i)*psi(ibeg-1+i)
    enddo

    call prod  (mxsize,f4,wk0)

    xnorm=ddot (mxsize,wgt2,ione,wk0,ione)
    orbNorm(iorb)=xnorm

  end subroutine norm94

  ! ### ortho ###
  !
  !     Orthogonalizes a given orbital using the Schmidt algorithm
  !
  subroutine ortho (iorb1)
    use params
    use commons
    use discrete
    use sharedMemory
    use utils
    use blas

    implicit none
    integer (KIND=IPREC) :: ibeg1,ibeg3,iorb1,iorb2,iorb3,ira,jor
    real (PREC) :: ano
    integer (KIND=IPREC),dimension(maxorb) :: index,istp
    real (PREC), dimension(maxorb) :: ovla,ovla1
    real (PREC), dimension(:), pointer :: f4,psi,psi1,psi3,wgt2,wk0

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr
    
    if (ifixorb(iorb1).ne.0) return

    jor=0
    ibeg1 = i1b (iorb1)
    psi1=>orbptr(ibeg1:)
    ! orbitals are stored in reverse order
    do iorb2=1,norb
       iorb3=norb-iorb2+1
       if (iorb3.le.iorb1) cycle
       if (mgx(6,iorb3).ne.mgx(6,iorb1)) cycle
       if (ige(iorb3).ne.ige(iorb1).and. .not.lbreakCi) cycle
       jor=jor+1
       ibeg3= i1b (iorb3 )
       istp(jor)=ibeg3
       psi3=>orbptr(ibeg3:)
       call prod2 (mxsize,psi1,psi3,wk0)
       call prod  (mxsize,f4,wk0)

       ovla(jor)=ddot (mxsize,wgt2,ione,wk0,ione)
       ovla1(jor)=abs(ovla(jor))

       ! print overlap integrals
#ifdef PRINT
! print= 20: ortho: checking overlap of orbitals
       if (iprint(20).ne.0) then
          write(*,60) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb3),bond(iorb3),gusym(iorb3),ovla(jor)
60        format(4x,' <',i2,1x,a5,1x,a1,1x,'|',1x,i2,1x,a5,1x,a1,1x,'> = ',1Pe10.2)
       endif
#endif
       
       ! print overlap integrals only for pairs of orbitals having
       ! non-zero off-diagonal Lagrange multipliers
#ifdef PRINT
! print= 21: ortho: overlap integrals for pairs of orbitals with non-zero off-diagonal Lagrange multipliers
       if (iprint(21).ne.0.and.nlm(iorb1,iorb3).ne.0) then
          write(*,60) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb3),bond(iorb3),gusym(iorb3),ovla(jor)
       endif
#endif
    enddo

    ! determine the worst nonorthogonality for a given orbital
    if (jor.gt.0) then
       call quicksort(jor,ovla1,index)
       wstorthog(iorb1)=ovla1(index(jor))
    endif

    ! start orthogonalization
    if (jor.eq.0) return

    ano=zero
    do ira=1,jor
       ano=ano+ovla(ira)*ovla(ira)
       psi3=>orbptr(istp(ira):)
       call daxpy (mxsize,-ovla(ira),psi3,ione,psi1,ione)
    enddo

    ! normalize the orthogonalized orbital isa  but be sure the
    ! unorthogonalized isa was normalized otherwise use <isa/isa>
    ! instead of 1. in ano
    ano=sqrt(one+ano)/(one-ano*ano)
    call dscal (mxsize,ano,psi1,ione)

    return
    
    ! Except for the lowest orbitals in any given symmetry the
    ! normalization factors are not equal to unity through the scf
    ! process. These orbitals aquire proper normalization upon completion
    ! of the scf process. In order to force the normalization to unity
    ! set iprint(35) to 1.

    ! if (iprint(35).eq.0) return

    ! call prod2 (mxsize,psi1,psi1,wk0)
    ! call prod  (mxsize,f4,wk0)

    ! ovla(jor)=ddot (mxsize,wgt2,ione,wk0,ione)
    ! ano=one/sqrt(ovla(jor))
    ! call dscal (mxsize,ano,psi1,ione)

  end subroutine ortho

  ! ### checkOrtho ###
  !
  !     Checks orthogonalization of a given orbital
  !
  subroutine checkOrtho (lprint,iorb1)
    use blas
    use params
    use commons
    use discrete
    use sharedMemory
    use utils

    implicit none
    logical :: lprint
    integer (KIND=IPREC) :: ibeg1,ibeg3,iorb1,iorb2,iorb3,jor,jor1

    integer (KIND=IPREC),dimension(maxorb) :: istp
    real (PREC), dimension(maxorb) :: ovla,ovla1
    real (PREC), dimension(:), pointer :: f4,psi,wgt2,wk0
    
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    f4=>supplptr(i4b(9):)
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr
    
    jor=0
    jor1=0
    ibeg1 = i1b (iorb1)

    ! orbitals are stored in reverse order
    do iorb2=1,norb
       iorb3=norb-iorb2+1
       if (iorb3.le.iorb1) cycle

       if (mgx(6,iorb3).ne.mgx(6,iorb1)) cycle
       if (ige(iorb3).ne.ige(iorb1).and. .not.lbreakCi) cycle
       
       jor=jor+1
       ibeg3= i1b (iorb3 )
       istp(jor)=ibeg3

       call prod2 (mxsize,psi(ibeg1:),psi(ibeg3:),wk0)
       call prod  (mxsize,f4,wk0)

       ovla(jor)=ddot(mxsize,wgt2,ione,wk0,ione)

       jor1=jor1+1
       ovla1(jor1)=abs(ovla(jor))

       if (lprint) then
          ! print overlap integrals
          write(*,1000) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb3),bond(iorb3),gusym(iorb3),ovla(jor)
          ! print overlap integrals only for pairs of orbitals having
          ! non-zero off-diagonal Lagrange multipliers
          if (nlm(iorb1,iorb3).ne.0) then
             write(*,1000) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb3),bond(iorb3),gusym(iorb3),ovla(jor)
          endif
       endif
    enddo
1000 format(4x,' <',i2,1x,a5,1x,a1,1x,'|',1x,i2,1x,a5,1x,a1,1x,'> = ',e8.2)
  end subroutine checkOrtho


  ! ### rotate ###
  !
  !     Rotate a pair of orbitals by theta
  !
  subroutine rotate (theta,iorb1,iorb2)
    use params
    use commons
    use discrete
    use sharedMemory
    use utils
    use blas

    implicit none
    integer (KIND=IPREC) :: i,ibeg1,ibeg2,iorb1,iorb2
    real (PREC) :: theta,thetaRad
    real (PREC), dimension(:), pointer :: psi1,psi2,wk1,wk2

    wk1  =>scratchptr(         1:   mxsize8)
    wk2  =>scratchptr( mxsize8+1: 2*mxsize8)
    
    ibeg1= i1b (iorb1 )
    ibeg2= i1b (iorb2 )    

    psi1=>orbptr(ibeg1:)
    psi2=>orbptr(ibeg2:)
    thetaRad=theta/180.0_PREC*pii
    do i=1,mxsize
       wk1(i)=cos(thetaRad)*psi1(i) - sin(thetaRad)*psi2(i)
       wk2(i)=sin(thetaRad)*psi1(i) + cos(thetaRad)*psi2(i)       
    enddo

    do i=1,mxsize
       psi1(i)=wk1(i)
       psi2(i)=wk2(i)
    enddo

    print *,"rotate: theta orb1 orb2",theta,iorb1,iorb2
  end subroutine rotate

  
end module normOrtho
