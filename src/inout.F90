! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2023  Jacek Kobus 

module inout
  implicit none
 
contains
  ! ### putin ###
  !
  !     Immerses FUN array into WORK and provides missing values along (1,i),
  !     (nni,i) and (j,1) boundaries according to orbital symmetry isym.
  !
  !     The SOR relaxations are carried out for grid points up to (j,nmi),
  !     j=1,nni+8 employing the boundary conditions at (j,nmi+1, (j,nmi+2),
  !     (j,nmi+3) and (j,nmi+4) grid points.
  !
  !     In case of orbitals the orbAsymptSet routine is used to update these
  !     values due to the known asymptotic behaviour of the functions. In case
  !     of Coulomb and exchange potentials the multipole expansion is used to
  !     provide the boundary values (see potAsympt).

  subroutine putin (nni,nmi,isym,fun,work)
    use params
    use blas
    
    implicit none
    integer (KIND=IPREC) :: i,j,isym,jj,nni,nmi

    real (PREC), dimension(nni,nmi) :: fun
    real (PREC), dimension(nni+8,nmi+8) :: work

    ! fill the interior of work array
    do i=1,nmi
       call dcopy(nni,fun(1,i),ione,work(5,i+4),ione)
    enddo

    ! do i=1,4
    !    do j=2,nni-1
    !       work(4+j,nmi+4+i)=fun(j,nmi)
    !    enddo
    ! enddo
    
    ! isym = 1 - even symmetry, isym =-1 - odd symmetry
    ! values over i=nmi boundary are determined from the asymptotic expansion

    ! mu=nmu+5....nmu+8
    ! the following code is necessary since the derivatives must be
    ! calculated up to i=nmi

    if (isym.eq.1) then
       ! mu=1...4
       do i=2,5
          call dcopy(nni,fun(1,i),ione,work(5,6-i),ione)
       enddo

       do j=2,5
          call dcopy(nmi,fun(j,1),nni,work(6-j,5),nni+8)
       enddo

       do j=1,4
          call dcopy(nmi,fun(nni-j,1),nni,work(nni+4+j,5),nni+8)
       enddo
    else
       do i=2,5
          call dcopyi(nni,fun(1,i),ione,work(5,6-i),ione)
       enddo

       do j=2,5
          call dcopyi(nmi,fun(j,1),nni,work(6-j,5),nni+8)
       enddo
       
       do j=1,4
          call dcopyi(nmi,fun(nni-j,1),nni,work(nni+4+j,5),nni+8)
       enddo
    endif

    ! do j=1,nni+8                                                                                        
    !    if (mod(j,10)==0) write(*,'(i5,5e16.8)')j,work(j,nmi+4),work(j,nmi+4+1),work(j,nmi+4+2)
    ! enddo

  end subroutine putin

  subroutine putin1 (nni,nmi,isym4nu,isym4mu,fun,work)
    use params
    use blas

    implicit none
    integer (KIND=IPREC) :: i,j,isym4nu,isym4mu,jj,nni,nmi

    real (PREC), dimension(nni,nmi) :: fun
    real (PREC), dimension(nni+8,nmi+8) :: work

    ! fill the interior of work array
    do i=1,nmi
       call dcopy(nni,fun(1,i),ione,work(5,i+4),ione)
    enddo

    do i=nmi+1,nmi+4
       do j=1,nni
          work(j+4,i+4)=fun(j,nmi)
       enddo
    enddo

    ! isym = 1 - even symmetry, isym =-1 - odd symmetry
    ! values over i=nmi boundary are determined from the asymptotic expansion

    ! mu=nmu+5....nmu+8
    ! the following code is necessary since the derivatives must be
    ! calculated up to i=nmi

    if (isym4mu.eq.1) then
       ! mu=1...4
       do i=2,5
          do j=1,nni
             work(j+4,6-i)= fun(j,i)
          enddo
       enddo
    else
       do i=2,5
          do j=1,nni
             work(j+4,6-i)=-fun(j,i)
          enddo
       enddo
    endif

    if (isym4nu.eq.1) then
       ! ni=1...4
       do i=1,nmi
          do j=2,5
             work(6-j,i+4)= fun(j,i)
          enddo
       enddo

       ! ni=ni+4...ni+8
       do i=1,nmi
          jj=0
          do j=nni-4,nni-1
             jj=jj+1
             work(nni+9-jj,i+4)= fun(j,i)
          enddo
       enddo
    else
       do i=1,nmi
          do j=2,5
             work(6-j,i+4)=-fun(j,i)
          enddo
       enddo

       do i=1,nmi
          jj=0
          do j=nni-4,nni-1
             jj=jj+1
             work(nni+9-jj,i+4)=-fun(j,i)
          enddo
       enddo
    endif

  end subroutine putin1

  ! ### putout ###
  !
  !     Reverses the action of putin and putin[2-4]
  !
  subroutine putout (nni,nmi,fun,work)
    use params
    use blas
    
    implicit none
    integer (KIND=IPREC) :: i,j,nni,nmi

    real (PREC), dimension(nni,nmi) :: fun
    real (PREC), dimension(nni+8,nmi+8) :: work

    ! refill the interior of work array

    do i=1,nmi
       ! do j=1,nni
       !    fun(j,i)=work(j+4,i+4)
       ! enddo
       call dcopy(nni,work(5,i+4),ione,fun(1,i),ione)
    enddo
  end subroutine putout

#if defined (PTHREAD) || defined (TPOOL)

  ! ### putinc ###
  !
  !     Immerses FUN array into WORK and adds boundary values from symmetry
  !
  !     Immerses FUN array into WORK and provides missing values along (1,i), (nni,i) and (j,1)
  !     boundaries according to orbital symmetry isym.
  !     Values along (j,nmi) boundary are set to zero as orbitals decay exponentially in this region.
  !
  subroutine putinc (nni,nmi,isym,fun,work) bind(c,name='putinc_')
    !subroutine putinc (nni,nmi,isym,fun,work)
    use params
    use blas
    use iso_c_binding
    
    implicit none
    integer :: i,j,isym,jj,nni,nmi
    
    real (PREC), dimension(nni,nmi) :: fun
    real (PREC), dimension(nni+8,nmi+8) :: work
    
    ! interface
    !    subroutine putinc(nni,nmi,isym,fun,work) 
    !      integer(c_int), intent(in)  :: nni,nmi,isym
    !      real(kind=c_double), intent(in)  :: fun
    !      real(kind=c_double), intent(out)  :: work
    !    end subroutine putinc
    ! end interface
    
    
    ! fill the interior of work array
    do i=1,nmi
       call dcopyc(nni,fun(1,i),ione,work(5,i+4),ione)
    enddo
    
    ! isym = 1 - even symmetry, isym =-1 - odd symmetry
    ! values over i=nmi boundary are determined from the asymptotic expansion
    
    ! mu=nmu+5....nmu+8
    ! the following code is necessary since the derivatives must be
    ! calculated up to i=nmi
    
    if (isym.eq.1) then
       ! mu=1...4
       do i=2,5
          call dcopyc(nni,fun(1,i),ione,work(5,6-i),ione)
       enddo
       
       do j=2,5
          call dcopyc(nmi,fun(j,1),nni,work(6-j,5),nni+8)
       enddo
       
       do j=1,4
          call dcopyc(nmi,fun(nni-j,1),nni,work(nni+4+j,5),nni+8)
       enddo
    else
       do i=2,5
          call dcopyic(nni,fun(1,i),ione,work(5,6-i),ione)
       enddo
       
       do j=2,5
          call dcopyic(nmi,fun(j,1),nni,work(6-j,5),nni+8)
       enddo
       
       do j=1,4
          call dcopyic(nmi,fun(nni-j,1),nni,work(nni+4+j,5),nni+8)
       enddo
    endif

    do i=1,4                                                                                               
       do j=1,nni+8                                                                                        
          work(j,nmi+4+i)=work(j,nmi+4)                                                                    
       enddo                                                                                               
    enddo             
  end subroutine putinc
  
  ! ### putoutc ###
  !
  !     Reverses the action of putin and putin[2-4]
  !
  subroutine putoutc (nni,nmi,fun,work) bind(c,name='putoutc_')
    !subroutine putoutc_ (nni,nmi,fun,work) 
    use params
    use blas
    use iso_c_binding
    
    implicit none
    integer :: i,j,nni,nmi
    
    real (PREC), dimension(nni,nmi) :: fun
    real (PREC), dimension(nni+8,nmi+8) :: work
    
    ! interface
    !    subroutine putoutc(nni,nmi,isym,fun,work)) bind(c,name='putoutc_')
    !      integer(c_int), intent(in)  :: nni,nmi,isym
    !      real(kind=c_double), intent(out)  :: fun
    !      real(kind=c_double), intent(in)  :: work
    !    end subroutine putoutc
    ! end interface
    
    ! refill the interior of work array
    do i=1,nmi
       call dcopyc(nni,work(5,i+4),ione,fun(1,i),ione)
    enddo
  end subroutine putoutc
  
  subroutine  dcopyc(n,dx,ix,dy,iy)
    use params
    implicit none
    integer :: i,ix,iy,n
    real (PREC), dimension(*) :: dx,dy
    
    do i=1,n
       dy((i-1)*iy+1)=dx((i-1)*ix+1)
    enddo
  end subroutine dcopyc
  
  
  subroutine  dcopyic(n,dx,ix,dy,iy)
    use params
    implicit none
    integer :: i,ix,iy,n
    real (PREC), dimension(*) :: dx,dy
    do i=1,n
       dy((i-1)*iy+1)=-dx((i-1)*ix+1)
    enddo
  end subroutine dcopyic
  
#endif
  
end module inout
