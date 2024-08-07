! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module sormcsor
  implicit none
contains
  ! ### sor ###
  !
  !     Performs one iteration of the sor process applied to
  !     the equation (diagonal part of the differentiation operator, i.e.
  !     the coefficient multiplying fun(ni,mi) is included in lhs array)
  !
  !    { d2(mi) + b(mi)*d1(mi) + d2(ni) + d(ni)*d1(ni)
  !		      + lhs(ni,mi) } fun(ni,mi) = rhs(ni,mi)
  !
  !     warning! this routine uses mesh defined by mesh
  !
  subroutine sor (isym,fun,lhs,rhs,b,d,indx1,indx2,indx6a,indx6b,indx7)
    use params
    use discrete
    use solver
    use commons
    
    implicit none
    integer (KIND=IPREC) :: i,ij,ik,itr2,isym
    real (PREC) :: cc,ddmi1,ddmi2,ddni1,ddni2
    !real (PREC) :: fcurr,fprev
    
    integer (KIND=IPREC),dimension(*) :: indx1,indx2,indx6a,indx6b,indx7
    real (PREC), dimension(*) :: fun,lhs,rhs,b,d

    ! indx1 (indx2) array elements point to the particular (consequtive)
    ! element of fun (lhs, rhs, b and d) to be relaxed

    do itr2=1,maxsor2
       do i=isstart,isstop,isstep
          ij=indx1(i)
          ik=indx2(i)

          ddmi2=                                              &
               dmu2(1) * ( fun(ij-nni4) + fun(ij+nni4) )+    &
               dmu2(2) * ( fun(ij-nni3) + fun(ij+nni3) )+    &
               dmu2(3) * ( fun(ij-nni2) + fun(ij+nni2) )+    &
               dmu2(4) * ( fun(ij-nni1) + fun(ij+nni1) )

          ddmi1=                                              &
               dmu1(1) * ( fun(ij-nni4) - fun(ij+nni4) )+    &
               dmu1(2) * ( fun(ij-nni3) - fun(ij+nni3) )+    &
               dmu1(3) * ( fun(ij-nni2) - fun(ij+nni2) )+    &
               dmu1(4) * ( fun(ij-nni1) - fun(ij+nni1) )

          ddni2=                                              &
               dni2(1) * ( fun(ij-   4) + fun(ij+   4) )+     &
               dni2(2) * ( fun(ij-   3) + fun(ij+   3) )+     &
               dni2(3) * ( fun(ij-   2) + fun(ij+   2) )+     &
               dni2(4) * ( fun(ij-   1) + fun(ij+   1) )

          ddni1=                                              &
               dni1(1) * ( fun(ij-   4) - fun(ij+   4) )+     &
               dni1(2) * ( fun(ij-   3) - fun(ij+   3) )+     &
               dni1(3) * ( fun(ij-   2) - fun(ij+   2) )+     &
               dni1(4) * ( fun(ij-   1) - fun(ij+   1) )

          cc=ddmi2 + b(ik)*ddmi1 + ddni2 + d(ik)*ddni1
          fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)

          !write(*,'("sor:",3i8,6e16.6)') i,ij,ik,b(ik),d(ik),rhs(ik),lhs(ik),cc,fun(ij)
          !if (mod(i,5000)==0) write(*,'("sor:",3i8,6e16.6)') i,ij,ik,b(ik),d(ik),rhs(ik),lhs(ik),cc,fun(ij)
          !fprev=fun(ij)
          !fcurr = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)        
          !fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)
          !if (abs(fcurr-fprev)>9.93e-5) write(*,'("sor:",2i8,4e16.6,e16.6)') &    
          !     indx1nu(ij),indx1mu(ij),lhs(ik),rhs(ik),cc,fcurr,fcurr-fprev
       enddo

       !   determine values at the boundary points from extrapolation
       !   in version 2.0 single grid implies ifill=1
       if (isym == 1) then
          do i=1,ngrid7
             ij=indx7(i)
             fun(ij)=exeven(1)*fun(ij+nni1)+exeven(2)*fun(ij+nni2)+   &
                  exeven(3)*fun(ij+nni3)+exeven(4)*fun(ij+nni4)+exeven(5)*fun(ij+nni5)
             !write(*,'(i5,6e12.4)') ij,fun(ij),fun(ij+nni1),fun(ij+nni2),fun(ij+nni3),fun(ij+nni4),fun(ij+nni5)
             !write(*,'(i5,6e12.4)') ij,exeven(1),exeven(2),exeven(3),exeven(4),exeven(5)
          enddo

          do i=1,ngrid6a
             ij=indx6a(i)
             fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+  &
                  exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
          enddo

          do i=1,ngrid6b
             ij=indx6b(i)
             fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+  &
                  exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
          enddo
       else
          do i=1,ngrid7
             ij=indx7(i)
             fun(ij)=0.0_PREC
          enddo

          do i=1,ngrid6a
             ij=indx6a(i)
             fun(ij)=0.0_PREC
          enddo

          do i=1,ngrid6b
             ij=indx6b(i)
             fun(ij)=0.0_PREC
          enddo
       endif
    enddo
  end subroutine sor

  ! ### mcsor ###
  !
  !     Performs one iteration of the sor process applied to the equation (diagonal
  !     part of the differentiation operator, i.e.  the coefficient multiplying
  !     fun(ni,mi) is included in lhs array)
  !
  !     { d2(mi) + b(mi)*d1(mi) + d2(ni) + d(ni)*d1(ni)
  !                     + lhs(ni,mi) } fun(ni,mi) = rhs(ni,mi)
  !
  subroutine mcsor (isym,fun,lhs,rhs,b,d,indx1,indx2,indx6a,indx6b,indx7)
    use params
    use discrete
    use solver
#ifdef OPENMP
    use omp_lib
#endif 
    implicit none
    integer (KIND=IPREC) :: i,ij,ik,itr2,isym
    real (PREC) :: cc,ddmi1,ddmi2,ddni1,ddni2

    integer (KIND=IPREC),dimension(*) :: indx1,indx2,indx6a,indx6b,indx7
    real (PREC), dimension(*) :: fun,lhs,rhs,b,d

    !     indx1 (indx2) array elements point to the particular (consequtive)
    !     element of fun (lhs, rhs, b and d) to be relaxed
    !
    !     first relax the inner points -- region i (formula a1)
    !
    !     ipoiss=2 -- inner variables are relaxed in 5 sweeps each relaxing
    !     every 5th variable.
    !
    !     Open MPI directives are used to force parallelization of the loops.

    do itr2=1,maxsor2
       !$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(SHARED) PRIVATE(i,ij,ik,ddmi2,ddmi1,ddni2,ddni1,cc)
       !$OMP DO SCHEDULE(STATIC) 
       do i=1,ngrid1,5
          ij=indx1(i)
          ik=indx2(i)
          ddmi2=dmu2(1) * ( fun(ij-nni4) + fun(ij+nni4) )+  &
               dmu2(2) * ( fun(ij-nni3) + fun(ij+nni3) )+   &
               dmu2(3) * ( fun(ij-nni2) + fun(ij+nni2) )+   &
               dmu2(4) * ( fun(ij-nni1) + fun(ij+nni1) )
          
          ddmi1=dmu1(1) * ( fun(ij-nni4) - fun(ij+nni4) )+  &
               dmu1(2) * ( fun(ij-nni3) - fun(ij+nni3) )+   &
               dmu1(3) * ( fun(ij-nni2) - fun(ij+nni2) )+   &
               dmu1(4) * ( fun(ij-nni1) - fun(ij+nni1) )
          
          ddni2=dni2(1) * ( fun(ij-   4) + fun(ij+   4) )+   &
               dni2(2) * ( fun(ij-   3) + fun(ij+   3) )+    &
               dni2(3) * ( fun(ij-   2) + fun(ij+   2) )+    &
               dni2(4) * ( fun(ij-   1) + fun(ij+   1) )
          
          ddni1=dni1(1) * ( fun(ij-   4) - fun(ij+   4) )+   &
               dni1(2) * ( fun(ij-   3) - fun(ij+   3) )+    &
               dni1(3) * ( fun(ij-   2) - fun(ij+   2) )+    &
               dni1(4) * ( fun(ij-   1) - fun(ij+   1) )
          
          cc=ddmi2 + b(ik)*ddmi1 + ddni2 + d(ik)*ddni1
          fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)
       enddo
       !$OMP END DO
       !$OMP FLUSH (FUN)
       !$OMP BARRIER

       !$OMP DO SCHEDULE(STATIC) 
       do i=2,ngrid1,5
          ij=indx1(i)
          ik=indx2(i)

          ddmi2=dmu2(1) * ( fun(ij-nni4) + fun(ij+nni4) )+   &
               dmu2(2) * ( fun(ij-nni3) + fun(ij+nni3) )+    &
               dmu2(3) * ( fun(ij-nni2) + fun(ij+nni2) )+    &
               dmu2(4) * ( fun(ij-nni1) + fun(ij+nni1) )
          
          ddmi1=dmu1(1) * ( fun(ij-nni4) - fun(ij+nni4) )+   &
               dmu1(2) * ( fun(ij-nni3) - fun(ij+nni3) )+    &
               dmu1(3) * ( fun(ij-nni2) - fun(ij+nni2) )+    &
               dmu1(4) * ( fun(ij-nni1) - fun(ij+nni1) )
          
          ddni2=dni2(1) * ( fun(ij-   4) + fun(ij+   4) )+    &
               dni2(2) * ( fun(ij-   3) + fun(ij+   3) )+     &
               dni2(3) * ( fun(ij-   2) + fun(ij+   2) )+     &
               dni2(4) * ( fun(ij-   1) + fun(ij+   1) )
          
          ddni1=dni1(1) * ( fun(ij-   4) - fun(ij+   4) )+    &
               dni1(2) * ( fun(ij-   3) - fun(ij+   3) )+     &
               dni1(3) * ( fun(ij-   2) - fun(ij+   2) )+     &
               dni1(4) * ( fun(ij-   1) - fun(ij+   1) )
          
          cc=ddmi2 + b(ik)*ddmi1 + ddni2 + d(ik)*ddni1
          fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)
       enddo
       !$OMP END DO
       !$OMP FLUSH (FUN)
       !$OMP BARRIER

       !$OMP DO SCHEDULE(STATIC) 
       do i=3,ngrid1,5
          ij=indx1(i)
          ik=indx2(i)

          ddmi2=dmu2(1) * ( fun(ij-nni4) + fun(ij+nni4) )+   &
               dmu2(2) * ( fun(ij-nni3) + fun(ij+nni3) )+    &
               dmu2(3) * ( fun(ij-nni2) + fun(ij+nni2) )+    &
               dmu2(4) * ( fun(ij-nni1) + fun(ij+nni1) )
          
          ddmi1=dmu1(1) * ( fun(ij-nni4) - fun(ij+nni4) )+   &
               dmu1(2) * ( fun(ij-nni3) - fun(ij+nni3) )+    &
               dmu1(3) * ( fun(ij-nni2) - fun(ij+nni2) )+    &
               dmu1(4) * ( fun(ij-nni1) - fun(ij+nni1) )
          
          ddni2=dni2(1) * ( fun(ij-   4) + fun(ij+   4) )+    &
               dni2(2) * ( fun(ij-   3) + fun(ij+   3) )+     &
               dni2(3) * ( fun(ij-   2) + fun(ij+   2) )+     &
               dni2(4) * ( fun(ij-   1) + fun(ij+   1) )
          
          ddni1=dni1(1) * ( fun(ij-   4) - fun(ij+   4) )+    &
               dni1(2) * ( fun(ij-   3) - fun(ij+   3) )+     &
               dni1(3) * ( fun(ij-   2) - fun(ij+   2) )+     &
               dni1(4) * ( fun(ij-   1) - fun(ij+   1) )
          
          cc=ddmi2 + b(ik)*ddmi1 + ddni2 + d(ik)*ddni1
          fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)
       enddo
       !$OMP END DO
       !$OMP FLUSH (FUN)
       !$OMP BARRIER

       !$OMP DO SCHEDULE(STATIC) 
       do i=4,ngrid1,5
          ij=indx1(i)
          ik=indx2(i)

          ddmi2=dmu2(1) * ( fun(ij-nni4) + fun(ij+nni4) )+    &
               dmu2(2) * ( fun(ij-nni3) + fun(ij+nni3) )+ &
               dmu2(3) * ( fun(ij-nni2) + fun(ij+nni2) )+ &
               dmu2(4) * ( fun(ij-nni1) + fun(ij+nni1) )
          
          ddmi1=dmu1(1) * ( fun(ij-nni4) - fun(ij+nni4) )+    &
               dmu1(2) * ( fun(ij-nni3) - fun(ij+nni3) )+     &
               dmu1(3) * ( fun(ij-nni2) - fun(ij+nni2) )+     &
               dmu1(4) * ( fun(ij-nni1) - fun(ij+nni1) )
          
          ddni2=dni2(1) * ( fun(ij-   4) + fun(ij+   4) )+     &
               dni2(2) * ( fun(ij-   3) + fun(ij+   3) )+      &
               dni2(3) * ( fun(ij-   2) + fun(ij+   2) )+      &
               dni2(4) * ( fun(ij-   1) + fun(ij+   1) )
          
          ddni1=dni1(1) * ( fun(ij-   4) - fun(ij+   4) )+     &
               dni1(2) * ( fun(ij-   3) - fun(ij+   3) )+      &
               dni1(3) * ( fun(ij-   2) - fun(ij+   2) )+      &
               dni1(4) * ( fun(ij-   1) - fun(ij+   1) )
          
          cc=ddmi2 + b(ik)*ddmi1 + ddni2 + d(ik)*ddni1
          fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)
       enddo
       !$OMP END DO 
       !$OMP FLUSH(FUN)       
       !$OMP BARRIER

       !$OMP DO SCHEDULE(STATIC) 
       do i=5,ngrid1,5
          ij=indx1(i)
          ik=indx2(i)

          ddmi2=dmu2(1) * ( fun(ij-nni4) + fun(ij+nni4) )+   &
               dmu2(2) * ( fun(ij-nni3) + fun(ij+nni3) )+    &
               dmu2(3) * ( fun(ij-nni2) + fun(ij+nni2) )+    &
               dmu2(4) * ( fun(ij-nni1) + fun(ij+nni1) )
          
          ddmi1=dmu1(1) * ( fun(ij-nni4) - fun(ij+nni4) )+   &
               dmu1(2) * ( fun(ij-nni3) - fun(ij+nni3) )+    &
               dmu1(3) * ( fun(ij-nni2) - fun(ij+nni2) )+    &
               dmu1(4) * ( fun(ij-nni1) - fun(ij+nni1) )
          
          ddni2=dni2(1) * ( fun(ij-   4) + fun(ij+   4) )+    &
               dni2(2) * ( fun(ij-   3) + fun(ij+   3) )+     &
               dni2(3) * ( fun(ij-   2) + fun(ij+   2) )+     &
               dni2(4) * ( fun(ij-   1) + fun(ij+   1) )
          
          ddni1=dni1(1) * ( fun(ij-   4) - fun(ij+   4) )+    &
               dni1(2) * ( fun(ij-   3) - fun(ij+   3) )+     &
               dni1(3) * ( fun(ij-   2) - fun(ij+   2) )+     &
               dni1(4) * ( fun(ij-   1) - fun(ij+   1) )
          
          cc=ddmi2 + b(ik)*ddmi1 + ddni2 + d(ik)*ddni1
          fun(ij) = omega * (rhs(ik)-cc)/lhs(ik)+ omega1 * fun(ij)
       enddo
       !$OMP END DO
       !$OMP FLUSH(FUN)
       !$OMP BARRIER

       !$OMP END PARALLEL 

       ! determine values at the bondary points from extrapolation 
       ! implies ifill=1

       !$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(SHARED) PRIVATE(i,ij)
       if (isym.eq.1) then
          !$OMP DO
          do i=1,ngrid7
             ij=indx7(i)
             fun(ij)=exeven(1)*fun(ij+nni1)+exeven(2)*fun(ij+nni2)+ &
                  exeven(3)*fun(ij+nni3)+exeven(4)*fun(ij+nni4)+exeven(5)*fun(ij+nni5)
          enddo
          !$OMP END DO

          !$OMP DO          
          do i=1,ngrid6a
             ij=indx6a(i)
             fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+    &
                  exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
          enddo
          !$OMP END DO

          !$OMP DO
          do i=1,ngrid6b
             ij=indx6b(i)
             fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+    &
                  exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
          enddo
          !$OM END DO
       else
          !$OMP DO
          do i=1,ngrid7
             ij=indx7(i)
             fun(ij)=0.0_PREC
          enddo
          !$OMP END DO
                    
          !$OMP DO
          do i=1,ngrid6a
             ij=indx6a(i)
             fun(ij)=0.0_PREC
          enddo
          !$OMP END DO

          !$OMP DO
          do i=1,ngrid6b
             ij=indx6b(i)
             fun(ij)=0.0_PREC
          enddo
          !$OMP END DO
       endif
       !$OMP FLUSH
       !$OMP BARRIER

       !$OMP END PARALLEL  
    enddo

  end subroutine mcsor
end module sormcsor
