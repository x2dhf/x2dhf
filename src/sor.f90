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

subroutine sor (fun,lhs,rhs,b,d,indx1,indx2,indx6a,indx6b,indx7)
  use params
  use discret
  use solver
  use commons8

  implicit none
  integer :: i,ij,ik
  real (PREC) :: cc,ddmi1,ddmi2,ddni1,ddni2

  integer, dimension(*) :: indx1,indx2,indx6a,indx6b,indx7
  real (PREC), dimension(*) :: fun,lhs,rhs,b,d



  !   indx1 (indx2) array elements point to the particular (consequtive)
  !   element of fun (lhs, rhs, b and d) to be relaxed

  !   first relax the inner points- region i (formula a1)

  do i=isstart,isstop,isstep
     ij=indx1(i)
     ik=indx2(i)

     ddmi2=                                              &
          dmu2t(1) * ( fun(ij-nni4) + fun(ij+nni4) )+    &
          dmu2t(2) * ( fun(ij-nni3) + fun(ij+nni3) )+    &
          dmu2t(3) * ( fun(ij-nni2) + fun(ij+nni2) )+    &
          dmu2t(4) * ( fun(ij-nni1) + fun(ij+nni1) )

     ddmi1=                                              &
          dmu1t(1) * ( fun(ij-nni4) - fun(ij+nni4) )+    &
          dmu1t(2) * ( fun(ij-nni3) - fun(ij+nni3) )+    &
          dmu1t(3) * ( fun(ij-nni2) - fun(ij+nni2) )+    &
          dmu1t(4) * ( fun(ij-nni1) - fun(ij+nni1) )

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
  enddo

  !   determine values at the bondary points from extrapolation
  !   in version 2.0 single grid implies ifill=1

  if (isym.eq.1) then
     !       if (ifill.eq.1.or.ifill.eq.2) then
     do i=1,ngrd7
        ij=indx7(i)
        fun(ij)=exeven(1)*fun(ij+nni1)+exeven(2)*fun(ij+nni2)+   &
             exeven(3)*fun(ij+nni3)+exeven(4)*fun(ij+nni4)+exeven(5)*fun(ij+nni5)
     enddo
     !       endif
     do i=1,ngrd6a
        ij=indx6a(i)
        fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+  &
             exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
     enddo

     do i=1,ngrd6b
        ij=indx6b(i)
        fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+  &
             exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
     enddo
  else
     !       if (ifill.eq.1.or.ifill.eq.2) then
     do i=1,ngrd7
        ij=indx7(i)
        fun(ij)=0.0_PREC
     enddo
     !       endif

     do i=1,ngrd6a
        ij=indx6a(i)
        fun(ij)=0.0_PREC
     enddo

     do i=1,ngrd6b
        ij=indx6b(i)
        fun(ij)=0.0_PREC
     enddo
  endif

end subroutine sor

