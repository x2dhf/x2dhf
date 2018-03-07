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
! ### vcoul ###
!     Evaluates and returns the boundary value of the Coulomb potential
!     for a given orbital at a given point ($\tilde{V}^a_C$))

module vcoul_m
  implicit none
contains
  function vcoul(iorb,i,j,costh)
    use params
    use discret
    use scf
    use commons8

    implicit none
    integer :: i,iorb,j,kxk,m,n
    real (PREC) :: vcoul
    real (PREC) :: costh,rr,rr1,rr2,pe
    real (PREC), dimension(10) :: dome

    rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
    rr1=1.0_PREC/(rr*r2)

    dome(1)=costh
    dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*5.d-01
    do n=2,mpole-1
       dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
    enddo

    pe=0.0_PREC
    rr2=rr1
    do m=1,mpole
       !         rr2=rr2*rr1
       kxk=iorb+(m-1)*norb
       pe=pe+cmulti(kxk)*dome(m)*(rr1**dble(m+1))
       !         pe=pe+cmulti(kxk)*dome(m)*rr2
    enddo
    vcoul=pe

  end function vcoul
end module vcoul_m
