! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
!  ### expw86sup ###

!     Calculates exchange energy according to a formula of Parr and Wang Yue
!     Phys. Rev. B 54, 16 533 (1996), Phys. Rev. B 45, 13 244 (1992).

!     J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson
!     and C. Fiolhais, Phys. Rev. B 46, 6671 (1992)

module expw91sup_m
  implicit none
contains
  function expw91sup (wgt2,rho,grho,wk0,wk1)
    use params
    use discret
    use commons8
    use blas_m

    implicit none
    integer :: i
    real (PREC) :: expw91sup
    real (PREC) :: a,a1,a2,a3,a4,ash,b1,const,const23,const43,rho43,s,s2

    real (PREC), dimension(*) :: wgt2,rho,grho,wk0,wk1

    parameter(a=7.795600_PREC,a1=0.1964500_PREC,a2=0.274300_PREC,a3=0.1508400_PREC, &
         a4=100.00_PREC,b1=0.00400_PREC,const23=2.0_PREC/3.0_PREC,const43=4.0_PREC/3.0_PREC)

    ! arcsinh
    ash(s)=log(s+sqrt(one+s*s))

    const=(24.00_PREC*pii*pii)**(-const23)

    !     grho = nabla rho  nabla rho
    !     |nabla rho| = sqrt(grho)

    do i=1,mxsize
       if (abs(rho(i)).lt.precis) then
          wk0(i)=0.0_PREC
       else
          rho43=rho(i)**const43
          s2=const*grho(i)/(rho43*rho43)
          s=sqrt(s2)
          wk0(i)=rho43*(one+a1*s*ash(a*s)+s2*(a2-a3*exp(-a4*s2)))/(one+a1*s*ash(a*s)+b1*s2*s2)
       endif
    enddo

    !     take care of F4 factor
    call multf4(wk0)

    expw91sup=dot(mxsize,wgt2,ione,wk0,ione)

  end function expw91sup
end module expw91sup_m
