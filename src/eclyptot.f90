! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### eclyptot ###

!     Calculates correlation energy according to Lee, Yang, Prr (PRB 37
!     (1988) 785-789)

module eclyptot_m
  implicit none
contains
  function eclyptot (psi,wgt2,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discret
    use commons8
    use blas_m

    implicit none
    real (PREC) :: eclyptot
    integer :: i,iborb,iorb,isiorb
    real (PREC) :: a,b,c,cf,d,const13,const83,const23,const53,const19,const118,ocdown,ocup,t1,t2

    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,  &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    parameter (const13=1.0_PREC/3.0_PREC,const83=8.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC, &
         const53=5.0_PREC/3.0_PREC,const19=1.0_PREC/9.0_PREC,const118=1.0_PREC/18.0_PREC)

    !     coefficients of the Colle-Salvetti formula
    parameter(a=0.049180_PREC,b=0.1320_PREC,c=0.25330_PREC,d=0.3490_PREC)


    cf=0.30_PREC*(three*pii*pii)**const23

    call zeroArray(mxsize,rhotup)
    call zeroArray(mxsize,rhotdown)

    !     calculate total densities for up and down spins
    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 10
       iborb=i1b(iorb)
       isiorb=i1si(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call scal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call scal (isiorb,ocdown,wk2,ione)

       !        store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
10     continue
    enddo

    !     calculate nabla rho nabla rho (up)
    call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhotup)

    !     calculate nabla rho nabla rho (down)
    call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhotdown)


    !     calculate nabla^2 rho (up)
    call n2f (rhotup,wk0,wk1,wk2,rhot)

    !     calculate nabla^2 rho (down)
    call n2f (rhotdown,wk0,wk1,wk2,grhot)

    !     tw (up)
    do i=1,mxsize
       if (abs(rhotup(i)).lt.precis) then
          wk1(i)=0.0_PREC
       else
          wk1(i)=(grhotup(i)/rhotup(i)-rhot(i))/8.0_PREC
       endif
    enddo

    !     tw (down)
    do i=1,mxsize
       if (abs(rhotdown(i)).lt.precis) then
          wk2(i)=0.0_PREC
       else
          wk2(i)=(grhotdown(i)/rhotdown(i)-grhot(i))/8.0_PREC
       endif
    enddo

    !     wk7 - total density, wk6 - tw(total)
    do i=1,mxsize
       wk7(i)=rhotup(i)+rhotdown(i)
       if (abs(wk7(i)).lt.precis) then
          ! FIXME
          wk6(i)=two
       else
          wk6(i)=two*(one-(rhotup(i)*rhotup(i)+rhotdown(i)*rhotdown(i))/(wk7(i)*wk7(i)) )
       endif
    enddo

    !     wk4 - Eq.22 terms from square brackets except for rhot tw
    do i=1,mxsize
       wk4(i)=two**const23*cf*(rhotup(i)**const83+rhotdown(i)**const83)  &
            +const19*(rhotup(i)*wk1(i)+rhotdown(i)*wk2(i))               &
            +const118*(rhotup(i)*rhot(i)+rhotdown(i)*grhot(i))
    enddo

    !     wk7->rhot
    call copy(mxsize,wk7,ione,rhot,ione)

    !     wk4->grhotup
    call copy(mxsize,wk4,ione,grhot,ione)

    !     calculate nabla rho nabla rho (up+down)
    call nfng (rhot,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhotup)

    !     calculate nabla^2 rho (up+down)
    call n2f (rhot,wk0,wk1,wk2,grhotdown)

    !     wk5 - tw
    !     wk6 - gamma (total)
    do i=1,mxsize
       if (rhot(i).lt.precis) then
          ! FIXME
          wk5(i)=0.0_PREC
          wk6(i)=two
       else
          wk5(i)=(grhotup(i)/rhot(i)-grhotdown(i))/8.0_PREC
          wk6(i)=two*(one-(rhotup(i)*rhotup(i)+rhotdown(i)*rhotdown(i))/(rhot(i)*rhot(i)))
       endif
    enddo

    do i=1,mxsize
       if (rhot(i).lt.precis) then
          wk0(i)=0.0_PREC
       else
          !     wk4/grhot +  rhot tw
          t1=grhot(i)-rhot(i)*wk5(i)
          t2=two*b*rhot(i)**(-const53)*t1*exp(-c*rhot(i)**(-const13))
          wk0(i)=-a*(rhot(i)+t2)*wk6(i)/(one+d*rhot(i)**(-const13))
       endif

    enddo

    !     take care of f4 factor
    call multf4(wk0)

    eclyptot=dot(mxsize,wgt2,ione,wk0,ione)

  end function eclyptot
end module eclyptot_m
