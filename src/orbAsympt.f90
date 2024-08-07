! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module orbAsympt
  implicit none
contains
  ! ### orbAsymptGet ###
  !
  !   Initializes array edecay which is used by orbAsymptSet to calculate
  !   boundary values of a given orbital at practical infinity.
  !
  subroutine orbAsymptGet (edecay,fa)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,jj
    real (PREC) :: raiq,raiq1
    real (PREC), dimension(*) :: fa
    real (PREC), dimension(nni,4) :: edecay

    jj=0
    do j=mxnmu-3,mxnmu
       jj=jj+1
       do i=1,nni
          raiq1=sqrt(vxisq(j-1)+vetasq(i)-1.0_PREC)
          raiq =sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
          edecay(i,jj)=(abs(fa(i+(j-1)*nni)/fa(i+(j-1)*nni-nni)))**0.250_PREC&
               *exp(sqrt(abs(fa(i+(j-1)*nni-nni)))*(raiq1-raiq))
       enddo
    enddo

  end subroutine orbAsymptGet

  ! ### orbAsymptSet ###
  !
  !     Recalculates asymptotic values of a given orbital at the practical
  !     infinity using exponential decay values prepared by calling orbAsymptDet.
  !
  subroutine orbAsymptSet (psi,edecay)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,jj
    real (PREC), dimension(*) :: psi
    real (PREC), dimension(nni,4) :: edecay

    jj=0
    do j=mxnmu-3,mxnmu
       jj=jj+1
       do i=1,nni
          psi(i+(j-1)*nni)=psi(i+(j-1)*nni-nni)*edecay(i,jj)
       enddo
    enddo

  end subroutine orbAsymptSet

end module orbAsympt
