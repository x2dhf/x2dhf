! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fdftpot ###

!     Calculates the coefficients appearing in the LDA potential formulea.
module fdftpot_m
  implicit none
contains
  function fdftpot (alpha)
    use params
    use commons8

    implicit none
    real (PREC) :: fdftpot
    real (PREC) :: alpha,const13,const23,const34

    parameter (const13=1.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC,const34=3.0_PREC/4.0_PREC)

    if (idftex.eq.1) then

       !         fldapot=-alphaf*three/two*(three/pii)**const13

       !        In case of the spin local density approximation where
       !        \rho_{total}^{1/3}=\sum_{spinorbital i}
       !                           [rho_{\alpha}^{1/3}(i) + \rho_{\beta}^{1/3}(i)]
       !        an additional factor 2**(-2/3) shows up

       fdftpot=-alpha*three/two*(three/pii)**const13/two**(const23)

       ! according to ACES
       !         fdftpot=-alpha*three/two*(three/pii)**const13*two**(const13)

    elseif (idftex.eq.2) then
       !        D.Becke Phys. Rev. A 38 (1988) 3098
       fdftpot=-alpha*three/two*(three/pii)**const13/two**(const23)

    elseif (idftex.eq.3) then
       !        generalized gradient approximation (GGA)
       !        Perdew & Wang Phys. Rev. B 33 (1986) 8800
       fdftpot=-const34*(three/pii)**const13

    else
       stop "Invalid argument idftex to fdftpot."

    endif
  end function fdftpot
end module fdftpot_m
