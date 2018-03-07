! ### zz2g ###
!
!     Evaluate the nuclear potential for Gaussian models.
!
!     This routine is based on nucpot routine from the GRASP2 package (ver. 1.00, 1992)

module zz2g_m
  implicit none
contains
  function zz2g(i,j)
    use params
    use discret
    use commons8

    implicit none

    integer :: i,j
    real (PREC) :: zz2g
    real (PREC) :: at3,atw,eta2,fmtoau,ri,rrms,rrmsfm

    fmtoau=1.0e-13_PREC/ainfcm

    if (z2.eq.0.0_PREC) then
       zz2g=0.0_PREC
       eta2=1.e35_PREC
       return
    endif

    !     set atomic weight

    atw=z2atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC

    !     change units from fm into bohr

    rrms = rrmsfm*fmtoau

    !     exponent of the Gaussian distribution $\rho_0 exp(-\eta r^2)$
    !     rrms=(3/(2\eta)^{1/2)

    eta2=3.0_PREC/(2.0_PREC*rrms*rrms)

    !     set the normalization constant $Z=rho_0 (pi/eta)^{3/2}$
    !      rho0=z2/(pii/eta2)**1.50_PREC

    ri = r*(vxi(i)-veta(j))/2.0_PREC

    !      x=eta2*ri*ri
    !      dg15= dgamit(1.50_PREC,x)*x**1.50_PREC*dgamma(1.50_PREC)
    !      zz2g = 2.00_PREC*pii*rho0*(dg15/eta2**1.50_PREC+ri*exp(-x)/eta2)

    !     making use of F77 intrinsic function
    !     erf(x)=2/(Pi)^(1/2)*int_0^x exp(-t^2)dt
    !     if available

    zz2g=z2*erf(sqrt(eta2)*ri)

  end function zz2g
end module zz2g_m
