!     gnuplot fit to omega_orb  (FH)
!        a =  1.07272         +/- 0.01652      (1.54%)
!        b = -0.145639        +/- 0.03261      (22.39%)

!     gnuplot fit to omega_orb  (LiH)
!        a = 1.26598          +/- 0.04347      (3.433%)
!        b = -0.528979        +/- 0.08577      (16.21%)

module setOmegaOrb_m
  implicit none
contains
  function setOmegaOrb()
    use params
    use discret
    use commons8

    implicit none
    real (PREC) :: setOmegaOrb
    real (PREC) :: a,b,rho

    parameter (a=1.073_PREC,b=-0.1456_PREC)
    !      parameter (a=1.266,b=-0.5290)

    rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
    setOmegaOrb=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b

  end function setOmegaOrb
end module setOmegaOrb_m
