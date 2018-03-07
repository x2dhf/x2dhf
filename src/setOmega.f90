! ### setOmega ###
!     Determines optimal omega overrelaxation parameter for relaxing
!     potentials according to a semiempirical formula due to B.Sobczak's
!     MSc thesis, Torun 2002

module setOmega_m
  implicit none
contains
  function setOmega()
    use params
    use discret
    use commons8

    implicit none
    real (PREC) :: setOmega
    real (PREC) :: a,b,rho

    parameter (a=0.603_PREC,b=0.79_PREC)

    rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
    setOmega=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b

  end function setOmega

  !     gnuplot fit to omega_orb  (FH)
  !        a =  1.07272         +/- 0.01652      (1.54%)
  !        b = -0.145639        +/- 0.03261      (22.39%)

  !     gnuplot fit to omega_orb  (LiH)
  !        a = 1.26598          +/- 0.04347      (3.433%)
  !        b = -0.528979        +/- 0.08577      (16.21%)
end module setOmega_m
