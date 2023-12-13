! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module initCBAllocArrays
  use params, only : IPREC, PREC
  real (PREC), parameter :: pii = 4.0_PREC*atan(1.0_PREC)
  integer (KIND=IPREC), dimension(4,10) :: ingr1
contains

  ! ****************************************************************************
  ! *                                                                          *
  ! * This routine is part of the GRASP2 package (ver. 1.00, 17 Jan 1992)      *
  ! *                                                                          *
  ! *                                                                          *
  ! *   ----------------   SECTION 01   SUBPROGRAM 14   ----------------       *
  ! *   EVALUATE THE SUM OF THE SERIES                                         *
  ! *                                                                          *
  ! *                       INFINITY      N              K                     *
  ! *              S  (F) =   SUM     (-1)  EXP (N*F) / N                      *
  ! *               K        N = 0                                             *
  ! *                                                                          *
  ! *   FOR K = 2, 3 TO MACHINE PRECISION. THIS IS A UTILITY SUBROUTINE,       *
  ! *   CALLED BY SUBROUTINES NUCPOT AND NCHARG.                               *
  ! *                                                                          *
  ! *   NO SUBROUTINES CALLED                                                  *
  ! *                                                                          *
  ! *   WRITTEN BY FARID A PARPIA, AT OXFORD  LAST REVISION: 25 JAN 1988       *
  ! *                                                                          *
  ! ****************************************************************************

  subroutine es (f,s2f,s3f)
    use params

    implicit none
    integer (KIND=IPREC) :: n
    real (PREC) :: en,enf,f,fase,obn,s2f,s3f,s2last,term2,term3

    n = 0
    s2f = zero
    s3f = zero
    fase = one
1   n = n+1
    en = dble(n)
    obn = one/en
    fase = -fase
    enf = exp (en*f)
    term2 = fase*enf*obn*obn
    term3 = term2*obn
    s2last = s2f
    s2f = s2f+term2
    s3f = s3f+term3
    !if (abs (s2f) .ne. abs (s2last)) goto 1
    if ( abs(s2f-s2last)>epsilon(zero)) goto 1
  end subroutine es
  
  ! ### potscoul2 ###
  ! This function returns a value of the model potential -V0/Sqrt(a^2+r^2)
  ! transformed into the cylindrical coordinates, i.e. the value
  !         -V0/Sqrt(a^2+z^2+s^2)
  ! plus a m-dependent correction resulting from the azimuthal variable
  ! being factored out.
  !
  ! PARAMETERS
  !   s - cylindrical coordinates (z is set to zero)
  !   m - magnetic quantum number
  !  V0 - potential depth
  !   a - potential core width
  function potCoul2(s, m, V0, a, precis)
    implicit none
    real (PREC), intent(in) :: s, V0, a, precis
    integer (KIND=IPREC),intent(in) :: m
    real (PREC) :: popr_cylind, potCoul2

    if (s .lt. precis) then
       popr_cylind=dble(m**2)/(2.0_PREC * precis**2)
    else
       popr_cylind=dble(m**2)/(2.0_PREC * s**2)
    end if
    potCoul2=-V0/sqrt(a**2 + s**2) + popr_cylind;

    return
  end function potCoul2

  ! ### potCoul3 ###
  ! This function returns a value of the model potential -V0/Sqrt(a^2+r^2)
  ! (the smoothed Coulomb poential) transformed into the cylindrical
  ! coordinates, i.e. the value -V0/Sqrt(a^2+z^2+s^2)
  ! plus a m-dependent correction resulting from the azimuthal variable
  ! being factored out.
  !
  ! PARAMETERS
  ! z,s - cylindrical coordinates
  !   m - magnetic quantum number
  !  V0 - potential depth
  !   a - potential core width
  !	 if a=0, V0=1 -> Coulomb potential
  !	 if a>0, V0>0 -> smoothed Coulomb potential
  function potCoul3(z, s, m, V0, a, precis)
    implicit none
    real (PREC), intent(in) :: z, s, V0, a, precis
    integer (KIND=IPREC),intent(in) :: m
    real (PREC) :: popr_cylind, potCoul3

    if (abs(s) .lt. precis .and. abs(z) .lt. precis) then
       popr_cylind = (m**2)/(2.0*precis**2)
       potCoul3 = -V0 / sqrt(a**2 + precis**2) + popr_cylind
    else
       if (s .lt. precis) then
          popr_cylind=(m**2)/(2.0_PREC * precis**2)
       else
          popr_cylind=(m**2)/(2.0_PREC * s**2)
       end if
       potCoul3=-V0/sqrt(a**2 + z**2 + s**2) + popr_cylind;
    end if

    return
  end function potCoul3
  
  ! ### potkh ###
  !
  !     This function returns a value of the Kramers-Henneberger potential
  !     at a point (z,s) in cylindrical coordinates
  !     for the smoothed Coulomb potential plus a m-dependent correction
  !     correction resulting from the azimuthal variable being factored
  !     out.
  !
  !     PARAMETERS
  !      z,s - cylindrical coordinates
  !        m - magnetic quantum number
  !      eps - laser field intensity
  !        w - laser cycle frequency
  !       V0 - original potential depth
  !        a - original potential core width (a>0)
  !        N - number of intervals in the Simpson quadrature
  !
  function potkh(z, s, m, eps, w, V0, a, N, precis)
    implicit none
    real (PREC), intent(in) :: z, s, eps, w, V0, a, precis
    integer (KIND=IPREC),intent(in) :: m, N
    real (PREC) :: popr_cylind, potkh

    if (abs(s) .lt. precis) then
       popr_cylind = (m**2)/(2.0_PREC * precis**2)
       potkh = popr_cylind + AM_HK_ZSimpson(z, precis, 0.0_PREC, eps, w, V0, a, N)
    else
       popr_cylind = (m**2)/(2.0_PREC * s**2);
       potkh = popr_cylind + AM_HK_ZSimpson(z, s, 0.0_PREC, eps, w, V0, a, N)
    end if

    return
  end function potkh

  ! ### am_hk ###
  !
  !     Auxiliary functions
  !     HK in integral form
  !
  function am_hk(t, x, y, z, e0, w, V0, a)
    implicit none
    real (PREC), intent(in) :: t, x, y, z, e0, w, V0, a
    real (PREC) :: alpha0, bx, AM_HK

    alpha0 = e0 / w**2
    bx = x + alpha0 + alpha0 * (cos(w*t)-1.0_PREC)
    AM_HK = -V0 / ( (2.0_PREC * pii / w) * sqrt(a**2 + bx**2 + y**2 + z**2) )

    return
  end function am_hk

  ! ### am_hk_zSimpson ###
  !
  ! Numerical integration by means of the composite Simpson quadrature
  !
  function AM_HK_ZSimpson(x, y, z, e0, w, V0, a, Nt)
    real (PREC), intent(in) :: x, y, z, e0, w, V0, a
    integer (KIND=IPREC),intent(in) :: Nt
    integer (KIND=IPREC) :: i

    ! Interval length
    real (PREC) :: dT
    ! Result
    real (PREC) :: AM_HK_ZSimpson

    ! Calculate interval
    dT = 2*pii/(w*Nt)

    ! Values at the end points
    AM_HK_ZSimpson = AM_HK(0.0_PREC, x, y, z, e0, w, V0, a) &
         + AM_HK(dble(Nt)*dT, x, y, z, e0, w, V0, a)
    ! Values at the intermediate points
    do i=1,Nt-1
       AM_HK_ZSimpson = AM_HK_ZSimpson + 2.0_PREC * AM_HK(dble(i)*dT,x,y,z,e0,w,V0,a)
    end do
    ! Values at the mid points
    do i=0,Nt-1
       AM_HK_ZSimpson = AM_HK_ZSimpson + 4.0_PREC * AM_HK((dble(i)+0.5_PREC)*dT,x,y,z,e0,w,V0,a)
    end do
    ! Normalization
    AM_HK_ZSimpson = AM_HK_ZSimpson * dt / 6.0_PREC

    return
  end function AM_HK_ZSimpson

  ! ### zz1 ###
  !
  !     Evaluate the nuclear potential for Fermi models.
  !
  !     This routine is based on nucpot routine from the GRASP2 package
  !     (ver. 1.00, 1992)
  !
  function zz1(i,j)
    use params
    use discrete
    use commons

    implicit none

    integer (KIND=IPREC) :: i,j
    real (PREC) :: zz1
    real (PREC) :: a,abc,abc2,abc3,at3,atw,c,cba,dmsas,en,facto1,fmtoau,h3,h3php,hpiac2,pi2,&
         rbc,ri,rmc,rmcba,rrms,rrmsfm,s2mcba,s3mcba,s2rcba,s3rcba,sabc3,t,tabc,tfm,thabc2,zbn

    !if (z1.eq.0.0_PREC) then
    if ( abs(z1)<epsilon(zero)) then
       zz1=0.0_PREC
       return
    endif

    fmtoau=1.0e-13_PREC/ainfcm

    ! potential energy is V(i,j)=-zz1(i,j)/r1 -z2/r2
    ! set atomic weight

    atw=z1atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC
    tfm = 2.300_PREC

    ! change units from fm into bohr
    rrms = rrmsfm*fmtoau
    t = tfm*fmtoau
    a = t/(4.00_PREC*log(3.00_PREC))
    facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
    c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)
    abc = a/c
    tabc = two*abc
    abc2 = abc*abc
    thabc2 = three*abc2
    abc3 = abc2*abc
    cba = c/a
    pi2 = pii*pii
    hpiac2 = half*pi2*abc2
    h3 = half*three
    h3php = h3+hpiac2

    call es (-cba,s2mcba,s3mcba)

    sabc3 = six*abc3
    dmsas = -sabc3*s3mcba
    en = one + abc2*pi2 + dmsas
    zbn = z1/en

    ri = r*(vxi(i)+veta(j))/2.0_PREC
    rmc = ri-c
    rmcba = rmc/a
    rbc = ri/c
    if (rbc .le. one) then
       call es (rmcba,s2rcba,s3rcba)
       zz1 = zbn*( dmsas + sabc3*s3rcba+rbc*( h3php-thabc2*s2rcba-half*rbc*rbc) )
    else
       call es (-rmcba,s2rcba,s3rcba)
       zz1 = z1 * ( one+thabc2 * ( rbc *s2rcba+tabc*s3rcba ) / en )
    endif

  end function zz1

  ! ### zz2 ###
  !
  !     Evaluate the nuclear potential for Fermi models.
  !
  !     This routine is based on nucpot routine from the GRASP2 package
  !     (ver. 1.00, 1992)
  !
  function zz2(i,j)
    use params
    use discrete
    use commons
    implicit none
    integer (KIND=IPREC) :: i,j
    real (PREC) :: zz2
    real (PREC) :: a,abc,abc2,abc3,at3,atw,c,cba,dmsas,en,facto1,fmtoau,h3,h3php,hpiac2,pi2,&
         rbc,ri,rmc,rmcba,rrms,rrmsfm,s2mcba,s3mcba,s2rcba,s3rcba,sabc3,t,tabc,tfm,thabc2,zbn

    if (abs(z2)<epsilon(zero)) then
       zz2=0.0_PREC
       return
    endif

    fmtoau=1.0e-13_PREC/ainfcm

    !     Fermi distribution

    !     potential energy is V(i,j)=-zz1(i,j)/r1 -z2/r2
    !     set atomic weight for

    atw=z2atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC
    tfm = 2.300_PREC

    !     change units from fm into bohr

    rrms = rrmsfm*fmtoau
    t = tfm*fmtoau
    a = t/(4.00_PREC*log(3.00_PREC))
    facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
    c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)
    abc = a/c
    tabc = two*abc
    abc2 = abc*abc
    thabc2 = three*abc2
    abc3 = abc2*abc
    cba = c/a
    pi2 = pii*pii
    hpiac2 = half*pi2*abc2
    h3 = half*three
    h3php = h3+hpiac2
    call es (-cba,s2mcba,s3mcba)
    sabc3 = six*abc3
    dmsas = -sabc3*s3mcba
    en = one + abc2*pi2 + dmsas
    zbn = z2/en

    ri = r*(vxi(i)-veta(j))/2.0_PREC
    rmc = ri-c
    rmcba = rmc/a
    rbc = ri/c
    if (rbc .le. one) then
       call es (rmcba,s2rcba,s3rcba)
       zz2 = zbn*( dmsas + sabc3*s3rcba+rbc*( h3php-thabc2*s2rcba-half*rbc*rbc) )
    else
       call es (-rmcba,s2rcba,s3rcba)
       zz2 = z2 * ( one+thabc2 * ( rbc *s2rcba+tabc*s3rcba ) / en )
    endif

  end function zz2

  ! ### zz1g ###
  !
  !     Evaluate the nuclear potential for Gaussian models.
  !
  !     This routine is based on nucpot routine from the GRASP2 package (ver. 1.00, 1992)
  !
  function zz1g(i,j)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j
    real (PREC) :: zz1g
    real (PREC) :: at3,atw,eta1,fmtoau,ri,rrms,rrmsfm

    fmtoau=1.0e-13_PREC/ainfcm

    if (abs(z1)<epsilon(zero)) then
       zz1g=0.0_PREC
       eta1=1.d35
       return
    endif

    !     set atomic weight

    atw=z1atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC

    !     change units from fm into bohr

    rrms = rrmsfm*fmtoau

    !     exponent of the Gaussian distribution $\rho_0 exp(-\eta r^2)$
    !     rrms=(3/(2\eta)^{1/2)

    eta1=3.0_PREC/(2.0_PREC*rrms*rrms)

    ri = r*(vxi(i)+veta(j))/2.0_PREC

    !     set the normalization constant $Z=rho_0 (pi/eta)^{3/2}$

    !     rho0=z1/(pii/eta1)**1.50_PREC
    !     x=eta1*ri*ri
    !     dg15= dgamit(1.50_PREC,x)*x**1.50_PREC*dgamma(1.50_PREC)
    !     zz1g = 2.00_PREC*pii*rho0*(dg15/eta1**1.50_PREC+ri*exp(-x)/eta1)

    !     making use of F77 intrinsic function
    !     erf(x)=2/(Pi)^(1/2)*int_)^x exp(-t^2)dt
    !     if available

    zz1g=z1*erf(sqrt(eta1)*ri)

  end function zz1g

  ! ### zz2g ###
  !
  !     Evaluate the nuclear potential for Gaussian models.
  !
  !     This routine is based on nucpot routine from the GRASP2 package (ver. 1.00, 1992)
  !
  function zz2g(i,j)
    use params
    use discrete
    use commons

    implicit none

    integer (KIND=IPREC) :: i,j
    real (PREC) :: zz2g
    real (PREC) :: at3,atw,eta2,fmtoau,ri,rrms,rrmsfm

    fmtoau=1.0e-13_PREC/ainfcm

    !if (z2.eq.0.0_PREC) then
    if ( abs(z2)>epsilon(zero)) then
       zz2g=0.0_PREC
       eta2=1.e35_PREC
       return
    endif

    ! set atomic weight

    atw=z2atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC

    ! change units from fm into bohr

    rrms = rrmsfm*fmtoau

    ! exponent of the Gaussian distribution $\rho_0 exp(-\eta r^2)$
    ! rrms=(3/(2\eta)^{1/2)

    eta2=3.0_PREC/(2.0_PREC*rrms*rrms)

    ! set the normalization constant $Z=rho_0 (pi/eta)^{3/2}$
    ! rho0=z2/(pii/eta2)**1.50_PREC

    ri = r*(vxi(i)-veta(j))/2.0_PREC

    !  x=eta2*ri*ri
    !  dg15= dgamit(1.50_PREC,x)*x**1.50_PREC*dgamma(1.50_PREC)
    !  zz2g = 2.00_PREC*pii*rho0*(dg15/eta2**1.50_PREC+ri*exp(-x)/eta2)

    ! making use of F77 intrinsic function
    ! erf(x)=2/(Pi)^(1/2)*int_0^x exp(-t^2)dt
    ! if available

    zz2g=z2*erf(sqrt(eta2)*ri)

  end function zz2g
  
  ! ### zgsz1 ###
  !
  !     Evaluate the model HF potential according to the formula derived by
  !     Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z1 centre
  !
  function zgsz1(i,j)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,izz1
    real (PREC) :: zgsz1
    real (PREC) hc,ri

    ri = r*(vxi(i)+veta(j))/two
    izz1=nint(z1)
    if (izz1.eq.0) then
       zgsz1=0.0_PREC
       return
    endif

    if (izz1.eq.1) then
       hc=dgsz(izz1)*(z1)**0.4_PREC
       zgsz1=(z1)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one
       return
    endif

    hc=dgsz(izz1)*(z1-one)**0.4_PREC
    zgsz1=(z1-one)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one

  end function zgsz1

  ! ### zgsz2 ###
  !
  !      Evaluate the model HF potential according to the formula derived by
  !      Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z1 centre
  !
  function zgsz2(i,j)
    use params
    use discrete
    use commons

    implicit none

    integer (KIND=IPREC) :: i,j,izz2
    real (PREC) :: zgsz2
    real (PREC) hc,ri

    ri = r*(vxi(i)-veta(j))/two
    izz2=nint(z2)
    if (izz2.eq.0) then
       zgsz2=0.0_PREC
       return
    endif

    if (izz2.eq.1) then
       hc=dgsz(izz2)*(z2)**0.4_PREC
       zgsz2=(z2)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
       return
    endif

    hc=dgsz(izz2)*(z2-one)**0.4_PREC
    zgsz2=(z2-one)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one

  end function zgsz2
  
  ! ### zgsz1g ###
  !
  !     Evaluate the model HF potential according to the formula derived by
  !     Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z1 centre
  !     with Z1 modified according to the Gauss nucleus model
  !
  function zgsz1g(i,j)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,izz1
    real (PREC) :: zgsz1g
    real (PREC) :: hc,ri,zz1t

    ri = r*(vxi(i)+veta(j))/two
    izz1=nint(z1)
    if (izz1.eq.0) then
       zgsz1g=0.0_PREC
       return
    endif

    if (izz1.eq.1) then
       hc=dgsz(izz1)*(z1)**0.4_PREC
       zgsz1g=(z1)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one
       return
    endif

    zz1t=zz1g(i,j)
    !if (zz1t.eq.zero) then
    if ( abs(zz1t)<epsilon(zero)) then
       zgsz1g=zero
       return
    endif

    ! A modified version of the GSZ model potential is used to make the
    ! comparison between QRHF and X2DHF/DIRAC2D calculations possible  
    
    !if (zz1t<one) then
    !   hc=dgsz(izz1)*(zz1t)**0.4_PREC
    !   zgsz1g=(zz1t)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one
    !else
    !   hc=dgsz(izz1)*(zz1t-one)**0.4_PREC
    !   zgsz1g=(zz1t-one)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one
    !endif

    hc=dgsz(izz1)*(zz1t)**0.4_PREC
    zgsz1g=(zz1t)/(hc*(exp(ri/dgsz(izz1))-one)+one)+one

  end function zgsz1g

  ! ### zgsz2g ###
  !
  !     Evaluate the model HF potential according to the formula derived by
  !     Green, Sellin, Zachor (Phys. Rev. 184 (1969) 1) using Z2 centre
  !     with Z2 modified according to the Gauss nucleus model
  !
  function zgsz2g(i,j)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,izz2
    real (PREC) :: zgsz2g
    real (PREC) :: hc,ri,zz2t

    ri = r*(vxi(i)-veta(j))/two
    izz2=nint(z2)
    if (izz2.eq.0) then
       zgsz2g=0.0_PREC
       return
    endif

    if (izz2.eq.1) then
       hc=dgsz(izz2)*(z2)**0.4_PREC
       zgsz2g=(z2)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
       return
    endif

    zz2t=zz2g(i,j)
    !if (zz2t.eq.zero) then
    if ( abs(zz2t)<epsilon(zero)) then
       zgsz2g=zero
       return
    endif

    !if (zz2t<one) then
    !   hc=dgsz(izz2)*(zz2t)**0.4_PREC
    !   zgsz1g=(zz2t)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
    !else
    !   hc=dgsz(izz2)*(zz2t-one)**0.4_PREC
    !   zgsz1g=(zz2t-one)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
    !endif

    hc=dgsz(izz2)*(zz2t)**0.4_PREC
    zgsz2g=(zz2t)/(hc*(exp(ri/dgsz(izz2))-one)+one)+one
    
  end function zgsz2g


  ! ### setOmega ###
  ! Determines optimal omega overrelaxation parameter for relaxing
  ! potentials according to a semiempirical formula due to B.Sobczak's
  ! MSc thesis, Torun 2002.
  function setOmega()
    use params
    use discrete
    use commons

    implicit none
    real (PREC) :: setOmega
    real (PREC) :: a,b,rho

    parameter (a=0.603_PREC,b=0.79_PREC)

    rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
    setOmega=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b

  end function setOmega

  ! ### setOmegaOrb ###
  !
  !     gnuplot fit to omega_orb  (FH)
  !        a =  1.07272         +/- 0.01652      (1.54%)
  !        b = -0.145639        +/- 0.03261      (22.39%)
  
  !     gnuplot fit to omega_orb  (LiH)
  !        a = 1.26598          +/- 0.04347      (3.433%)
  !        b = -0.528979        +/- 0.08577      (16.21%)
  
  function setOmegaOrb()
    use params
    use discrete
    use commons

    implicit none
    real (PREC) :: setOmegaOrb
    real (PREC) :: a,b,rho

    parameter (a=1.073_PREC,b=-0.1456_PREC)
    ! parameter (a=1.266,b=-0.5290)

    rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
    setOmegaOrb=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b
  end function setOmegaOrb

  ! ### setOmegaPot ###  
  !     gnuplot fit to omega_pot  (FH)
  !        a = 0.568667         +/- 0.007234     (1.272%)
  !        b = 0.863499         +/- 0.01428      (1.654%)
  
  !     gnuplot fit to omega_pot  (LiH)
  !        a = 0.52864          +/- 0.005446     (1.03%)
  !        b = 0.942668         +/- 0.01075      (1.14%)
  
  function setOmegaPot()
    use params
    use discrete
    use commons
    
    implicit none
    real (PREC) :: setOmegaPot
    real (PREC) :: a,b,rho
    
    parameter (a=0.5687_PREC,b=0.8635_PREC)
    ! parameter (a=0.5286,b=0.9427)
    
    rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
    setOmegaPot=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b
  end function setOmegaPot

  ! ### setOmegaPot2 ###    
  ! optimal omega for 2nd-order stencil

  function setOmegaPot2()
    use params
    use discrete
    use commons

    implicit none
    real (PREC) :: setOmegaPot2
    real (PREC) :: a,b,rho

    parameter (a=1.0_PREC,b=0.0_PREC)

    rho=0.5_PREC*(cos(pii/dble(mxnmu))+cos(pii/dble(nni)))
    setOmegaPot2=2.0_PREC*a/(1.0_PREC+sqrt(1.0_PREC-rho*rho))+b

  end function setOmegaPot2
  
  ! ### initCBlocks ###
  !
  !     Initializes arrays within common blocks and prepares grids.
  !
  subroutine initCBlocks
    use commons
    use discrete
    use params
    use scfshr
    use solver

    implicit none

    character*8 :: sigma,pi,delta,phi,spac,ger,uger

    integer (KIND=IPREC) :: i,ial,ib,iba,ibo,ie,idi,idu,ipi,ips,ipsu,ipu,isi,isu

    real (PREC) :: heta,hxibeg,hxiend,ovf,rinfig,xmi0
    real (PREC), parameter :: zFactor=0.001_PREC
    data sigma/'sigma'/,pi/'pi'/,delta/'delta'/,phi/'phi'/
    data spac/' '/,ger/'g'/,uger/'u'/

    ! homonuclear case is treated as a heteronuclear one; if ihomon=2, i.e
    ! if HOMO label is used, the g/u symmetry is imposed explicitely

    ! ihomon=0 --  heteronuclear case
    ! ihomon=1 --  homonuclear case |z1-z2|<1.d-6
    ! ihomon=2 --  homonuclear case with forced symmetry

    ! When counting the maximum number of points in mu variable do not forget
    ! that the first point of the next grid is the last one of the previous grid.

    ! Arrays ibmu and iemu contain the first and last addresses of
    ! particular grid within the combined grid

    ibmu(1)=1
    iemu(1)=nmu(1)
    mxnmu=nmu(1)
    ioffs(1)=0
    ! do i=2,ngrids
    !    mxnmu=mxnmu+nmu(i)-1
    !    ibmu(i)=iemu(i-1)
    !    iemu(i)=mxnmu
    !    ioffs(i)=nni*(iemu(i-1)-1)
    ! enddo

    ! determine the size of each grid

    do i=1,ngrids
       ngsize(i)=nni*nmu(i)
    enddo

    ! total no of grid points

    mxsize=nni*mxnmu
    if ( norb.gt.maxorb ) then
       write(*,*) 'Error: number of orbitals cannot exceed',maxorb
       stop 'initCBlocks'
    endif

    nel=0
    do i=1,norb
       nn(i)=mgx(1,i)
       ll(i)=mgx(3,i)
       mm(i)=mgx(3,i)
       iocc(i)=NINT(occ(i)+.000001_PREC)
       nel=nel+iocc(i)

       if (mod(ll(i),2).eq.0) then
          isymOrb(i)= 1
       else
          isymOrb(i)=-1
       endif
    enddo
    nelInput=nel

    if (mm(1)==2.and.mpole<4) then
       mpole=4
       write(iout6,'(/2x,"Warning: mpole set to 4 since delta-type orbital(s) detected."/)')       
    endif

    if (mm(1)==3.and.mpole<=4) then
       mpole=8
       write(iout6,'(/2x,"Warning: mpole set to 8 since phi-type orbital(s) detected."/)')       
    endif

    ! hni - step in ni variable
    ! hmu - step in mu variable

    ! determine step size in ni variable

    hni=pii/dble(nni-1)

    ! determine step size in mu variable for each grid

    ! when ngrids=1 hmu can be calculated from the practical infinity
    ! variable if its value is provided in the input;

    xmi0=2.0_PREC*rgrid(1)/r
    xmi0=log(xmi0+sqrt(xmi0*xmi0-one))
    hmu(1)=xmi0/dble(mxnmu-1)
    rgrid(1)=hmu(1)
    ! initialize mu and xi arryas

    vmu(1) =0.0_PREC
    vxi(1) =cosh(vmu(1))
    vxisq(1)=vxi(1)*vxi(1)
    vxi1(1)=sqrt(vxi(1)*vxi(1)-one)
    vxi2(1)=0.0_PREC

    ib=2

    ie=iemu(1)
    do i=ib,ie
       vmu(i)=vmu(i-1)+hmu(1)
       vxi(i)=cosh(vmu(i))
       vxisq(i)=vxi(i)*vxi(i)
       
       vxi1(i)=sqrt(vxi(i)*vxi(i)-one)
       vxi2(i)=vxi(i)/vxi1(i)
    enddo
    ib=ie+1

    ! initialize ni and eta arrays
    do i=1,nni
       vni(i)=dble((i-1))*hni
       veta(i)=cos(vni(i))
       vetasq(i)=veta(i)*veta(i)

       ! sqrt(1-veta*veta) and  veta/sqrt(1-veta*veta)
       veta1(i)=sqrt(1.0_PREC-veta(i)*veta(i))
       if (veta1(i).lt.precis) then
          veta2(i)=0.0_PREC
       else
          veta2(i)=veta(i)/veta1(i)
       endif
    enddo

    ! determine rinf
#ifdef PRINT
! print=180: subgrid   nni  nmu     hni     hmu      rb    heta     hxib     hxie    
    if (iprint(180).ne.0) then
       write(iout6,*)
       write(iout6,*) 'subgrid   nni  nmu     hni     hmu      rb    heta     hxib     hxie'
    endif
#endif
    
    rinf=r*vxi(iemu(ngrids))/2._PREC
    heta=abs(veta(nni)-veta(nni-1))
    ib=1
    rinfig=r*vxi(iemu(1))/2._PREC
    ie=iemu(1)
    hxibeg=vxi(ib+1)-vxi(ib)
    hxiend=vxi(ie)-vxi(ie-1)

#ifdef PRINT
! print=180: nni,nmu(1),hni,hmu(1),rinfig,heta,hxibeg,hxiend
    if (iprint(180).ne.0) write(iout6,1000) nni,nmu(1),hni,hmu(1),rinfig,heta,hxibeg,hxiend
#endif

#ifdef PRINT
! print=181: vni, vmi
    if (iprint(181).ne.0) then
       write(iout6,'("vni: ")') 
       write(iout6,'(5e16.8)') (vni(i),i=1,nni)
       write(iout6,'("vmu: ")') 
       write(iout6,'(5e16.8)') (vmu(i),i=1,mxnmu)
    endif
#endif
    
    ib=ie

01000 format(4x,i2,3x,i4,2x,i4,2f9.5,f8.3,3f9.5)

    ! ige(i)=1 ungerade
    ! ige(i)=2 gerade or heteronuclear case

    if (lbreakCi) then
       do i=1,norb
          gusym(i)=spac
       enddo
    endif

    ! give a number within a symmetry

    isu=0
    isi=0
    ipu=0
    ipi=0
    idu=0
    idi=0
    ipsu=0
    ips=0

    do ibo=1,norb
       iba=norb+1-ibo
       if (bond(iba).eq.sigma) then
          if (gusym(iba).eq.ger) then
             isi=isi+1
             iorn(iba)=isi
          endif

          if (gusym(iba).eq.spac) then
             isi=isi+1
             iorn(iba)=isi
          endif

          if (gusym(iba).eq.uger) then
             isu=isu+1
             iorn(iba)=isu
          endif
          goto 100
       endif

       if (bond(iba).eq.pi) then
          if (gusym(iba).eq.ger) then
             ipi=ipi+1
             iorn(iba)=ipi
          endif

          if (gusym(iba).eq.spac) then
             ipi=ipi+1
             iorn(iba)=ipi
          endif

          if (gusym(iba).eq.uger) then
             ipu=ipu+1
             iorn(iba)=ipu
          endif
          goto 100
       endif

       if (bond(iba).eq.delta) then
          if (gusym(iba).eq.ger) then
             idi=idi+1
             iorn(iba)=idi
          endif

          if (gusym(iba).eq.spac) then
             idi=idi+1
             iorn(iba)=idi
          endif

          if (gusym(iba).eq.uger) then
             idu=idu+1
             iorn(iba)=idu
          endif
          goto 100
       endif

       if (bond(iba).eq.phi) then
          if (gusym(iba).eq.ger) then
             ips=ips+1
             iorn(iba)=ips
          endif

          if (gusym(iba).eq.spac) then
             ips=ips+1
             iorn(iba)=ips
          endif

          if (gusym(iba).eq.uger) then
             ipsu=ipsu+1
             iorn(iba)=ipsu
          endif
       endif
100    continue
    enddo

    ! ihsym = 1 - g symmetry 
    ! ihsym =-1 - u symmetry
    do i=1,norb
       ihomo(i)=1
       if (gusym(i).eq.'u'.and.orbsym(i).eq.sigma) ihomo(i)=-1
       if (gusym(i).eq.'u'.and.orbsym(i).eq.delta) ihomo(i)=-1
       if (gusym(i).eq.'g'.and.orbsym(i).eq.pi)    ihomo(i)=-1
       if (gusym(i).eq.'g'.and.orbsym(i).eq.phi)   ihomo(i)=-1
    enddo

    do ial=1,norb
       ige(ial)=2
       if (gusym(ial).eq.'u') ige(ial)=1
    enddo

    ! set (more or less) optimal value of omega for orbitals ovforb< 0
    ! by scaling the optimal value for potentials

    !  set optimal value of omega for potentials if ovfcoul(1)< 0
    if (ovfcoul.lt.0.0_PREC) then
       ovfcoul=setOmega()
       ovfexch=ovfcoul
    endif

    ! in case of convergence problems decrease this value explicitely or
    ! change the omega scaling factor accordingly
    if (ovforb.lt.0.0_PREC) then
       ! it turns out that omegasf factor is grid-dependent and the new
       ! formula accounts for that dependence fairly well
       
       !ovforb=omegasf*setOmega()
       ovf=setOmega()
       ovforb=(ovf-one)*ovf

       ! For heavier species (like Kr or Rn) the automatically calculated
       ! overrelaxation parameter for orbital turns out to be a bit too
       ! large. As a result the convergence of the SCF process is not
       ! monotonic. This label can be used to decrease the overrelaxation
       ! factor by Z-dependent term zFactor*max(Z1,Z2).
       
       if (lomegaz) then
          ovforb=ovforb-zFactor*max(z1,z2)
       endif
    endif


    ! omegaopt label
    if (iomega.eq.1) then
       ovforb=omegasfOrb*setOmegaOrb()
       ovfcoul=omegasfPot*setOmegaPot()
       ovfexch=ovfcoul
    endif

    if (iomega.eq.2) then
       ovfcoul=omegasfPot*setOmegaPot2()
       ovforb=omegasfOrb*setOmegaPot2()
       ovfexch=ovfcoul
    endif

    ! initialize array calp of coefficients of associated Legendre
    ! functions needed by mulex and multi routines

    ! data c/1.25d-01,       5.d-01,         5.59016994d-01,
    !        4.33012702d-01, 1.2247448710_PREC, 3.95284708d-01,
    !        1.369306394_PREC,  6.12372436d-01, 1.479019946_PREC,
    !        5.59016994d-01, 5.22912517d-01, 7.07106781d-01/

    calp( 1)=0.125_PREC
    calp( 2)=0.5_PREC
    calp( 3)=sqrt( 5._PREC/ 16._PREC)
    calp( 4)=sqrt( 3._PREC/ 16._PREC)
    calp( 5)=sqrt( 3._PREC/  2._PREC)
    calp( 6)=sqrt( 5._PREC/ 32._PREC)
    calp( 7)=sqrt(15._PREC/  8._PREC)
    calp( 8)=sqrt( 3._PREC/  8._PREC)
    calp( 9)=sqrt(35._PREC/ 16._PREC)
    calp(10)=sqrt( 5._PREC/ 16._PREC)
    calp(11)=sqrt(35._PREC/128._PREC)
    calp(12)=sqrt( 1._PREC/  2._PREC)
  end subroutine initCBlocks

  ! ### initAddr ###
  !
  !     Determines dimensions of arrays and addresses of particular
  !     orbitals and coulomb potentials in cw_orb and cw_coul arrays
  !     addresses of exchange potentials v(i,j) in the cw_exch array.
  !
  !     With each cw_sth array is associated corresponding address array
  !     defining its division into subarrays holding individual
  !     orbitals/potentials:
  !
  !     cw_orb   - i1xx(iorb)
  !     cw_coul  - i2xx(iorb)
  !     cw_exch  - i3xx(k)
  !
  !     cw_suppl - iaddr4(isuppl)
  !     cw_sctch - iaddr5(iscratch)
  !
  !
  !     i1b(iorb)  - the address of the first element of orbital iorb within
  !                  cw_orb array (starting address of the iorb-th orbital)
  !
  !     i1e(iorb)  - the address of the last element of orbital iorb within
  !                  cw_orb array (ending address of orbital iorb)
  !
  !     i1si(iorb) - the size of the cw_orb subarray holding values of the
  !                  iorb-th orbital
  !
  !     i1ng(iorb) - number of grids used by the iorb-th orbital (must be
  !                  1 in x2dhf ver. 2.0)
  !
  !
  !
  !     i2b(iorb)  - the address of the first element of Coulomb potential
  !                  iorb within cw_coul array (starting address of the iorb-th
  !                  potential)
  !
  !     i2e(iorb)  - the address of the last element of Coulomb potential
  !                  iorb within cw_coul array (ending address of the iorb-th
  !                  potential)
  !
  !     i2si(iorb) - the size of the cw_coul subarray holding values of
  !                  the iorb-th Coulomb pptential
  !
  !     i2ng(iorb) - number of grids used by the iorb-th potential (must
  !                  be 1 in x2dhf ver. 2.0)
  !
  !     i3b(k)     - the address of the first element of the exchange potential for
  !                  the k-th pair of orbitals (iorb1,iorb2): k=iorb1+iorb2*(iorb2-1)/2
  !
  !     i3e(k)     - the address of the last element of the exchange potential for
  !                  the k-th pair of orbitals (iorb1,iorb2)
  !
  !     i3si(k)    - the size of the cw_exch subarray holding values of
  !                  exchange potential for the k-th pair of orbitals
  !
  !     i3ng(k)    - number of grids needed to define the k-th exchange potential
  !
  !     ilc(k)     - 0, 1 or 2, the number of exchange potentials for a given
  !                  pair of orbitals

  !     ilc2(iorb1,iorb2)
  !                - 0, 1 or 2, the number of exchange potentials for
  !                  a given pair of orbitals
  !
  !     k2(iorb1,iorb2)
  !                - the index k of exchange potential corresponding to the pair of
  !                  orbitals iorb1 and iorb2
  !
  !     Subarrays of cw_suppl and cw_sctch are defined accordingly by i4xx
  !     and i5xx arrays.
  !
  !     We assume no more than maxorb orbitals and Coulomb potentials and
  !     consequently no more than (maxorb*(maxorb+1)/2) exchange
  !     potentials.
  !
  subroutine initAddr
    use params
    use discrete
    use memory
    use sharedMemory
    use commons
    use solver

    implicit none

    integer (KIND=IPREC) :: i,i3beg,ibeg,iorb,iorb1,iorb2,irec,k,&
         l,l1cur,l2cur,l3cur,l5cur,ngrid,ngds,nons
    integer (KIND=IPREC8) :: iend

    integer (KIND=IPREC) :: idel,in1,in2,ifirst,ipc,nexchpot,nexchpotsTot

    integer (KIND=IPREC) :: ins1c,ins2c,i1bc,ipcsc,ibexcpc,deltam4potc,isymsc,nexchpotsc,maxpotsc
    common /c_interface_1/ ins1c(maxorb,2*maxorb)
    common /c_interface_2/ ins2c(maxorb,2*maxorb)
    common /c_interface_3/ ibexcpc(maxorb,2*maxorb)
    common /c_interface_4/ deltam4potc(maxorb,2*maxorb)
    !common /c_interface_5/ ipcsc(maxorb,2*maxorb)
    common /c_interface_6/ isymsc(maxorb,2*maxorb)
    common /c_interface_7/ i1bc(maxorb)
    common /c_interface_8/ nexchpotsc(maxorb),maxpotsc    
    
    i1b(1)=1
    ngrid=nni*iemu(i1ng(1))
    i1e(1)=ngrid
    i1si(1)=ngrid
    i1mu(1)=iemu(i1ng(1))

    ! current value of cw_orb length is stored in l1cur
    l1cur=ngrid
    do iorb=2,norb
       i1b(iorb)=i1e(iorb-1)+1
       ngrid=nni*iemu(i1ng(iorb))
       i1e(iorb)=i1b(iorb)+ngrid-1
       i1si(iorb)=ngrid
       i1mu(iorb)=iemu(i1ng(iorb))
       l1cur=l1cur+ngrid
    enddo

    ! it is assumed that potential functions are defined on the same grid
    ! as the orbitals 

    i2b(1)=1
    !ngrid=nni*iemu(i2ng(1))
    i2e(1)=ngrid
    i2si(1)=ngrid
    i2mu(1)=iemu(i1ng(1))

    ! current value of cw_coul length is stored in l2cur
    l2cur=ngrid
    do iorb=2,norb
       i2b(iorb)=i2e(iorb-1)+1
       !ngrid=nni*iemu(i2ng(iorb))
       i2e(iorb)=i2b(iorb)+ngrid-1
       i2si(iorb)=ngrid
       i2mu(iorb)=iemu(i2ng(iorb))
       l2cur=l2cur+ngrid
       ! now excp is used for storing also Coulomb potentials
       !i3b(iorb)=i2b(iorb)
    enddo

    ! dimension of i3xx arrays is norb*(norb+1)/2 to allow for two
    ! exchange potentials between non-sigma orbitals

    ! (obsolete) exchange potential corresponding to orbitals iorb1 and iorb2
    ! defined on 1ing(iorb1) and i1ng(iorb2) grids respectively, is
    ! defined on smaller of these grids, i.e
    ! i3ng=min(i1ng(iorb1),i1ng(iorb2))

    ! Since excp is also used for storing Coulomb potentials shift iend
    ! length2.
    if (lcoulexch) then
       iend=length2
    else
       iend=0
    endif

    l=0
    l3cur=0
    do  iorb1=1,norb
       do  iorb2=iorb1,norb
          k=iorb1+iorb2*(iorb2-1)/2
          ilc(k)=0
          ilc2(iorb1,iorb2)=0
          ilc2(iorb2,iorb1)=0          
          i3si(k)=0
          if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) cycle
          ilc(k)=1
          l=l+1
          ilc2(iorb1,iorb2)=1
          ilc2(iorb2,iorb1)=1          

          k2(iorb1,iorb2)=k
          k2(iorb2,iorb1)=k          
          ngds=min(i1ng(iorb1),i1ng(iorb2))
          i3ng(k)=ngds
          ngrid=mxsize
          iend=iend+ngrid
          i3e(k)=iend
          i3b(k)=i3e(k)-ngrid+1
          i3si(k)=ngrid
          i3mu(k)=iemu(i3ng(k))

          i3orb1(k)=iorb1
          i3orb2(k)=iorb2

          l3cur=l3cur+ngrid
          if (iorb1.eq.iorb2) cycle
          if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) cycle
          ilc(k)=2
          ilc2(iorb1,iorb2)=2
          ilc2(iorb2,iorb1)=2          
          l=l+1
          iend=iend+ngrid
          l3cur=l3cur+ngrid

          if (k.gt.maxorb*(maxorb+1)/2) then
             write(*,*) 'initAddr:'
             write(*,*) 'address array for exchange potentials is too short'
             stop 'initAddr'
          endif
       enddo
    enddo

    ! is working array cw_suppl large enough?

    if (nsuppl*mxsize.gt.length4) then
       write(iout6,*) 'cw_suppl is too short!'
       write(iout6,*) 'declared length is:',length4
       write(iout6,*) 'needed length is nsuppl*mxsize'
       stop 'initAddr'
    endif

    do i=1,nsuppl
       i4b(i)=(i-1)*mxsize+1
       i4barr(i)=(i-1)*mxsize+1
       i4e(i)=iorb*mxsize
       i4si(i)=mxsize
       i4ng(i)=ngrids
    enddo

    ! current value of cw_sctch length is stored in l5cur
    mxsize8=(nni+8)*(mxnmu+8)
    do i=1,nsctch
       i5b(i)=(i-1)*mxsize8+1
       i5e(i)=i*mxsize8
       i5si(i)=mxsize8
       i5ng(i)=ngrids
    enddo
    l5cur=i5e(20)

    ! are working arrays large enough?

    if (l1cur.gt.length1) then
       write(iout6,*) 'cw_orb is too short!'
       write(iout6,*) 'declared length is:',length1
       write(iout6,*) 'needed length is:',l1cur
       stop 'initAddr'
    endif

    if (l2cur.gt.length2) then
       write(iout6,*) 'cw_coul is too short!'
       write(iout6,*) 'declared length is:',length2
       write(iout6,*) 'needed length is:',l2cur
       stop 'initAddr'
    endif

    ! l3cur is calculated as if all exchange potential functions were
    ! to be kept in core.
    if (HF) then
       if (l3cur.gt.length3) then
          write(iout6,*) 'cw_exch is too short!'
          write(iout6,*) 'declared length is:',length3
          write(iout6,*) 'needed length is:',l3cur
          stop 'initAddr'
       endif
    endif

    if (l5cur.gt.length5) then
       write(iout6,*) 'cw_sctch is too short!'
       write(iout6,*) 'declared length is:',length5
       write(iout6,*) 'needed length is',nsctch*mxsize8
    endif


01200  format(2i4,i6,3i12)
01210  format(i4,5i12)

    nexchpotsTot=0
    do iorb=1,norb
       ! By default recalculate only those exchange potentials which depend
       ! on orbitals modified so far (note the reverse order of relaxation
       ! in scf)
       ifirst=iorb
       if (.not.lfastexch) then
          ifirst=1
       endif
           
       nexchpot=0
       do iorb1=ifirst,norb
          if (iorb.eq.iorb1.and.mgx(6,iorb).eq.0 ) cycle
          !if ((iorb.eq.iorb1).and.(ilc(iorb*(iorb+1)/2).lt.1)) cycle
          if ((iorb.eq.iorb1).and.(ilc2(iorb,iorb)==0)) cycle          
          nexchpot=nexchpot+1
          
          if (iorb.lt.iorb1) then
             in1=iorb
             in2=iorb1
          else
             in1=iorb1
             in2=iorb
          endif
       
          ins1(iorb,nexchpot)=in1
          ins2(iorb,nexchpot)=in2

         
          idel=abs(mgx(6,in1)-mgx(6,in2))
          if (in1.eq.in2) idel=2*mgx(6,in1)
          
          deltam4pot(iorb,nexchpot)=idel
          
          ipc=in1+in2*(in2-1)/2
          ipcs(iorb,nexchpot)=ipc
          
          if (mod(idel,itwo).eq.0) then
             isyms(iorb,nexchpot)= 1
          else
             isyms(iorb,nexchpot)=-1
          endif
          
          ibexcp(iorb,nexchpot)=i3b(ipc)
          
          !if (ilc(ipc).eq.2) then
          if (ilc2(in1,in2)==2) then          
             nexchpot=nexchpot+1
             deltam4pot(iorb,nexchpot)=mgx(6,in1)+mgx(6,in2)
             
             ins1(iorb,nexchpot)=in1
             ins2(iorb,nexchpot)=in2
             
             ipcs(iorb,nexchpot)=in1+in2*(in2-1)/2
             ibexcp(iorb,nexchpot)=i3b(ipc)+mxsize
             idel=deltam4pot(iorb,nexchpot)
             if (mod(idel,itwo).eq.0) then
                isyms(iorb,nexchpot)= 1
             else
                isyms(iorb,nexchpot)=-1
             endif
          endif
       enddo
       nexchpots(iorb)=nexchpot
       if (nexchpot>maxpots) maxpots=nexchpot
       if (maxpots>nthreads4coulexch) nthreads4coulexch=maxpots
       nexchpotsTot=nexchpotsTot+nexchpot
#ifdef PRINT
! print= 10: iorb,nexchpots(iorb)       
       if (iprint(10).ne.0) write(*,'(/2x,"iorb nexchpots(iorb)",2i5)') iorb,nexchpots(iorb)
#endif
    enddo
    nexch=nexchpotsTot
    
#ifdef PRINT
! print= 11: nexchpots/norb
    if (iprint(11).ne.0) write(*,'(/2x,"nexchpots/norb=",f6.1)') &
         dble(nexchpotsTot)/dble(norb)
#endif
    
    if (lcoulexch) then
#ifdef PRINT
! print= 12: lcoulexch
       if (iprint(12).ne.0) write(*,'(/2x,"lcoulexch=",l1)') lcoulexch
#endif
       do iorb=1,norb
          if (DFT.or.HFS.or.SCMC.or.OED.or.TED) nexchpots(iorb)=0
          nexchpot=nexchpots(iorb)+1
          nexchpots(iorb)=nexchpot
          in1=iorb
          in2=iorb
          
          ins1(iorb,nexchpot)=in1
          ins2(iorb,nexchpot)=in2
         
          deltam4pot(iorb,nexchpot)=0
          isyms(iorb,nexchpot)=1
          
          ibexcp(iorb,nexchpot)=i2b(iorb)

#ifdef PRINT
! print= 12: iorb,nexchpots(iorb)
          if (iprint(12).ne.0) write(*,'(/2x,"iorb nexchpots(iorb)",2i5)') iorb,nexchpots(iorb)
#endif
       enddo
    endif
    
    maxpots=0    
    do iorb=1,norb
       i1bc(iorb)=i1b(iorb)
       nexchpotsc(iorb)=nexchpots(iorb)
       if (nexchpots(iorb)>maxpots) maxpots=nexchpots(iorb)
       !if (maxpots>nthreads4coulexch) maxpots=nthreads4coulexch

       do nexchpot=1,nexchpots(iorb)
          ins1c(iorb,nexchpot)=ins1(iorb,nexchpot)
          ins2c(iorb,nexchpot)=ins2(iorb,nexchpot)
          ibexcpc(iorb,nexchpot)=ibexcp(iorb,nexchpot)
          deltam4potc(iorb,nexchpot)=deltam4pot(iorb,nexchpot)
          isymsc(iorb,nexchpot)=isyms(iorb,nexchpot)
       enddo
    enddo
    maxpotsc=maxpots
#ifdef PRINT
! print= 10: max number of exchange potentials per orbital
    if (iprint(10).ne.0) write(*,'(/2x, "max number of Coulomb and exchange potentials per orbital=",i2)') maxpots
#endif
  end subroutine initAddr

  ! ### initArrays ###
  !
  !     Setup various arrays and variables
  !
  subroutine initArrays
    use params
    use commons

    ! initialize arrays within common blocks (except for differentiating
    ! arrays which are defined in prepdiff) and check input data
    call initCBlocks

    ! employ input data to determine size and division of cw arrays
    ! initialize address arrays
    call initAddr

    ! Initialises various supplementary arrays of case-dependent lengths
    ! i4b(1)=borb
    ! i4b(2)=bpot
    ! i4b(3)=d
    ! i4b(4)=e
    ! i4b(5)=f0
    ! i4b(6)=f1
    ! i4b(7)=f2
    ! i4b(8)=f3
    ! i4b(9)=f4
    ! i4b(10)=g
    ! i4b(11)=wjac1
    ! i4b(12)=wjac2
    ! i4b(13)=wgt1
    ! i4b(14)=wgt2

    call initSuppl

    ! prepare meshes for all subgrids
    call initMesh

    ! calculate weights of exchange contributions to the total energy
    ! expression for the given open/closed shell scf case

    call initExWeights

  end subroutine initArrays


  ! ### initLP ###
  !
  !     Initialise Legendre polynomials
  !
  subroutine initLP
    use params
    use discrete
    use scfshr
    use commons
    use utils
    use sharedMemory
    
    implicit none
    integer (KIND=IPREC) :: i,ibeg,idel,iorb,iorb1,iorb2,ipc,mu,n,ni
    real (PREC) :: costh,rr,xr,xw
    real (PREC), dimension(10) ::  dome
    real (PREC), dimension(:), pointer :: dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8
    
    dd1 => legendreptr(         1:   mxsize)
    dd2 => legendreptr(  mxsize+1: 2*mxsize)
    dd3 => legendreptr(2*mxsize+1: 3*mxsize)
    dd4 => legendreptr(3*mxsize+1: 4*mxsize)
    dd5 => legendreptr(4*mxsize+1: 5*mxsize)
    dd6 => legendreptr(5*mxsize+1: 6*mxsize)
    dd7 => legendreptr(6*mxsize+1: 7*mxsize)
    dd8 => legendreptr(7*mxsize+1: 8*mxsize)

    do mu=1,mxnmu
       do ni=1,nni
          i=(mu-1)*nni+ni
          
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xr=r2*rr
          
          ! xr=r2*rr  r2=R/2  costh == xi*eta/rr
          ! r==xr=(R/2)rr see eq.13 (CPC 98 (1996) 346)
          
          if (abs(rr).lt.precis) then
             costh=0.0_PREC
          else
             costh=vxi(mu)*veta(ni)/rr
          endif
          
          ! domei=P_k (Legendre polynomial of order k)
          dome(1)=costh
          dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
          do n=2,(mpole-1)
             dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
          enddo
          
          xw=xr
          dd1(i)=dome(1)*xw
          xw=xw*xr
          dd2(i)=dome(2)*xw
          
          if (mpole.ge.3) then
             xw=xw*xr
             dd3(i)=dome(3)*xw
          endif
          
          if (mpole.ge.4) then
             xw=xw*xr
             dd4(i)=dome(4)*xw
          endif
          
          if (mpole.ge.5) then
             xw=xw*xr
             dd5(i)=dome(5)*xw
          endif
          if (mpole.ge.6) then
             xw=xw*xr
             dd6(i)=dome(6)*xw
          endif
          
          if (mpole.ge.7) then
             xw=xw*xr
             dd7(i)=dome(7)*xw
          endif
          
          if (mpole.ge.8) then
             xw=xw*xr
             dd8(i)=dome(8)*xw
          endif
       enddo
    enddo

    
    ! Init some arrays needed for exchMomPT
    nexchmm=0

    do  iorb1=1,norb
       do  iorb2=iorb1,norb
          if (iorb1.eq.iorb2.and.mgx(6,iorb2).eq.0) cycle
          nexchmm=nexchmm+1
          iorb1mm(nexchmm)=iorb1
          iorb2mm(nexchmm)=iorb2

          idel=abs(mgx(6,iorb2)-mgx(6,iorb1))
          if (iorb1.eq.iorb2) idel=2*mgx(6,iorb2)
          idelmm(nexchmm)=idel
          
          ipc=k2(iorb1,iorb2)
          if (ipc.eq.0) cycle
          ipcmm(nexchmm)=ipc
          if (ilc2(iorb1,iorb2)==2) then          
             nexchmm=nexchmm+1
             idel=mgx(6,iorb2)+mgx(6,iorb1)
             idelmm(nexchmm)=idel             
             iorb1mm(nexchmm)=iorb1
             iorb2mm(nexchmm)=iorb2
             ipcmm(nexchmm)=ipc+norb*(norb+1)/2
          endif
          !write(*,'("initLP2:",5i5)') nexchmm,iorb1mm(nexchmm),iorb2mm(nexchmm),idelmm(nexchmm),ipcmm(nexchmm)
       enddo
    enddo
  end subroutine initLP

  ! ### initSupplSM ###
  !
  !     Initialises various supplementary arrays of case-dependent lengths
  !     supported by cw_suppl (one-electron potentials, Jacobians,
  !     integration and differentiation weights, Legendre polynomials, etc)
  !
  subroutine initSuppl
    use params
    use discrete
    use scfshr
    use commons
    use sapPotential
    use blas
    use sharedMemory
    use utils
    implicit none

    integer (KIND=IPREC) :: i,ib,ie,ig,ii,imu,in,izz1,izz2,j,k

    real (PREC) :: atw1,atw2,costh,precis1,r1,rr,rr2,w1,w2,w4,w5,w6,w7,w8,w9,w10,wj1,wj2,wmup,xrr,wktmp,&
         xxplusyy,xxplusyy2,xy2,z,z1t,z2t,zcm

    real (PREC) :: r1t,r2t
    !real (PREC), dimension(nni,mxnmu) :: borb,bpot,d,e,f0,f1,f2,f3,f4,g,wjac1,wjac2
    !real (PREC), dimension(mxsize) :: wgt1,wgt2
    
    real (PREC), pointer :: borb(:,:),bpot(:,:),d(:,:),e(:,:),&
         f0(:,:),f1(:,:),f2(:,:),f3(:,:),f4(:,:),g(:,:),wjac1(:,:),wjac2(:,:)
    
    real (PREC), dimension(9) :: aa1,aa2,a1,a2
    real (PREC), dimension(9,7,9) :: dc1,dc2

    !(borb,bpot,d,e,f0,f1,f2,f3,f4,g,wjac1,wjac2,wgt1,wgt2)
    real (PREC), dimension(:), pointer :: borbpt,bpotpt,dpt,ept,&
         f0pt,f1pt,f2pt,f3pt,f4pt,gpt,wjac1pt,wjac2pt

    real (PREC), dimension(:), pointer :: wgt1pt,wgt2pt

    !     Coefficients of the first and second derivatives taken from from
    !     the 8th-order Stirling interpolation formula
    real (PREC) dmu1c,dmu2c,dnu1c,dnu2c,exevenc,omegac,omega1c
    !common /c_interface_18/ dmu1c(4),dmu2c(4),dnu1c(4),dnu2c(4),exevenc(5),omegac,omega1c

    common /c_interface_40/ dmu1c(4)
    common /c_interface_41/ dmu2c(4)
    common /c_interface_42/ dnu1c(4)
    common /c_interface_43/ dnu2c(4)
    common /c_interface_44/ exevenc(5)
    common /c_interface_45/ omegac,omega1c

    data aa1/ 3.0_PREC, -32.0_PREC, 168.0_PREC, -672.0_PREC, 0.0_PREC, &
         672.0_PREC, -168.0_PREC, 32.0_PREC, -3.0_PREC /
    data aa2/ -9.0_PREC,  128.0_PREC, -1008.0_PREC, 8064.0_PREC,-14350.0_PREC,&
         8064.0_PREC,-1008.0_PREC, 128.0_PREC, -9.0_PREC /

    precis1=precis

    ! iaddr4(1)=borb
    ! iaddr4(2)=bpot
    ! iaddr4(3)=d
    ! iaddr4(4)=e
    ! iaddr4(5)=f0
    ! iaddr4(6)=f1  R^2/2 (cosh^2(mu) - cos^2(nu))=R^2/2 (xi^2 -eta^2)
    ! iaddr4(7)=f2
    ! iaddr4(8)=f3+diag
    ! iaddr4(9)=f4  R/2 cosh(mu)
    ! iaddr4(10)=g
    ! iaddr4(11)=wjac1
    ! iaddr4(12)=wjac2
    ! iaddr4(13)=wgt1
    ! iaddr4(14)=wgt2

    ! wjac1 - Jacobian for Rayleigh quotient
    ! wjac2 - Jacobian for two-electron integrals

    ! initialize arrays needed for the poisson equations

    ! f1 - all potentials must be multiplied by this factor times -1

    ! f4 - array of factors multiplying Coulomb and exchange
    !  potentials (to make the tilda counterparts), see eq.9 and 10

    borbpt =>supplptr(i4b( 1):)
    bpotpt =>supplptr(i4b( 2):)
    dpt    =>supplptr(i4b( 3):)
    ept    =>supplptr(i4b( 4):)                    
    f0pt   =>supplptr(i4b( 5):)
    f1pt   =>supplptr(i4b( 6):)
    f2pt   =>supplptr(i4b( 7):)
    f3pt   =>supplptr(i4b( 8):)
    f4pt   =>supplptr(i4b( 9):)
    gpt    =>supplptr(i4b(10):)
    wjac1pt=>supplptr(i4b(11):)
    wjac2pt=>supplptr(i4b(12):)
    wgt1pt =>supplptr(i4b(13):)        
    wgt2pt =>supplptr(i4b(14):)

    borb (1:nni,1:mxnmu) => borbpt(:)
    bpot (1:nni,1:mxnmu) => bpotpt(:)
    d    (1:nni,1:mxnmu) => dpt(:)        
    e    (1:nni,1:mxnmu) => ept(:)
    f0   (1:nni,1:mxnmu) => f0pt(:)
    f1   (1:nni,1:mxnmu) => f1pt(:)
    f2   (1:nni,1:mxnmu) => f2pt(:)
    f3   (1:nni,1:mxnmu) => f3pt(:)
    f4   (1:nni,1:mxnmu) => f4pt(:)
    g    (1:nni,1:mxnmu) => gpt(:)        
    wjac1(1:nni,1:mxnmu) => wjac1pt(:)
    wjac2(1:nni,1:mxnmu) => wjac2pt(:)                
    
    if (lfefield) then
       izz1=nint(z1)
       izz2=nint(z2)
       atw1=atweight(izz1)
       atw2=atweight(izz2)
       zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC
    else
       zcm=0.0_PREC
    endif

    if (lpotSAP) then
       izz1=nint(z1)
       izz2=nint(z2)
    endif

    wj1 =-r*pii/2.0_PREC
    wj2 = r*r*pii/2.0_PREC
    xrr=half*r
    
    diag=-14350.0_PREC/(5040.0_PREC)*( 1.0_PREC/(hmu(1)*hmu(1))+1.0_PREC/(hni*hni) )
    
    do i=1,mxnmu
       w1=vxi1(i)
       w2=vxi(i)*vxi(i)
       w9=xrr*vxi(i)

       w4=vxi2(i)
       if(i.eq.1) then
          w5=0.0_PREC
          w6=0.0_PREC
       else
          w5= vxi(i)/vxi1(i)-2.0_PREC*vxi1(i)/vxi(i)
          w6= 1.0_PREC/(vxi1(i)*vxi1(i))
       endif
       
       w8=-r*r*r*pii*vxi(i)/2.0_PREC
       
       do j=1,nni
          w10=veta1(j)
          borb(j,i)=w4
          bpot(j,i)=w5
          d(j,i)=veta2(j)
          if (w10.lt.precis) then
             e(j,i)= zero
          else
             e(j,i)=-w6 - one/(w10*w10)
          endif
          w7=w2-veta(j)*veta(j)
          f1(j,i)= r*r*w7/two

          if (lpotCoulomb) then
             ! point nuclei - Coulomb potential
             f0(j,i)=r*(vxi(i)*(z1+z2)+veta(j)*(z2-z1))
          elseif (lpotFermi) then
             ! finite nuclei - Fermi model
             z1t=zz1(i,j)
             z2t=zz2(i,j)
             f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
          elseif (lpotGauss) then
             ! finite nuclei - Gauss model
             z1t=zz1g(i,j)
             z2t=zz2g(i,j)
             f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
          elseif (lpotCoul2) then
             ! OED - 3D smoothed, z-independent Coulomb potential reduced
             ! to 2D external potential
             ! must be multiplied by -f1 
             z=(r/two)*vxi(i)*veta(j)
             xxplusyy2=(r/two)*vxi1(i)*veta1(j)
             f0(j,i)=-f1(j,i)*potCoul2(xxplusyy2,mpot,v0pot,apot,precis1)
          elseif (lpotCoul3) then
             ! OED - 3D smoothed Coulomb potential reduced to 2D external potential
             ! must be multiplied by -f1
             z=(r/two)*vxi(i)*veta(j)
             xxplusyy2=(r/two)*vxi1(i)*veta1(j)
             f0(j,i)=-f1(j,i)*potCoul3(z,xxplusyy2,mpot,v0pot,apot,precis1)
          elseif (lpotKH) then
             ! OED - Kramers-Henneberger potential
             z=(r/two)*vxi(i)*veta(j)
             xxplusyy2=(r/2.0_PREC)*vxi1(i)*veta1(j)
             f0(j,i)=-f1(j,i)*potkh(z,xxplusyy2,mpot,epspot,ompot,v0pot,apot,nsimp,precis1)
          elseif (lpotGSZ) then
             ! OED - Green, Sellin, Zachor model potential
             z1t=zgsz1(i,j)
             z2t=zgsz2(i,j)
             f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
          elseif (lpotGSZG) then
             ! OED - Green, Sellin, Zachor model potential + finite Gauss nucleus model
             z1t=zgsz1g(i,j)
             z2t=zgsz2g(i,j)
             f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
          elseif (lpotSAP) then
             ! OED - SAP due to Susi Lehtola
             r1t=(r/two)*(vxi(i)+veta(j))
             r2t=(r/two)*(vxi(i)-veta(j))
             ! Effective charges Z1(r1) and Z2(r2)
             z1t=effective_charge(izz1,r1t)
             z2t=effective_charge(izz2,r2t)
             f0(j,i)=r*(vxi(i)*(z1t+z2t)+veta(j)*(z2t-z1t))
          elseif (lpotHooke) then
             ! harmonic potential (Hooke's atom, harmonium) located at the geometrical centre
             z=(r/two)*vxi(i)*veta(j)
             xy2=(r*r/4.0_PREC)*(vxisq(i)-one)*(one-vetasq(j))
             f0(j,i) =-half*hooke*hooke*(xy2+z*z)*f1(j,i)
          elseif (lintracule) then
             ! intracular Hamiltonian  E=5/8 
             r1t=(r/two)*(vxi(i)+veta(j))
             f0(j,i) =-one/eight*hooke*hooke*(r1t*r1t)*f1(j,i) - half*r*(vxi(i)-veta(j))
          elseif (lextracule.or.lpotHarm3) then
             ! 3D harmonic potential located at the geometrical centre
             ! T+\frac{1}{2} \omega^2 R^2, R=r_1-r_2 (extracular Hamiltonian) hooke=1/2 E=3/4 (OK)
             z=(r/two)*vxi(i)*veta(j)
             xy2=(r*r/4.0_PREC)*(vxisq(i)-one)*(one-vetasq(j))
             f0(j,i) =-half*hooke*hooke*(xy2+z*z)*f1(j,i) 
          elseif (lpotHarm2) then
             ! harmonic potential 1/2 k^2 (x^2+y^2)
             xy2=(r*r/4.0_PREC)*(vxisq(i)-one)*(one-vetasq(j))
             !f0(j,i) =  -half*hooke*hooke*(xy2)*r*r*(vxisq(i)-vetasq(j))
             f0(j,i) =  -half*hooke*hooke*(xy2)*f1(j,i)
          endif
          
          if (magfield.eq.1) then
             ! external static magnetic field along z-axis for the
             ! first orbital on the input list
             rr=(r/two)*sqrt(vxisq(i)+vetasq(j)-one)
             rr2=rr*rr
             xxplusyy=((r/2.0_PREC)*vxi1(i)*veta1(j))**2
             f0(j,i)=f0(j,i)+(half*gammaf*dble(mgx(3,1))+one/eight*gammaf*gammaf*xxplusyy)*(-one)*f1(j,i)
          endif

          if (lfefield) then
             ! when external electric field is present (z-zcm)*E has to be added (times -f1)
             z=(r/2.0_PREC)*vxi(i)*veta(j)
             if (abs(z).lt.zcutoff) then
                f0(j,i)=f0(j,i)-ffield*(z-zcm)*f1(j,i)
             else
                f0(j,i)=f0(j,i)-ffield*(z-zcm)*exp(zcutoff-abs(z))*f1(j,i)
             endif
          endif
          
          if (iharm2xy.eq.1) then
             ! external three-dimensional harmonic oscilator in x-y-z coordinates
             ! potential must be multiplied by -f1 and added to f0
             imu=i
             in=j
             f0(j,i)=f0(j,i)-f1(in,imu)*harm2xy*harm2xy/2.0_PREC*xrr*xrr*(vxi(i)*vxi(i)+veta(j)*veta(j)-1.0_PREC)
             
          elseif (iharm2xy.eq.2) then
             ! external two-dimensional harmonic oscilator in x-y coordinates
             ! potential must be multiplied by -f1 and added to f0
             imu=i
             in=j
             f0(j,i)=f0(j,i)-harm2xy*harm2xy/2.0_PREC*xrr*xrr*f1(in,imu)*vxi1(imu)*vxi1(imu)*veta1(in)*veta1(in)
          endif
          
          if(i.eq.1) then
             f2(j,i)=0.0_PREC
          else
             f2(j,i)=-r*w7/vxi(i)
          endif
          ! coulSOR and exchSOR routines use f3+diag values
          f3(j,i)=-2.0_PREC/w2+diag 
          g (j,i)= w8*w7
          f4(j,i)= w9
          wjac1(j,i)=wj1*w1*w10
          wjac2(j,i)=wj2*w1*w10*(w1*w1+w10*w10)/vxi(i)
       enddo
    enddo

    ! if (lpotHooke) then
    !    w9=3.0/2.0*(hooke/4.0)**(2.0/3.0)+sqrt(3.0)/4.0*hooke
    !    print *,"hooke energy=",w9,two*w9
    ! endif
    
    ! initialize vector of integration weights
    ! (7-point formula	Abramovic 25.4.16.)

    wmup=0.0_PREC
    do ig=1,ngrids
       ib=ibmu(ig)
       ie=iemu(ig)
       w1 =hmu(ig)/140.0_PREC
       wmu(ib)=wmup+w1*41.0_PREC
       do i=ib+1,ie,6
          wmu(i  )=w1*216.0_PREC
          wmu(i+1)=w1* 27.0_PREC
          wmu(i+2)=w1*272.0_PREC
          wmu(i+3)=w1* 27.0_PREC
          wmu(i+4)=w1*216.0_PREC
          if (i+5.eq.ie) then
             wmu(i+5)=w1*41.0_PREC
          else
             wmu(i+5)=2.0_PREC*w1*41.0_PREC
          endif
       enddo
       wmup=wmu(ie)
    enddo
    
    w1    =hni/140.0_PREC
    wni(1)=w1*41.0_PREC
    do i=2,nni,6
       wni(i  )=w1*216.0_PREC
       wni(i+1)=w1* 27.0_PREC
       wni(i+2)=w1*272.0_PREC
       wni(i+3)=w1* 27.0_PREC
       wni(i+4)=w1*216.0_PREC
       if (i+5.eq.nni) then
          wni(i+5)=w1*41.0_PREC
       else
          wni(i+5)=2.0_PREC*w1*41.0_PREC
       endif
       ! enddo over i
    enddo

    ! multiply Jacobians for one and two-electron integrals by elements
    ! of wmu array
    do i=1,mxnmu
       w1=wmu(i)
       do j=1,nni
          wjac1(j,i)=w1*wjac1(j,i)
          wjac2(j,i)=w1*wjac2(j,i)
       enddo
    enddo

    ! prepare superarrays with the integration weights over ni variable
    ! and with jacobians included
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          wgt1pt(ii+j)=wni(j)*wjac1(j,i)
          wgt2pt(ii+j)=wni(j)*wjac2(j,i)
       enddo
    enddo

    ! Initialize arrays needed for differentiation based on the 9-point
    ! interpolation formula

    ! Initialize dni1, dni2, dmu1 and dmu2

    w1=1.0_PREC/(840.0_PREC)
    w2=1.0_PREC/(5040.0_PREC)

    do  i=1,4
       dmu2(i)=w2*aa2(i)/(hmu(1)*hmu(1))
       dmu1(i)=w1*aa1(i)/hmu(1)
    enddo

    do i=1,4
       dni2(i)=w2*aa2(i)/(hni*hni)
       dni1(i)=w1*aa1(i)/hni
    enddo


    do  i=1,4
       dmu2c(i)=w2*aa2(i)/(hmu(1)*hmu(1))
       dmu1c(i)=w1*aa1(i)/hmu(1)
    enddo

    do i=1,4
       dnu2c(i)=w2*aa2(i)/(hni*hni)
       dnu1c(i)=w1*aa1(i)/hni
    enddo

    
    ! Initialize arrays needed for fast (matrix times matrix) evaluation
    ! of the first and second dervivative terms in the Laplasian (dmu, dnu) and
    ! the first derivatives for gradients (d1mu, d1nu)

    ! derivatives over mu variable
    do k=1,9
       a2(k)=w2*aa2(k)/(hmu(1)*hmu(1))
       a1(k)=w1*aa1(k)/hmu(1)
    enddo

    do i=1,mxnmu
       do k=1,9
          dmu(k,i)=a2(k)+borb(1,i)*a1(k)
          d1mu(k,i)=a1(k)
       enddo
    enddo

    ! derivatives over nu variable
    do k=1,9
       a2(k)=w2*aa2(k)/(hni*hni)
       a1(k)=w1*aa1(k)/hni
    enddo

    do i=1,nni
       do k=1,9
          dni(k,i)=a2(k)+d(i,2)*a1(k)
          d1ni(k,i)=a1(k)
       enddo
    enddo

    ! initialize extrapolation coefficients for even functions
    ! boundary values for odd functions must be all zero

    exeven(1)= 210._PREC/126._PREC
    exeven(2)=-120._PREC/126._PREC
    exeven(3)=  45._PREC/126._PREC
    exeven(4)= -10._PREC/126._PREC
    exeven(5)=   1._PREC/126._PREC

    exevenc(1)= 210._PREC/126._PREC
    exevenc(2)=-120._PREC/126._PREC
    exevenc(3)=  45._PREC/126._PREC
    exevenc(4)= -10._PREC/126._PREC
    exevenc(5)=   1._PREC/126._PREC

    call initLP
    
  end subroutine initSuppl

  subroutine initMesh
    use sharedMemory
    use commons
    use discrete
    use memory
    use params
    use scfshr
    use solver
    use sharedMemory

    implicit none
    integer (KIND=IPREC) :: ib1,ib2,ib3,ib4,ib5,ibeg1,ibeg2,ibeg3,ibeg4,ibeg5, &
         igs,imut,init,is34,is5,mxnmu16,mxsmu,ngrid16

    integer (KIND=IPREC) :: error,i,icol,j,threadID,startRange,step,stopRange
    integer (KIND=IPREC),dimension(:), allocatable :: itemp1,itemp2
    !integer (KIND=IPREC),dimension(length6) :: cw_sor
    integer (KIND=IPREC),dimension(:), pointer :: cw_sor

    integer (KIND=IPREC),parameter :: ig=1
    integer (KIND=IPREC),parameter :: imcase=1

    integer (KIND=IPREC) :: mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c
    
    integer (KIND=IPREC) :: nstartc,nstopc    
    common /c_interface_9/ nstartc(0:4,0:maxthreads4mcsor-1)
    common /c_interface_10/ nstopc(0:4,0:maxthreads4mcsor-1)    

    integer (KIND=IPREC) :: nstart6ac,nstop6ac,nstart6bc,nstop6bc,nstart7c,nstop7c        
    real (PREC) omegac,omega1c
    
    common /c_interface_11/ nstart6ac(0:maxthreads4mcsor-1)
    common /c_interface_12/ nstop6ac(0:maxthreads4mcsor-1)    
    common /c_interface_13/ nstart6bc(0:maxthreads4mcsor-1)
    common /c_interface_14/ nstop6bc(0:maxthreads4mcsor-1)    
    common /c_interface_15/ nstart7c(0:maxthreads4mcsor-1)
    common /c_interface_16/ nstop7c(0:maxthreads4mcsor-1)    
    common /c_interface_17/ mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c
    common /c_interface_45/ omegac,omega1c

    !  integer (KIND=IPREC),dimension(*) :: cw_sor

    ! Array cw_sor is used to store index arrays needed to perform
    ! sor relaxation on normal and extended grids. Two cases have to be
    ! considered since all ngrids-1 subgrids can either be interior or
    ! last ones depending on the particular orbital or potential
    ! function relaxed (outer diffused orbitals would need, say, a grid
    ! consisting of the 3 subgrids while for the proper description of the
    ! the core, very much confined orbitals only one subgrid would sufice).

    ! Division of the cw_sor array is as folows:

    ! ind1   - relaxation ordering on the extended grid
    ! ind2   - relaxation ordering on the normal grid
    ! ind1ex - extrapolation ordering on the extended grid
    ! ind2ex - extrapolation ordering on the normal grid
    ! ind3   - number of points for relaxation and extrapolation
    ! mxsize8==(nni+8)*(mxnmu+8)
    ! mxnmu8=mxnmu+8

    ! case 1 - all ngrids subgrids are considered as the last ones;
    !          in effect the last 4 columns of orbitals and potentials
    !      are not relaxed two arrays of the indx1 type

    ! case 2 - all (ngrids-1) subgrids are considered as the interior ones

    ! cw_sor consists of
    ! - ngrids     arrays of indx1 type (case 1);
    ! - (ngrids-1) arrays of indx1 type (case 2);
    !              the length of each array is determined by ngsize
    !
    ! - ngrids     arrays of indx2 type (case 1);
    ! - (ngrids-1) arrays of indx2 type (case 2);
    !          the length of each array is determined by ngsize
    !
    ! - ngrids	   arrays of indx6a (case 1)
    ! - (ngrids-1) arrays of indx6a (case 2)
    ! - ngrids	   arrays of indx6b (case 1)
    ! - (ngrids-1) arrays of indx6b (case 2)
    ! - ngrids	   arrays of indx7 (case 1)
    ! - (ngrids-1) arrays of indx7 (case 2)
    !          the length of each array is determined by nmu

    ! ingr1(k,ig), ig=1,ngrids, k=1,..,4, contain number of
    ! grids points to be relaxed in each subgrid (k=1) and number of points
    ! for which extrapolation is to be carried out (k=2 - ngrid6a, k=3 - ngrid6b,
    ! k=4 - ngrid7); ingr1 corresponds to the extended (imcase=2) and normal
    ! case (imcase=1), respectively.

    ! i6b(1) - begining of ind1
    ! i6b(2) - begining of ind2
    ! i6b(3) - begining of index1 (indx6a)
    ! i6b(4) - begining of index2 (indx6b)
    ! i6b(5) - begining of index3 (indx7)

    ! more space is reserved for the ordering arrays to avoid problems when
    ! more than 3 subgrid would be used in the future

    ! determine the largest subgrid length in mu variable

    cw_sor=>sorptr(1:)
    
    mxsmu=0
    if (mxsmu.lt.nmu(ig)) mxsmu=nmu(ig)

    is34=ngrids*mxsmu
    is5 =ngrids*nni
    mxnmu16=mxnmu+16
    ngrid16=nni*mxnmu16

    i6b(1)=1
    i6b(2)=i6b(1)+2*ngrid16
    i6b(3)=i6b(2)+2*ngrid16
    !i6b(4)=i6b(3)+2*is34
    !i6b(5)=i6b(4)+2*is5

    i6b(4)=i6b(3)+is34
    i6b(5)=i6b(4)+is34

    
    ! i6b(4)=i6b(3)+is34
    ! i6b(5)=i6b(4)+is5

    ! determine ordering of grid points assuming the subgrids are the last
    ! ones (grid points in the last 4 columns are excluded from relaxation
    ! process

    ibeg1=i6b(1)
    ibeg2=i6b(2)
    ibeg3=i6b(3)
    ibeg4=i6b(4)
    ibeg5=i6b(5)
    imut=0
    init=0
    igs=0

    ib1=ibeg1+igs
    ib2=ibeg2+igs
    ib3=ibeg3+imut
    ib4=ibeg4+imut
    ib5=ibeg5+init
    iadext(ig)=ib1
    iadnor(ig)=ib2
    iadex1(ig)=ib3
    iadex2(ig)=ib4
    iadex3(ig)=ib5

    iadextc=ib1
    iadnorc=ib2
    iadex1c=ib3
    iadex2c=ib4
    iadex3c=ib5
    
    ! call mesh(imcase,ig,ib1,ib2,ib3,ib4,ib5,cw_sor(ib1),cw_sor(ib2),cw_sor(ib3),cw_sor(ib4),cw_sor(ib5))
    call meshsize(imcase,ig)
    !call mesh(imcase,ig,cw_sor(ib1:),cw_sor(ib2:),cw_sor(ib3:),cw_sor(ib4:),cw_sor(ib5:))
    call mesh(imcase,ig)

    imut=imut+nmu(ig)
    init=init+nni
    igs=igs+ngsize(ig)

    ! initialize colours array used in multithreaded version of MCSOR
    ! reorder grid points according to their colours
    ! single grid is assumed

    if (lmcsorpt) then
       allocate(itemp1(ngrid1),stat=error)
       allocate(itemp2(ngrid1),stat=error)       
       !print *,"ngrid1=",ngrid1
       
       do i=1,ngrid1
          itemp1(i)=cw_sor(iadnor(1)+i-1)
          itemp2(i)=cw_sor(iadext(1)+i-1)
       enddo

       j=0
       icolours(1)=1
       do icol=1,5
          icolours(icol)=j+1
          do i=icol,ngrid1,5
             j=j+1
             !print *,'initmesh: icol,j,itemp1',icol,j,itemp1(i),itemp2(i)
             cw_sor(iadnor(1)+j-1)=itemp1(i)
             cw_sor(iadext(1)+j-1)=itemp2(i)
          enddo
       enddo
#ifdef PRINT
       ! print=185: icolours
       if (iprint(185)/=0) write(*,'(" icolours: ",5i10)') (icolours(i),i=1,5)
#endif
       ! calculate ranges of grid points to be assigned to a given thread
       do icol=0,4
          if (icol == 4 ) then
             startRange=icolours(icol+1)
             stopRange=ngrid1
          else 
             startRange=icolours(icol+1)
             stopRange=icolours((icol+1)+1)-1
          endif

          step=(stopRange-startRange+1)/nthreads

          do threadID=0,nthreads-1
             if ( nthreads == 1 ) then
                nstart(icol,threadID)=startRange
                nstop(icol,threadID)=stopRange
             elseif ( (threadID + 1)  == nthreads ) then
                ! last thread
                nstart(icol,threadID)=startRange+threadID*step
                nstop(icol,threadID)=stopRange
             else
                nstart(icol,threadID)=startRange+threadID*step
                nstop(icol,threadID)=nstart(icol,threadID)+step-1
             endif
#ifdef PRINT
! print=186: icol threadID start stop
             if (iprint(186)/=0) write(*,'("icol threadID start stop",8i8)') &
                  icol,threadID,nstart(icol,threadID),nstop(icol,threadID)
#endif
          enddo
       enddo


       ! ngrid6a: calculate ranges of boundary grid points to be assigned to a given thread       
       startRange=1
       stopRange=ngrid6a
       step=(stopRange-startRange+1)/nthreads

       do threadID=0,nthreads-1
          if ( nthreads == 1 ) then
             nstart6a(threadID)=startRange
             nstop6a(threadID)=stopRange
          elseif ( (threadID + 1)  == nthreads ) then
                ! last thread
             nstart6a(threadID)=startRange+threadID*step
             nstop6a(threadID)=stopRange
          else
             nstart6a(threadID)=startRange+threadID*step
             nstop6a(threadID)=nstart6a(threadID)+step-1
          endif
#ifdef PRINT
! print=186: threadID nstart6a nstop6a
          if (iprint(186)/=0) then
             write(*,'("threadID nstart6a nstop6a",3i8)') &
                  threadID,nstart6a(threadID),nstop6a(threadID)
          endif
#endif          
       enddo

       ! ngrid6b: calculate ranges of boundary grid points to be assigned to a given thread       
       startRange=1
       stopRange=ngrid6b
       step=(stopRange-startRange+1)/nthreads

       do threadID=0,nthreads-1
          if ( nthreads == 1 ) then
             nstart6b(threadID)=startRange
             nstop6b(threadID)=stopRange
          elseif ( (threadID + 1)  == nthreads ) then
             ! last thread
             nstart6b(threadID)=startRange+threadID*step
             nstop6b(threadID)=stopRange
          else
             nstart6b(threadID)=startRange+threadID*step
             nstop6b(threadID)=nstart6b(threadID)+step-1
          endif
#ifdef PRINT
! print=186: threadID nstart6b nstop6b
          if (iprint(186)/=0) then
             write(*,'("threadID nstart6b nstop6b",3i8)') &
                  threadID,nstart6b(threadID),nstop6b(threadID)
          endif
#endif
       enddo

       ! ngrid7: calculate ranges of boundary grid points to be assigned to a given thread       
       icol=0
       startRange=1
       stopRange=ngrid7
       step=(stopRange-startRange+1)/nthreads
       
       do threadID=0,nthreads-1
          if ( nthreads == 1 ) then
             nstart7(threadID)=startRange
             nstop7(threadID)=stopRange
          elseif ( (threadID + 1)  == nthreads ) then
             ! last thread
             nstart7(threadID)=startRange+threadID*step
             nstop7(threadID)=stopRange
          else
             nstart7(threadID)=startRange+threadID*step
             nstop7(threadID)=nstart7(threadID)+step-1
          endif
          
#ifdef PRINT
! print=186: threadID nstart7  nstop7
          if (iprint(186)/=0) then
             write(*,'("threadID nstart7  nstop7 ",3i8)') threadID,nstart7(threadID),nstop7(threadID)
          endif
#endif          
       enddo
       
       do threadID=0,nthreads-1
          nstart6ac(threadID)=nstart6a(threadID)
          nstop6ac(threadID)=nstop6a(threadID)
          nstart6bc(threadID)=nstart6b(threadID)
          nstop6bc(threadID)=nstop6b(threadID)
          nstart7c(threadID)=nstart7(threadID)
          nstop7c(threadID)=nstop7(threadID)

          do icol=0,4             
             nstartc(icol,threadID)=nstart(icol,threadID)
             nstopc(icol,threadID)=nstop(icol,threadID)
          enddo
       enddo
    endif

#ifdef PRINT
! print=187: dmu1t,dmu2t,dni1,dni2
    if (iprint(187)/=0) then
       do i=1,4
          write(*,'("dmu1t,dmu2t,dni1,dni2",i5,4e12.4)') i,dmu1(i),dmu2(i),dni1(i),dni2(i)  
       enddo
    endif
#endif
    
    ngrid1 =ingr1(1,1)
    ngrid6a=ingr1(2,1)
    ngrid6b=ingr1(3,1)
    ngrid7 =ingr1(4,1)

    ! Some variables are needed when C versions of (mc)sor routines are used
    ! i.e. when MCSORPT directive is on
    
    nnu1c=nni+8
    nnu2c=2*nnu1c
    nnu3c=3*nnu1c
    nnu4c=4*nnu1c
    nnu5c=5*nnu1c
    mxsizec=mxsize
    mxnmuc=mxnmu

    maxsor1c=maxsororb(1)    
    maxsor2c=maxsororb(2)
    omegac=ovforb
    omega1c=1.0_PREC-omegac

    ngrid6ac=ngrid6a
    ngrid6bc=ngrid6b
    ngrid7c =ngrid7

  end subroutine initMesh

  subroutine mesh(imcase,ig)
    use params
    use discrete
    use commons
    use sharedMemory
    use utils

    implicit none
    integer (KIND=IPREC) :: i,i2,i0,iend1,iend2,ii,imcase,ig,itype,j,j2,jj,nmu8,nmut,nmux,nnix

    ! integer (KIND=IPREC),dimension(ngrid1) :: indx1
    ! integer (KIND=IPREC),dimension(ngrid2) :: indx2
    ! integer (KIND=IPREC),dimension(ngrid6a) :: indx6a
    ! integer (KIND=IPREC),dimension(ngrid6b) :: indx6b
    ! integer (KIND=IPREC),dimension(ngrid7) :: indx7
    
    integer (KIND=IPREC) :: ib1,ib2,ib3,ib4,ib5
    integer (KIND=IPREC),dimension(:), pointer :: indx1
    integer (KIND=IPREC),dimension(:), pointer :: indx2
    integer (KIND=IPREC),dimension(:), pointer :: indx6a
    integer (KIND=IPREC),dimension(:), pointer :: indx6b
    integer (KIND=IPREC),dimension(:), pointer :: indx7


    ib1=iadext(ig)
    ib2=iadnor(ig)
    ib3=iadex1(ig)
    ib4=iadex2(ig)
    ib5=iadex3(ig)
    
    indx1 =>sorptr(ib1:)
    indx2 =>sorptr(ib2:)
    indx6a=>sorptr(ib3:)
    indx6b=>sorptr(ib4:)
    indx7 =>sorptr(ib5:)            
    
    ! initialize variables needed in sor and mcsor

    nni1=nni+8
    nni2=nni1+nni1
    nni3=nni2+nni1
    nni4=nni3+nni1
    nni5=nni4+nni1

    nni8=nni+8
    nmut=nmu(ig)
    nmu8=nmut+8

    if (imcase.eq.1) then
       iend1=nmut
       iend2=nmut-4
    elseif (imcase.eq.2) then
       iend1=nmut+4
       iend2=nmut
    else
       stop "imcase has invalid value in mesh"
    endif

    if (meshOrdering=="col-wise") then    

       ! Natural column-wise ordering, boundary points are excluded,
       ! the upper limit of mu variable is nmu.

       ! ?  values for points (j,nmu+1), (j,nmu+2) and (j,nmu+3) are calculated
       ! ?  from the asymptotic exspansion (the indexes refer to the larger
       ! ?  grid, see routine fill and refill).
       !
       ! region i
       !
       ! in case of more than one grid the upper limit in not nmut-4
       ! but nmut

       ngrid1=0
       do i=6,iend1
          ii=(i-1)*nni8
          do j=6,nni+3
             ngrid1=ngrid1+1
             indx1(ngrid1)=ii+j
          enddo
       enddo

       ngrid2=0
       do i=2,iend2
          ii=(i-1)*nni
          do j=2,nni-1
             ngrid2=ngrid2+1
             indx2(ngrid2)=ii+j
          enddo
       enddo

       goto 1000
    endif

    if (meshOrdering=="middle") then        

       ! The following sequence of mesh points is a modification of
       ! Laaksonen's "middle" type of sweeps

       ! mu changes from mumax/2 to 1 and from mumax/2 to mumax; for
       ! every mu value nu changes from pi/2 to 0 and then from pi/2 to pi

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do i=6,iend1
          if (i.le.nmux) then
             ii=nmux-i+6
          else
             ii=i
          endif
          i2=ii-4

          do j=6,nni+3
             if(j.le.nnix) then
                jj=nnix-j+6
             else
                jj=j
             endif
             j2=jj-4

             i0=i0+1
             indx2(i0)=(i2-1)*nni  +j2
             indx1(i0)=(ii-1)*nni8 +jj
             indx1nu(i0)=j2
             indx1mu(i0)=i2             
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case
       goto 1000
    endif


    if (meshOrdering=="row-wise") then        
       ! natural row-wise ordering
       ! region i
       ngrid1=0
       do j=6,nni+3
          do i=6,iend1
             ii=(i-1)*nni8
             ngrid1=ngrid1+1
             indx1(ngrid1)=ii+j
          enddo
       enddo

       ngrid2=0
       do j=2,nni-1
          do i=2,iend1-4
             ii=(i-1)*nni
             ngrid2=ngrid2+1
             indx2(ngrid2)=ii+j
          enddo
       enddo
       ! boundary points like in the itype=1 case
       goto 1000
    endif

    if (meshOrdering=="rrow-wis") then            

       ! reverse ordering: mi=nmut-4 sweep in ni variable, mi=nmut-5
       ! followed by another sweep in ni etc.. it should allow the
       ! boundary conditions to propagate quicker into the solution.

       ! region i

       ngrid1=0
       do j=nni+3,6,-1
          do i=iend1,6,-1
             ii=(i-1)*nni8
             ngrid1=ngrid1+1
             indx1(ngrid1)=ii+j
          enddo
       enddo

       ngrid2=0
       do j=nni-1,2,-1
          do i=iend1-4,2,-1
             ii=(i-1)*nni
             ngrid2=ngrid2+1
             indx2(ngrid2)=ii+j
          enddo
       enddo
       ! boundary points like in the itype=1 case
       goto 1000
    endif

    if (meshOrdering=="test1") then                
       ! modification of the 2nd ordering (for testing purposes)
       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=6,iend1
             if (i.le.nmux) then
                ii=nmux-i+6
             else
                ii=i
             endif
             i2=ii-4
             i0=i0+1
             indx2(i0)=(i2-1)*nni  +j2
             indx1(i0)=(ii-1)*nni8 +jj
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case

       goto 1000

    endif

    if (meshOrdering=="test2") then                    
       ! another modification of the 2nd ordering (for testing purposes)

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=iend1,6,-1
             ii=i
             i2=ii-4
             i0=i0+1
             indx2(i0)=(i2-1)*nni  +j2
             indx1(i0)=(ii-1)*nni8 +jj
          enddo
       enddo


       ngrid1=i0
       ngrid2=i0
       ! boundary points like in the itype=1 case
       goto 1000
    endif

    if (meshOrdering=="test3") then                        
       ! nu-mu sweeps
       ! nu changes from pi/2 to pi and then from pi/2 to 0; for every
       ! nu value mu changes from mumax/2 to 1 and from mumax/2 to mumax

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=6,iend1
             if (i.le.nmux) then
                ii=nmux-i+6
             else
                ii=i
             endif
             i2=ii-4

             i0=i0+1
             indx2(i0)=(i2-1)*nni  +j2
             indx1(i0)=(ii-1)*nni8 +jj
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case
       goto 1000

    endif

    if (meshOrdering=="test4") then                        
       ! mu changes from mumax to 1; for every mu value nu changes from
       ! pi/2 to pi and then from pi/2 to 0

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do i=iend1,6,-1
          ii=i
          i2=ii-4
          do j=6,nni+3
             if(j.le.nnix) then
                jj=nnix-j+6
             else
                jj=j
             endif
             j2=jj-4

             i0=i0+1
             indx2(i0)=(i2-1)*nni  +j2
             indx1(i0)=(ii-1)*nni8 +jj
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case
       goto 1000
    endif

    if (meshOrdering=="middle2") then
       ! The following sequence of mesh points is a variant of
       ! Laaksonen's "middle" type of sweeps

       ! nu changes from pi/2 to 0 and then from pi/2 to pi.  For every
       ! nu value mu changes from mumax/2 to 1 and from mumax/2 to mumax

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=6,iend1
             if (i.le.nmux) then
                ii=nmux-i+6
             else
                ii=i
             endif
             i2=ii-4

             i0=i0+1
             indx2(i0)=(i2-1)*nni  +j2
             indx1(i0)=(ii-1)*nni8 +jj
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0
       ! boundary points like in the itype=1 case
       goto 1000
    endif

    write(*,*) 'Error: undefined type of ordering'
    stop 'mesh'
    return

01000 continue

    ! region vi
    ! boundary values from extrapolation
    ! values for points (j,nmut-3),...,(j,nmut) are calculated from
    ! asymptotic expansion (orbitals) or multipole expansion (potentials)

    i0=0
    do i=5,iend1
       i0=i0+1
       indx6a(i0)=(i-1)*nni8+5
    enddo
    ngrid6a=i0

    i0=0
    do i=5,iend1
       i0=i0+1
       indx6b(i0)=(i-1)*nni8+nni+4
    enddo
    ngrid6b=i0

    i0=0
    do j=6,nni+3
       i0=i0+1
       indx7(i0)=4*nni8+j
    enddo
    ngrid7=i0

    ingr1(1,ig)=ngrid1
    ingr1(2,ig)=ngrid6a
    ingr1(3,ig)=ngrid6b
    ingr1(4,ig)=ngrid7
 
    ! Orbitals and potentials are stored in one-dimensional arrays. When
    ! they are printed as a two-dimensional ones (\nu=0,\mu_1) element
    ! coresponds to the centre A and (\nu=\pi,\mu_1) -- B.

#ifdef PRINT
! print= 40: ngrid1,ngdr2,ngrid6a,ngrid6b,ngrid7,indx1,indx2,indx6a,indx6b,indx7
    if (iprint(40)/=0) then
       write(*,'(/a8,(10i8))') 'ngrid1 ',ngrid1
       write(*,'(/a8,(10i8))') 'ngrid2 ',ngrid2
       write(*,'(/a8,(10i8))') 'ngrid6a ',ngrid6a
       write(*,'(/a8,(10i8))') 'ngrid6b ',ngrid6b
       write(*,'(/a8,(10i8))') 'ngrid7 ',ngrid7
       write(*,'(/a8/,(10i8))') 'indx1 ',(indx1 (i),i=1,ngrid1)
       write(*,'(/a8/,(10i8))') 'indx2 ',(indx2 (i),i=1,ngrid2)
       write(*,'(/a8/,(10i8))') 'indx6a',(indx6a(i),i=1,ngrid6a)
       write(*,'(/a8/,(10i8))') 'indx6b',(indx6b(i),i=1,ngrid6b)
       write(*,'(/a8/,(10i8))') 'indx7 ',(indx7 (i),i=1,ngrid7 )

    endif
#endif    

  end subroutine mesh

  
  ! ### meshsize ####
  !
  subroutine meshsize(imcase,ig)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,i2,i0,iend1,iend2,ii,imcase,ig,j,j2,jj,nmu8,nmut,nmux,nnix

    ! initialize variables needed in sor and mcsor

    nni1=nni+8
    nni2=nni1+nni1
    nni3=nni2+nni1
    nni4=nni3+nni1
    nni5=nni4+nni1

    nni8=nni+8
    nmut=nmu(ig)
    nmu8=nmut+8

    if (imcase.eq.1) then
       iend1=nmut
       iend2=nmut-4
    elseif (imcase.eq.2) then
       iend1=nmut+4
       iend2=nmut
    else
       stop "imcase has invalid value in meshsize"
    endif

    if (meshOrdering=="col-wise") then        

       ! Natural column-wise ordering, boundary points are excluded,
       ! the upper limit of mu variable is nmu.

       ! ?  values for points (j,nmu+1), (j,nmu+2) and (j,nmu+3) are calculated
       ! ?  from the asymptotic exspansion (the indexes refer to the larger
       ! ?  grid, see routine fill and refill).
       !
       ! region i
       !
       ! in case of more than one grid the upper limit in not nmut-4
       ! but nmut

       ngrid1=0
       do i=6,iend1
          ii=(i-1)*nni8
          do j=6,nni+3
             ngrid1=ngrid1+1
          enddo
       enddo

       ngrid2=0
       do i=2,iend2
          ii=(i-1)*nni
          do j=2,nni-1
             ngrid2=ngrid2+1
          enddo
       enddo

       goto 1000
    endif

    if (meshOrdering=="middle") then            
       ! The following sequence of mesh points is a modification of
       ! Laaksonen's "middle" type of sweeps

       ! mu changes from mumax/2 to 1 and from mumax/2 to mumax; for
       ! every mu value nu changes from pi/2 to 0 and then from pi/2 to pi

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do i=6,iend1
          if (i.le.nmux) then
             ii=nmux-i+6
          else
             ii=i
          endif
          i2=ii-4

          do j=6,nni+3
             if(j.le.nnix) then
                jj=nnix-j+6
             else
                jj=j
             endif
             j2=jj-4

             i0=i0+1
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the col-wise case
       goto 1000
    endif

    if (meshOrdering=="row-wise") then                
       ! natural row-wise ordering

       ngrid1=0
       do j=6,nni+3
          do i=6,iend1
             ii=(i-1)*nni8
             ngrid1=ngrid1+1
          enddo
       enddo

       ngrid2=0
       do j=2,nni-1
          do i=2,iend1-4
             ii=(i-1)*nni
             ngrid2=ngrid2+1
          enddo
       enddo

       ! boundary points like in the col-wise case

       goto 1000

    endif

    if (meshOrdering=="rrow-wis") then                    

       ! reverse ordering: mi=nmut-4 sweep in ni variable, mi=nmut-5
       ! followed by another sweep in ni etc.. it should allow the
       ! boundary conditions to propagate quicker into the solution.

       ! region i

       ngrid1=0
       do j=nni+3,6,-1
          do i=iend1,6,-1
             ii=(i-1)*nni8
             ngrid1=ngrid1+1
          enddo
       enddo

       ngrid2=0
       do j=nni-1,2,-1
          do i=iend1-4,2,-1
             ii=(i-1)*nni
             ngrid2=ngrid2+1
          enddo
       enddo

       ! boundary points like in the col-wise case

       goto 1000
    endif

    if (meshordering=="test1") then
       ! modification of the 2nd ordering (for testing purposes)

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=6,iend1
             if (i.le.nmux) then
                ii=nmux-i+6
             else
                ii=i
             endif
             i2=ii-4
             i0=i0+1
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the col-wise case

       goto 1000

    endif

    if (meshordering=="test2") then    
       ! another modification of the 2nd ordering (for testing purposes)

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=iend1,6,-1
             ii=i
             i2=ii-4
             i0=i0+1
          enddo
       enddo


       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case

       goto 1000
    endif

    if (meshordering=="test3") then    
       ! nu-mu sweeps

       ! nu changes from pi/2 to pi and then from pi/2 to 0; for every
       ! nu value mu changes from mumax/2 to 1 and from mumax/2 to mumax

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=6,iend1
             if (i.le.nmux) then
                ii=nmux-i+6
             else
                ii=i
             endif
             i2=ii-4

             i0=i0+1
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=col-wise case

       goto 1000

    endif

    if (meshOrdering=="test4") then    

       ! mu changes from mumax to 1; for every mu value nu changes from
       ! pi/2 to pi and then from pi/2 to 0

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do i=iend1,6,-1
          ii=i
          i2=ii-4
          do j=6,nni+3
             if(j.le.nnix) then
                jj=nnix-j+6
             else
                jj=j
             endif
             j2=jj-4

             i0=i0+1
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case
       goto 1000
    endif

    if (meshordering=="middle2") then        

       ! The following sequence of mesh points is a variant of
       ! Laaksonen's "middle" type of sweeps

       ! nu changes from pi/2 to 0 and then from pi/2 to pi.  For every
       ! nu value mu changes from mumax/2 to 1 and from mumax/2 to mumax

       i0=0
       nmux=(nmu8-1)/2
       nnix=(nni8-1)/2

       do j=6,nni+3
          if(j.le.nnix) then
             jj=nnix-j+6
          else
             jj=j
          endif
          j2=jj-4

          do i=6,iend1
             if (i.le.nmux) then
                ii=nmux-i+6
             else
                ii=i
             endif
             i2=ii-4

             i0=i0+1
          enddo
       enddo

       ngrid1=i0
       ngrid2=i0

       ! boundary points like in the itype=1 case
       goto 1000
    endif


    write(*,*) 'Error: undefined type of ordering'
    stop 'meshsize'
    return

01000 continue

    ! region vi
    ! boundary values from extrapolation
    ! values for points (j,nmut-3),...,(j,nmut) are calculated from
    ! asymptotic expansion (orbitals) or multipole expansion (potentials)

    i0=0
    do i=5,iend1
       i0=i0+1
    enddo
    ngrid6a=i0

    i0=0
    do i=5,iend1
       i0=i0+1
    enddo
    ngrid6b=i0

    i0=0
    do j=6,nni+3
       i0=i0+1
    enddo
    ngrid7=i0

    ingr1(1,ig)=ngrid1
    ingr1(2,ig)=ngrid6a
    ingr1(3,ig)=ngrid6b
    ingr1(4,ig)=ngrid7
 
    ! Orbitals and potentials are stored in one-dimensional arrays. When
    ! they are printed as a two-dimensional ones (\nu=0,\mu_1) element
    ! coresponds to the centre A and (\nu=\pi,\mu_1) -- B.
#ifdef PRINT
! print= 41: ngrid1,ngrid2,ngrid6a,ngrid6b,ngrid7
    if (iprint(41).ne.0) then
       write(*,'(40a/,10i5)') 'ngrid1,ngrid2,ngrid6a,ngrid6b,ngrid7',  ngrid1,ngrid2,ngrid6a,ngrid6b,ngrid7
    endif
#endif
  end subroutine meshsize
  
  ! ### initExWeights ###
  !
  !     Calculates the weights of exchange contributions to the total
  !     energy expression for an open/closed shell configuration
  !
  subroutine initExWeights
    use params
    use scfshr
    use commons
    implicit none

    integer (KIND=IPREC) :: i, ial, ibe, idiv, ii, imag, imu, isu2, isumi, ix, ixx
    integer (KIND=IPREC) :: j, jj, jx, jxx, k, kk, match, mmi, mpl

    character*8 :: alpha,beta
    data alpha/'+'/, beta/'-'/

    do i=1,2*maxorb*maxorb
       gec(i)=0.0_PREC
    enddo

    do i=1,norb
       do j=1,norb
          offDiagLM=.false.
       enddo
    enddo

    do i=1,4*norb
       spn(i)=spin(i)
    enddo

    do i=1,norb
       isumi=4*(i-1)
       ! gec(i)=occ(i)-1.0
       iloop(i)=2
       iqm(1+isumi)=1
       iqm(2+isumi)=1
       if (mgx(3,i).gt.0) then
          iloop(i)=4
          iqm(3+isumi)=-1
          iqm(4+isumi)=-1
       endif
    enddo

    ! initialize spn and spin arrays
    ! throw out the electrons + = alpha and - = beta
    ! the orbital is locked = 1 (open shell case)
    !                       = 0 (close shell case)

    do i=1,norb
       if (lock(i).eq.0) then
          jj=4*(i-1)
          spn(jj+1)=alpha
          spn(jj+2)=beta
          if (mgx(3,i).gt.0) then
             spn(jj+3)=alpha
             spn(jj+4)=beta
          endif
       endif
    enddo

    do i=1,4*norb
       spin(i)=spn(i)
    enddo

    ! loop now through all electrons

    ! gec-a=gec(kk) gec-b=gec(kk+isu2)
#ifdef PRINT
! print=170: gec-a  gec-b  div
    if (iprint(170).ne.0) then
       write(6,1101)
01101  format(3x,'orb1    ',6x,'orb2       gec-a  gec-b  div')
    endif
#endif
    
    isu2=norb*norb
    do i=1,norb
       ixx=4*(i-1)
       div(i)=zero
       do ii=1,iloop(i)
          ix=ii+ixx
          if ((spn(ix).ne.alpha).and.(spn(ix).ne.beta)) goto 12
          div(i)=div(i)+one

          do j=1,norb
             jxx=4*(j-1)
             kk=i+norb*(j-1)

             mpl=abs(mgx(3,i)-mgx(3,j))
             mmi=abs(mgx(3,i)+mgx(3,j))

             do jj=1,iloop(j)
                jx=jxx+jj
                if ((spn(jx).ne.alpha).and.(spn(jx).ne.beta)) goto 15

                ! div(i)=div(i)+1.00_PREC
                if ((mpl.eq.mmi).and.(spn(ix).eq.spn(jx)).and.(i.ne.j)) then
                   gec(kk)=gec(kk)+1.0_PREC
                endif

                if ((mpl.ne.mmi).and.(spn(ix).eq.spn(jx))) then
                   if ((i.eq.j).and.(iqm(ix).ne.iqm(jx))) then
                      gec(kk)=gec(kk)+1.0_PREC
                   endif

                   if ((i.ne.j).and.(iqm(ix).eq.iqm(jx))) then
                      gec(kk)=gec(kk)+1.0_PREC
                   endif

                   if ((i.ne.j).and.(iqm(ix).ne.iqm(jx))) then
                      gec(kk+isu2)=gec(kk+isu2)+1.0_PREC
                   endif
                endif
15              continue
             enddo
          enddo
12        continue
       enddo

       do idiv=1,norb
          kk=i+norb*(idiv-1)
          if (abs(div(i)).gt.precis) then
             gec(kk)=gec(kk)/div(i)
             gec(kk+isu2)=gec(kk+isu2)/div(i)
          else
             gec(kk)=zero
             gec(kk+isu2)=zero
          endif
#ifdef PRINT
! print=172: gec(kk),gec(kk+isu2),div(i)
          if (iprint(172).ne.0) then
             write(6,1102)  iorn(i),bond(i),gusym(i),iorn(idiv),bond(idiv),gusym(idiv),gec(kk),gec(kk+isu2),div(i)
01102        format(i4,1x,a8,a1,i4,1x,a8,a1,3f7.3)
          endif
#endif
       enddo
    enddo

    ial=0
    ibe=0
    imag=0
    do i=1,norb
       isumi=4*(i-1)
       if (lock(i).eq.0) goto 100

       do ii=1,iloop(i)
          imu=ii+isumi
          if (spn(imu).eq.alpha) ial=ial+1
          if (spn(imu).eq.beta) ibe=ibe+1
          if (spn(imu).eq.alpha.or.spn(imu).eq.beta) then
             imag=imag+iqm(imu)
          endif
       enddo
100    continue
    enddo

    ! One can make the multipliers to appear in Fock equations by
    ! specifying the pair of orbitals in question via largra label.
    ! This array is always set in case an open-shell case is detected.
    
    ! Non-zero entries of nlm array indicate pairs of orbitals having
    ! non-zero off-diagonal Lagrange multipliers. Each pair of orbitals of
    ! the same symmetry is examined to see if an orthogonal transformation
    ! can be applied and off-diagonal Lagrange multipliers set to zero. 

    do j=1,(norb-1)
       jxx=4*(j-1)
       do i=(j+1),norb
          if (mgx(6,j).ne.mgx(6,i)) cycle
          if (gusym(j).ne.gusym(i)) cycle

          if (iloop(j).ne.iloop(i)) cycle
          ixx=4*(i-1)
          match=0
          do ii=1,iloop(i)
             ix=ii+ixx
             jx=ii+jxx
             if (spn(ix).eq.spn(jx)) match=match+1
          enddo

          if (match.ne.iloop(i)) then
             offDiagLM(j,i)=.true.
             offDiagLM(i,j)=.true.
             loffDiagLM=.true.
          endif
       enddo
    enddo

    ! check if off-diagonal LM should be included

    do j=1,(norb-1)
       do i=(j+1),norb
          if (nlm(j,i)==1) then
             offDiagLM(j,i)=.true.
             offDiagLM(i,j)=.true.
          endif
       enddo
    enddo
  end subroutine initExWeights
end module initCBAllocArrays
