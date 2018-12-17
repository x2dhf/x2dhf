module params
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0d0)
  !  integer, parameter :: DP16 = SELECTED_REAL_KIND (1.d308 )
  integer, parameter :: PREC = DP
  integer, parameter :: PREC16 = DP

! In order to make the 2DHF program the array dimensions in various routines have to be
! set according to the values of the following variables:

! maxnu    = maximum number of grid points in \nu direction
! maxmu    = maximum number of grid points in \mu direction

! maxbasis = maximum number of primitive gaussian basis functions used to construct
!            numerical molecular orbitals out of data retrieved from GAUSSIAN outputs

  integer, parameter :: maxnu    = 1500
  integer, parameter :: maxmu    = 2500
  integer, parameter :: maxbasis = 650

  ! the following values should be left unchanged, since the code
  ! still has the values hardcoded in some routines!
  integer, parameter :: maxgrids = 1
  integer, parameter :: maxorb   = 60
  integer, parameter :: maxmpole = 8

  character*10 :: formint,formfp,formfp64,formfp128
  data formint /'(10i15)  '/
  data formfp64,formfp128/'(5e25.16) ','(5e41.32) '/

  integer, parameter :: izero=0
  integer, parameter :: ione=1
  integer, parameter :: itwo=2
  integer, parameter :: ithree=3
  integer, parameter :: ifour=4

  integer*8, parameter :: i8zero=0
  integer*8, parameter :: i8one=1
  integer*8, parameter :: i8two=2
  integer*8, parameter :: i8three=3
  integer*8, parameter :: i8four=4

  integer, parameter :: iinp5=5
  integer, parameter :: iinp11=11
  integer, parameter :: iinp12=12
  integer :: iinp13=13   ! This can be reset in the code(!)
  integer, parameter :: iinp14=14

  integer*8, parameter :: i8inp5=5
  integer*8, parameter :: i8inp11=11
  integer*8, parameter :: i8inp12=12
  integer*8 :: i8inp13=13 ! This can be reset in the code(!)
  integer*8, parameter :: i8inp14=14

  integer, parameter :: iout6=6
  integer, parameter :: iout21=21
  integer, parameter :: iout22=22
  integer, parameter :: iout23=23
  integer, parameter :: iout24=24

  integer*8, parameter :: i8out6=6
  integer*8, parameter :: i8out21=21
  integer*8, parameter :: i8out22=22
  integer*8, parameter :: i8out23=23
  integer*8, parameter :: i8out24=24

  real (PREC), parameter :: zero=0.0_PREC
  real (PREC), parameter :: half=0.5_PREC
  real (PREC), parameter :: one=1.0_PREC
  real (PREC), parameter :: two=2.0_PREC
  real (PREC), parameter :: three=3.0_PREC
  real (PREC), parameter :: four=4.0_PREC
  real (PREC), parameter :: five=5.0_PREC
  real (PREC), parameter :: six=6.0_PREC
  real (PREC), parameter :: seven=7.0_PREC
  real (PREC), parameter :: eight=8.0_PREC
  real (PREC), parameter :: nine=9.0_PREC
  real (PREC), parameter :: ten=10.0_PREC

  ! 100 digits from Maple
  real (PREC), parameter :: pii=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_PREC
  !real (PREC), parameter :: pii=3.1415926535897932_PREC

  ! Numerical precision
  real (PREC) :: precis

   integer :: incrni,incrmu,maxunit,lengthfp,lengthint,lengthintin,lengthfpin,iprint16,idat,&
       iord_nu_orb,iord_nu_coul,iord_nu_exch,iord_mu_orb,iord_mu_coul,iord_mu_exch,&
       lexchrecl,nsuppl,nsctch,inon

   integer, dimension(1000) :: idbg,iprint

end module params
