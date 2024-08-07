module params
  integer, parameter :: IPREC = 4
  integer, parameter :: IPREC8 = 8

  integer, parameter :: PREC = kind(1.0d0)
  integer, parameter :: PREC16 = kind(1.0d0)
  ! integer (KIND=IPREC),parameter :: PREC16 = SELECTED_REAL_KIND (1.d308 )
  
  ! In order to make the 2DHF program the array dimensions in various routines have to be
  ! set according to the values of the following variables:

  ! maxnu    = maximum number of grid points in \nu direction
  ! maxmu    = maximum number of grid points in \mu direction

  ! maxbasis = maximum number of primitive gaussian basis functions used to construct
  !            numerical molecular orbitals out of data retrieved from GAUSSIAN outputs

  character*15, parameter :: version = "x2dhf-v3"
  integer (kind=IPREC), parameter :: maxnu    = 2500
  integer (kind=IPREC), parameter :: maxmu    = 5500
  integer (kind=IPREC), parameter :: maxbasis = 650

  ! The following values should be left unchanged, since the code still has
  ! the values hardcoded in some routines! In case of any changes see also
  ! sorpt.h to adjust maxorb and  max_threads4mcsor variables for C routines as
  ! well.
  integer (kind=IPREC), parameter :: maxgrids = 1
  integer (kind=IPREC), parameter :: maxorb   = 36
  integer (kind=IPREC), parameter :: maxmpole = 8
  integer (kind=IPREC), parameter :: maxthreads4mcsor = 16
  
  character*10 :: formint,formfp,formfp64,formfp128
  data formint /'(10i15)  '/
  data formfp64,formfp128/'(5e25.16) ','(5e41.32) '/

  integer (kind=IPREC), parameter :: izero=0
  integer (kind=IPREC), parameter :: ione=1
  integer (kind=IPREC), parameter :: itwo=2
  integer (kind=IPREC), parameter :: ithree=3
  integer (kind=IPREC), parameter :: ifour=4

  integer (kind=IPREC8), parameter :: i8zero=0
  integer (kind=IPREC8), parameter :: i8one=1
  integer (kind=IPREC8), parameter :: i8two=2
  integer (kind=IPREC8), parameter :: i8three=3
  integer (kind=IPREC8), parameter :: i8four=4

  integer (kind=IPREC), parameter :: iinp5=5
  integer (kind=IPREC), parameter :: iinp11=11
  integer (kind=IPREC), parameter :: iinp12=12
  integer (kind=IPREC) :: iinp13=13   ! This can be reset in the code(!)
  integer (kind=IPREC), parameter :: iinp14=14

  integer (kind=IPREC8), parameter :: i8inp5=5
  integer (kind=IPREC8), parameter :: i8inp11=11
  integer (kind=IPREC8), parameter :: i8inp12=12
  integer (kind=IPREC8) :: i8inp13=13 ! This can be reset in the code(!)
  integer (kind=IPREC8), parameter :: i8inp14=14

  integer (kind=IPREC), parameter :: iout6=6
  integer (kind=IPREC), parameter :: iout21=21
  integer (kind=IPREC), parameter :: iout22=22
  integer (kind=IPREC), parameter :: iout23=23
  integer (kind=IPREC), parameter :: iout24=24

  integer (kind=IPREC8), parameter :: i8out6=6
  integer (kind=IPREC8), parameter :: i8out21=21
  integer (kind=IPREC8), parameter :: i8out22=22
  integer (kind=IPREC8), parameter :: i8out23=23
  integer (kind=IPREC8), parameter :: i8out24=24

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
  !real (PREC), parameter :: pii=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_PREC
  real (PREC), parameter :: pii=3.14159265358979323846264338327950_PREC 
  !real (PREC), parameter :: pii=3.1415926535897932_PREC

  ! Numerical precision
  real (PREC) :: precis

  !integer (KIND=IPREC8) :: lengthfp  
  integer (KIND=IPREC) :: lengthfp
  integer (KIND=IPREC) :: lengthint,lengthintin,lengthfpin

  integer (KIND=IPREC) :: incrni,incrmu,maxunit,iprint16,&
       iord_nu_orb,iord_nu_coul,iord_nu_exch,iord_mu_orb,iord_mu_coul,iord_mu_exch,&
       lexchrecl,nsuppl,nsctch,inon
  integer (KIND=IPREC),dimension(1000) :: idebug,iprint  
 end module params
