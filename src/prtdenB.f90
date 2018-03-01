! ### prtdenB ###

subroutine prtdenB (m,n,a,ioutmat)
  use params
  use discret
  use commons8

  implicit none
  integer :: in,imu,ioutmat,m,n
  real (PREC) :: r2t
  real (PREC), dimension(m,n) :: a

  write(ioutmat,'(10x,"r(au)",12x,"total electronic density")')
  in=1
  do imu=1,n
     r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))
     write(ioutmat,'(2E25.16)') r2t,a(in,imu)
  enddo
  
end subroutine prtdenB
