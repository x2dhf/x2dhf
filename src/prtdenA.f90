! ### prtdenA ###

subroutine prtdenA (m,n,a,ioutmat)
  use params
  use discret
  use commons8

  implicit none
  integer :: in,imu,ioutmat,m,n
  real (PREC) :: r1t
  real (PREC), dimension(m,n) :: a
      
  write(ioutmat,'(10x,"r(au)",12x,"total electronic density")')
  in=m
  do imu=1,n
     r1t=(half*r)*(vxi(imu)+veta(in))
     write(ioutmat,'(2E25.16)') r1t,a(in,imu)
  enddo
  
end subroutine prtdenA
