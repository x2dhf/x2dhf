module scheduler
  use params

  integer :: ibonus,iepoch,nepochs
  parameter (ibonus=3,iepoch=10,nepochs=5000)
  integer, dimension(maxorb) :: iorbiter
  real (PREC), dimension(nepochs,maxorb) ::  orbenergy,orbnorm
end module scheduler
