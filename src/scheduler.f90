module scheduler
  use params

  integer :: ibonus,iepoch,nepochs
  parameter (ibonus=3,iepoch=10,nepochs=5000)
  integer, dimension(60) :: iorbiter
  real (PREC), dimension(nepochs,60) ::  orbenergy,orbnorm
end module scheduler
