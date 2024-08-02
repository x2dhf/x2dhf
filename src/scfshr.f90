! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2024  Jacek Kobus 

module scfshr
  use params
  implicit none
  real (PREC) :: diag
  real (PREC), dimension(1,1) :: h
  real (PREC), dimension(maxorb+(maxmpole-1)*maxorb) :: cmulti
  real (PREC), dimension(maxorb*(maxorb+1)/2) :: pnx
  
  real (PREC), dimension(maxorb*(maxorb+1)) :: excdi,excqu,excoc,exche,exc5,exc6,exc7,exc8
  real (PREC), dimension(2*maxorb*maxorb) :: gec,gca,gxa,gcb,gxb,gc,gx
  real (PREC) :: ovforb,ovfcoul,ovfexch  
  real (PREC), dimension(maxorb,maxorb) :: sflagrat
end module scfshr
