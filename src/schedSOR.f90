! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### schedSOR ####

!     Assigns maxsor for each orbital every iepoch iterations

subroutine schedsor
  use params
  use commons8
  use scheduler

  implicit none
  integer :: i,iorb,inde,indemin,maxlevel,minlevel,nlevels
  integer, dimension(60) :: index
  real (PREC), dimension(60) :: endiff

  if (mod(iscf,iepoch).eq.0) then
     i=iscf/iepoch
     do iorb=1,norb
        orbenergy(i,iorb)=eng(iorb)
        orbnorm(i,iorb)=area(iorb)
     enddo
     
     if (i.eq.1) then
        do iorb=1,norb
           iorbiter(iorb)=iepoch*maxsororb(iorb)
        enddo
        return
     endif
     
     do iorb=1,norb
        iorbiter(iorb)=iorbiter(iorb)+iepoch*maxsororb(iorb)
     enddo
     
     inde=0
     do iorb=1,norb
        endiff(iorb)=abs(orbenergy(i-1,iorb)-orbenergy(i,iorb))
     enddo
     
     call qsortf(norb,endiff,index)
     
     !      every orbital gets maxsororb number ranging from maxsor2-ibonus
     !      to maxsor2+ibonus according to its position in the index array
     
     inde=index(norb)
     indemin=index(ione)
     
     minlevel=maxsor2-ibonus
     maxlevel=maxsor2+ibonus
     nlevels=maxlevel-minlevel
     
     do i=1,norb
        iorb=index(i)
        maxsororb(iorb)=(norb-i)*minlevel/(norb-1)+(i-1)*maxlevel/(norb-1)
        
        maxsorpot(iorb)=maxsororb(iorb)
        print *, i,iorb,maxsororb(iorb)
     enddo
     
     write(*,'(2i6,5i6)') iscf,i,(iorbiter(iorb),iorb=norb,1,-1)
     
     write(*,'(8x,5e15.5)') (abs(orbenergy(i-1,iorb)-orbenergy(i,iorb)),iorb=norb,1,-1)
     write(*,'(8x,5e15.5)') (abs(orbnorm(i-1,iorb)-orbnorm(i,iorb)),iorb=norb,1,-1)
  endif
  
end subroutine schedsor
