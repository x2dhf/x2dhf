! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

module exchContribs
  implicit none
contains
  ! ### excont ###
  !
  !     Calculates the number of exchange integrals between two given
  !     (open) shells.
  subroutine excont (iorb1,iorb2,ox1,ox2)
    use params
    use commons

    implicit none
    integer (KIND=IPREC) :: i1,i2,i1end,i2end,io1,io2,iorb1,iorb2,ip1,ip2,ipe1,ipe2
    real (PREC) :: ox1,ox2

    character*8 :: alpha,beta

    data alpha,beta /'+','-'/

    ipe1 =mgx(6,iorb1)
    ipe2 =mgx(6,iorb2)

    io1=iorb1
    io2=iorb2

    ox1=0.0_PREC
    ox2=0.0_PREC

    ip1=4*(io1-1)
    ip2=4*(io2-1)

    if (ipe1.eq.0) then
       icase=1
       i1end=2
       if(ipe2.eq.0) then
          i2end=2
       else
          i2end=4
       endif
    else
       i1end=4
       if(ipe2.eq.0) then
          icase=1
          i2end=2
       else
          icase=2
       endif
    endif

    if (icase.eq.1) then

       ! interaction between sigma-sigma or sigma-nonsigma shells

       do i1=1,i1end
          if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 10
          do i2=1,i2end
             if(spin(ip2+i2).ne.alpha.and.spin(ip2+i2).ne.beta) goto 12
             if(spin(ip1+i1).eq.spin(ip2+i2)) ox1=ox1+1.0_PREC
00012        continue
          enddo
00010     continue
       enddo
    elseif(icase.eq.2) then

       ! interaction between nonsigma-nonsigma shells
       ! lambda positive for both orbitals

       do i1=1,2
          if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 20
          do i2=1,2
             if(spin(ip2+i2).ne.alpha.and.spin(ip2+i2).ne.beta) goto 22
             if(spin(ip1+i1).eq.spin(ip2+i2)) ox1=ox1+1.0_PREC
00022        continue
          enddo
00020     continue
       enddo

       ! lambda negative for both orbitals

       do i1=3,4
          if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 30
          do i2=3,4
             if(spin(ip2+i2).ne.alpha.and.spin(ip2+i2).ne.beta) goto 32
             if(spin(ip1+i1).eq.spin(ip2+i2)) ox1=ox1+1.0_PREC
00032        continue
          enddo
00030     continue
       enddo

       ! lambda positive and negative

       do i1=1,2
          if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 40
          do i2=3,4
             if(spin(ip2+i2).ne.alpha.and.spin(ip2+i2).ne.beta) goto 42
             if(spin(ip1+i1).eq.spin(ip2+i2)) ox2=ox2+1.0_PREC
00042        continue
          enddo
00040     continue
       enddo

       ! lambda negative positive

       do i1=3,4
          if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 50
          do i2=1,2
             if(spin(ip2+i2).ne.alpha.and.spin(ip2+i2).ne.beta) goto 52
             if(spin(ip1+i1).eq.spin(ip2+i2)) ox2=ox2+1.0_PREC
00052        continue
          enddo
00050     continue
       enddo
    endif
  end subroutine excont

  ! ### exint ###
  !
  !     Calculates the number of exchange integrals within one open nonsigma
  !     shell.
  !
  subroutine exint (iorb1,ox1)
    use params
    use commons

    implicit none
    integer (KIND=IPREC) :: i1,i2,iorb1,ip1
    real (PREC) :: ox1

    character*8 :: alpha,beta

    data alpha,beta /'+','-'/

    ox1=0.0_PREC
    ip1=4*(iorb1-1)

    ! interaction within nonsigma shell

    do i1=1,2
       if(spin(ip1+i1).ne.alpha.and.spin(ip1+i1).ne.beta) goto 10
       do i2=3,4
          if(spin(ip1+i2).ne.alpha.and.spin(ip1+i2).ne.beta) goto 12
          if(spin(ip1+i1).eq.spin(ip1+i2)) ox1=ox1+1.0_PREC
00012     continue
       enddo
00010  continue
    enddo
  end subroutine exint

end module exchContribs
