! ### inCard ###
!
!     This routine reads (and echoes) a single data card and
!
!       - changes upper case letters into a lower case ones
!       - changes tabs (if any) into spaces
!       - removes leading spaces
!
!     Then it scans the line for nonspace fields. The number of fields
!     is stored in jump, the starting point of a given field in istrt(i)
!     and the number of characters in that field in inumb(i). If the
!     routine finds an exclamation mark or a hash then anything on the
!     same line that follows these characters is treated as a comment.

module inCard_m
  implicit none
contains
  subroutine inCard

    use params
    use card

    implicit none

    character*1 :: blnk,exm,hash
    character*1, dimension(80) :: iatmp

    integer :: i,ike,ile,isw,it

    data blnk/' '/,exm/'!'/,hash/'#'/

    do i=1,i80
       ia(i)=' '
    enddo

    !     read a line

    jump = 0
    jrec = 0
    isw = 0
1   read(iinp5,100,err=110,end=40) ia
100 format(80a1)
101 format(2x,80a1)
110 continue

    do i=1,i80
       iatmp(i)=ia(i)
    enddo

    !     set into lower case
    !     FC3, ifort: error while passing i80 to a subroutine
    it=i80
    call lowcase(iatmp,it)

    !     check for tabs, and replace them by blanks
    call tabchk(iatmp,it)

    !     check for blanks at left and remove them, fill with blanks at the end
    call lftpos(iatmp,it)

    do i=1,i80
       ia(i)=iatmp(i)
    enddo


    !     check for an exclamation mark or a hash
    do ike=1,i80
       ile=ike
       if (ia(ile).eq.exm.or.ia(ile).eq.hash) go to 910
    enddo

910 if(ile.eq.1) go to 1
    !     echo the slightly modified input data
    write(iout6,101) ia

    if(ia(ile).eq.exm.or.ia(ile).eq.hash) ile=ile-1

    do i = 1,ile
       if (ia(i).eq.blnk) then
          isw=0
          goto 30
       else
          if (isw.le.0) then
             jump = jump +1
             istrt(jump) = i
             inumb(jump) = 0
             isw=1
             inumb(jump) = inumb(jump) + 1
          else
             inumb(jump) = inumb(jump) + 1
          endif
       endif
30     continue
    enddo

    if (jump.eq.0) goto 1
    return

40  write(iout6,45)
45  format(//1x,'Error: unexpected end of input data. Is STOP label missing?'//)
    stop 'inCard'
  end subroutine inCard


  ! ### inCardh ###
  !
  !     This routine reads a single data card (80 characters) and returns
  !     the data from columns 6-80.

  subroutine inCardh(header80)

    use params
    use card

    implicit none

    character*1 exm,blnk,hash
    character*80 :: header80,h80

    character*1, dimension(80) :: header1,iatmp

    integer :: i,ike,ile,isw,it

    data blnk/' '/,exm/'!'/,hash/'#'/
    equivalence(h80,header1(1))

    do i=1,i80
       ia(i)=' '
    enddo

    !     read a line

    jump = 0
    jrec = 0
    isw = 0

1   read(iinp5,100,err=110,end=40) ia
100 format(80a1)
101 format(2x,80a1)
110 continue

    do i=1,i80
       iatmp(i)=ia(i)
    enddo

    !     set into lower case
    !     FC3, ifort: error while passing i80 to a subroutine
    !     do not chane to lower case for header
    it=5
    call lowcase(iatmp,it)

    it=i80
    !     check for tabs, and replace them by blanks
    call tabchk(iatmp,it)

    !     check for blanks at left and remove them, fill with blanks at the end
    call lftpos(iatmp,it)

    do i=1,i80
       ia(i)=iatmp(i)
    enddo


    !     check for an exclamation mark or a hash
    do ike=1,i80
       ile=ike
       if (ia(ile).eq.exm.or.ia(ile).eq.hash) go to 910
    enddo

910 if(ile.eq.1) go to 1
    !     echo the slightly modified input data
    write(iout6,101) ia

    if(ia(ile).eq.exm.or.ia(ile).eq.hash) ile=ile-1

    do i = 1,ile
       if (ia(i).eq.blnk) then
          isw=0
          goto 30
       else
          if (isw.le.0) then
             jump = jump +1
             istrt(jump) = i
             inumb(jump) = 0
             isw=1
             inumb(jump) = inumb(jump) + 1
          else
             inumb(jump) = inumb(jump) + 1
          endif
       endif
30     continue
    enddo

    if (jump.eq.0) goto 1

    do i=6,i80
       header1(i-5)=ia(i)
    enddo

    do i=i80,i80+5
       header1(i-5)=blnk
    enddo

    header80=h80

    return

40  write(iout6,45)
45  format(//1x,'Error: unexpected end of input data. Is STOP label missing?'//)
    stop 'inCardh'
  end subroutine inCardh

  ! ### lftpos ###
  !
  ! Eliminates blanks to the left and left position chararcter string card.
  !
  subroutine lftpos(line,length)
    use params
    implicit none

    character*1, dimension(length) :: line

    integer :: ieff,ipos,length,ntest

    ieff = 0
    do ipos = 1, length
       if(ieff.gt.0) then
          ieff=ieff+1
          line(ieff) = line(ipos)
       end if
       if(line(ipos).ne.' '.and.ieff.eq.0) then
          ieff=1
          line(ieff) = line(ipos)
       end if
    end do

    !     fill end with trailing blanks

    do ipos = ieff+1,length
       line(ipos) = ' '
    end do

    ntest = 0
    if(ntest.ne.0) then
       write(6,*) ' Left adjusted character string '
       write(6,'(1H ,A)') line
    end if

    return
  end subroutine lftpos

  ! ### lowcase ###
  !
  ! Converts letters in a character string line to the lower case.

  subroutine lowcase(line,length)

    use params

    implicit none

    character*1, dimension(26) :: lower,upper
    character*1, dimension(length) :: line

    integer :: i,icha,length

    data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/
    data upper/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

    do icha = 1, length
       do i = 1,26
          if ( line(icha).eq.upper(i) ) line(icha) = lower(i)
       end do
    end do

  end subroutine lowcase


  ! ### tabchk ###
  !
  !     Searches for a tab in the string line and replace it by a space.

  subroutine tabchk(line,length)

    use params

    implicit none

    character*1 :: itab
    character*1, dimension(length) :: line

    integer :: i,length

    itab=char(9)
    do i=1,length
       if(line(i).eq.itab) then
          line(i)=' '
       end if
    end do

  end subroutine tabchk





  ! ### inCardG ###
  !
  !     This routine reads  a single data card and
  !
  !       - changes upper case letters into a lower case ones
  !       - changes tabs (if any) into spaces
  !       - removes leading spaces
  !
  !     Then it scans the line for nonspace fields. The number of fields
  !     is stored in jump, the starting point of a given field in istrt(i)
  !     and the number of characters in that field in inumb(i). If the
  !     routine finds an exclamation mark or a hash then anything on the
  !     same line that follows these characters is treated as a comment.

  !     If iecho is nonzero the (modified) card is written to the standard output.

  subroutine inCardG(iunit,iecho)

    use params
    use card

    implicit none

    character*1 :: blnk,exm,hash
    character*1, dimension(80) :: iatmp

    integer :: i,iecho,ike,ile,isw,it,iunit

    data blnk/' '/,exm/'!'/,hash/'#'/

    do i=1,i80
       ia(i)=' '
    enddo

    !     read a line

    jump = 0
    jrec = 0
    isw = 0
1   read(iunit,100,err=110,end=40) ia
100 format(80a1)
101 format(2x,80a1)
110 continue

    do i=1,i80
       iatmp(i)=ia(i)
    enddo

    !     set into lower case
    !     FC3, ifort: error while passing i80 to a subroutine
    it=i80
    call lowcase(iatmp,it)

    !     check for tabs, and replace them by blanks
    call tabchk(iatmp,it)

    !     check for blanks at left and remove them, fill with blanks at the end
    call lftpos(iatmp,it)

    do i=1,i80
       ia(i)=iatmp(i)
    enddo


    !     check for an exclamation mark or a hash
    do ike=1,i80
       ile=ike
       if (ia(ile).eq.exm.or.ia(ile).eq.hash) go to 910
    enddo

910 if(ile.eq.1) go to 1
    !     echo the slightly modified input data
    if (iecho.ne.0) write(iout6,101) ia

    if(ia(ile).eq.exm.or.ia(ile).eq.hash) ile=ile-1

    do i = 1,ile
       if (ia(i).eq.blnk) then
          isw=0
          goto 30
       else
          if (isw.le.0) then
             jump = jump +1
             istrt(jump) = i
             inumb(jump) = 0
             isw=1
             inumb(jump) = inumb(jump) + 1
          else
             inumb(jump) = inumb(jump) + 1
          endif
       endif
30     continue
    enddo

    if (jump.eq.0) goto 1
    return

40  write(iout6,45)
45  format(//1x,'Error: unexpected end of input data. Is STOP label missing?'//)
    stop 'inCard'
  end subroutine inCardG
end module inCard_m
