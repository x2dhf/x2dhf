! ### inStr ###
!
!    Examines the contents of IA and extracts a character string of up
!    to 8 characters. This string is stored in IBUF.  The remaining
!    non-blank characters (if any) are ignored.

subroutine inStr(guf)

  use params
  use card

  implicit none

  character*1 exm,blnk,hash
  character*80 :: header80,h80

  character*1, dimension(80) :: iatmp
  
  integer :: i,n,nstrt

  character*1 iblnk
  character*8 :: ibufr,guf
  
  character*1, dimension(8) ::  ibuf
  character*1, dimension(60) :: iall

  equivalence (ibuf(1),ibufr)
  
  data iblnk/' '/
  data iall/'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',&
       'u', 'v', 'w', 'x', 'y', 'z', '-', '-', '-', ' ',&
       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J','K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',&
       'U', 'V', 'W', 'X', 'Y', 'Z', '-', '-', '-', ' '/
  
  do i=1,8
     ibuf(i)=iblnk
  enddo
  jrec = jrec + 1
  if(jrec .gt. jump) goto 11
  n = inumb(jrec)
  nstrt = istrt(jrec)
  if (n.gt.8) n=8
  do i = 1,n
     ibuf(i) = ia(nstrt)
     nstrt = nstrt + 1
  enddo
  
11 guf=ibufr
  
end subroutine inStr
