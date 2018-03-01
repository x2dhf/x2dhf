! ### inFloat ###
!
!     Extracts a floating point number from an input card.

subroutine inFloat(buf)
  use params
  use card

  implicit none

  integer :: i,ist,j,n,nstrt
  integer :: iexp,iexpdig,limit

  real (PREC) :: buf,exponent,fact,fact2
  character*1 :: itemp
  character*1, dimension(17) :: ichar

  parameter (limit=17)
  data ichar/'0','1','2','3','4','5','6','7','8','9','+','&','^','-','.','e','d'/

  jrec = jrec + 1
  buf=0.0_PREC
  if(jrec .gt. jump) return
  buf = 0.0_PREC
  fact2 = 0.0_PREC
  iexp=0
  n = inumb(jrec)
  fact = 1.0_PREC
  ist=istrt(jrec)
  nstrt = ist + n - 1
  iexpdig=0
  do i=1,n
     itemp = ia(nstrt)
     do j = 1,limit
        if(ichar(j).eq.itemp) goto 5
     enddo

     write(iout6,*) 'Error detected in inFloat'
     stop 'inFloat'

5    continue
     if (j.eq.16.or.j.eq.17) goto 12
     if (j.lt.11) goto 7
     if (j.le.14) goto 6
     if (iexp.eq.1) then
        fact2 = dble(i-1-iexpdig)
     else
        fact2 = dble(i-1)
     endif
     go to 9

12   continue
     if (iexp.eq.0) then
        buf=(0.10_PREC**fact2)*buf
        exponent=buf
        iexp=1
        buf=0
        fact=1.0
        fact2=0.0
        iexpdig=i
     endif
     goto 9

     buf = buf + dble(j-1) * fact
     fact=fact*10.0_PREC
     goto 9

6    continue
     !         if(nstrt.ne.ist.and.iexp.eq.0) go to 4
     if(j.eq.14) buf=-buf
     goto 9

7    buf = buf + dble(j-1) * fact
     fact=fact*10.0_PREC
9    continue
     nstrt = nstrt - 1
  enddo

  if (iexp.eq.1) then
     buf=(0.10_PREC**fact2)*buf
     buf=buf*10.0_PREC**exponent
     return
  endif
  buf=(0.10_PREC**fact2)*buf

end subroutine inFloat
