subroutine checkLabel

  use input_m

  implicit none
  
  do i=1, nlabels
    if (clabel.eq.labellc(i)) return
  end do

  write(iout6,100) clabel,(labellc(i),i=1,nlabels)
100 format(/,'Error: label "',a8,'" is not supported. ',/,/,&
         'Try one of the following:',/,10(8x,a8,2x,a8,2x,a8,2x,a8,2x,a8/)//)
  stop 'checkLabel'
  
end subroutine checkLabel
