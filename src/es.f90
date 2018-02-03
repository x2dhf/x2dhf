! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1992 F.A. Parpia and I.P. Grant                         *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************

! ****************************************************************************
! *                                                                          *
! * This routine is part of the GRASP2 package (ver. 1.00, 17 Jan 1992)      *
! *                                                                          *
subroutine es (f,s2f,s3f)
! *                                                                          *
! *   ----------------   SECTION 01   SUBPROGRAM 14   ----------------       *
! *   EVALUATE THE SUM OF THE SERIES                                         *
! *                                                                          *
! *                       INFINITY      N              K                     *
! *              S  (F) =   SUM     (-1)  EXP (N*F) / N                      *
! *               K        N = 0                                             *
! *                                                                          *
! *   FOR K = 2, 3 TO MACHINE PRECISION. THIS IS A UTILITY SUBROUTINE,       *
! *   CALLED BY SUBROUTINES NUCPOT AND NCHARG.                               *
! *                                                                          *
! *   NO SUBROUTINES CALLED                                                  *
! *                                                                          *
! *   WRITTEN BY FARID A PARPIA, AT OXFORD  LAST REVISION: 25 JAN 1988       *
! *                                                                          *
! ****************************************************************************

  use params

  implicit none
  integer :: n
  real (PREC) :: en,enf,f,fase,obn,s2f,s3f,s2last,term2,term3

  n = 0
  s2f = zero
  s3f = zero
  fase = one
1 n = n+1
  en = dble(n)
  obn = one/en
  fase = -fase
  enf = exp (en*f)
  term2 = fase*enf*obn*obn
  term3 = term2*obn
  s2last = s2f
  s2f = s2f+term2
  s3f = s3f+term3
  if (abs (s2f) .ne. abs (s2last)) goto 1

end subroutine es
