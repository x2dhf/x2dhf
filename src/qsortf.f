C From Leonard J. Moss of SLAC:

C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.  

c      SUBROUTINE SORTRX(N,DATA,INDEX)
      SUBROUTINE QSORTF(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      integer*4 n,index(n)
      real*8    data(n)
 
      integer*4 lstk(31),rstk(31),istk
      integer*4 l,r,i,j,p,indexp,indext
      real*8    datap
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      integer*4  m
      parameter (m=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      do 50 i=1,n
         index(i)=i
   50    continue
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      if (n.le.m) goto 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      istk=0
      l=1
      r=n
 
  200 continue
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      i=l
      j=r
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      p=(l+r)/2
      indexp=index(p)
      datap=data(indexp)
 
      if (data(index(l)) .gt. datap) then
         index(p)=index(l)
         index(l)=indexp
         indexp=index(p)
         datap=data(indexp)
      endif
 
      if (datap .gt. data(index(r))) then
         if (data(index(l)) .gt. data(index(r))) then
            index(p)=index(l)
            index(l)=index(r)
         else
            index(p)=index(r)
         endif
         index(r)=indexp
         indexp=index(p)
         datap=data(indexp)
      endif
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 continue
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         i=i+1
         if (data(index(i)).lt.datap) goto 300
 
  400 continue
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         j=j-1
         if (data(index(j)).gt.datap) goto 400
 
C Q5: Have the two scans collided?
 
      if (i.lt.j) then
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         indext=index(i)
         index(i)=index(j)
         index(j)=indext
         goto 300
      else
 
c Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         if (r-j .ge. i-l .and. i-l .gt. m) then
            istk=istk+1
            lstk(istk)=j+1
            rstk(istk)=r
            r=i-1
         else if (i-l .gt. r-j .and. r-j .gt. m) then
            istk=istk+1
            lstk(istk)=l
            rstk(istk)=i-1
            l=j+1
         else if (r-j .gt. m) then
            l=j+1
         else if (i-l .gt. m) then
            r=i-1
         else
C Q8: Pop the stack, or terminate QuickSort if empty
            if (istk.lt.1) goto 900
            l=lstk(istk)
            r=rstk(istk)
            istk=istk-1
         endif
         goto 200
      endif
 
  900 continue
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      do 950 i=2,n
         if (data(index(i-1)) .gt. data(index(i))) then
            indexp=index(i)
            datap=data(indexp)
            p=i-1
  920       continue
               index(p+1) = index(p)
               p=p-1
               if (p.gt.0) then
                  if (data(index(p)).gt.datap) goto 920
               endif
            index(p+1) = indexp
         endif
  950    continue
 
c===================================================================
C
C     All done
 
      end
