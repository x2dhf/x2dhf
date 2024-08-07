c     jkob 8-04-2021
c     in dint d1mach is called instead if r1mach

c     in I1MACH the following lines were added
c      IMACH(10) = 2
c      IMACH(14) = 53
c      IMACH(15) = -1021
c      IMACH(16) = 1024
c      I1MACH = IMACH(I)
      
#ifdef REAL16
      real*16 function dgamit (a, x)
#else
      real*8 function dgamit (a, x)
#endif
      implicit integer*4 (i-n)

c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.

c
c     evaluate tricomi-s incomplete gamma function defined by
c
c     gamit = x**(-a)/gamma(a) * integral t = 0 to x of exp(-t) *
c     t**(a-1.)

c     and analytic continuation for a .le. 0.0.  gamma(x) is the
c     complete gamma function of x.  gamit is evaluated for arbitrary
c     real values of a and for non-negative values of x (even though
c     gamit is defined for x .lt. 0.0).


c     a slight deterioration of 2 or 3 digits accuracy will occur when
c     gamit is very large or very small in absolute value, because
c     logarithmic variables are used.  also, if the parameter a is very
c     close to a negative integer*4 (but not a negative integer*4), there is
c     a loss of accuracy, which is reported if the result is less than
c     half machine precision.

c
c     ref. -- w. gautschi, an evaluation procedure for incomplete gamma
c     functions, acm trans. math. software.
c
#ifdef REAL16
      real*16 a, x, aeps, ainta, algap1, alneps, alng, alx,
     1  bot, h, sga, sgngam, sqeps, t, d1mach, dgamr, d9gmit, d9lgit,
     2  dlngam, d9lgic, dint
#else
      real*8 a, x, aeps, ainta, algap1, alneps, alng, alx,
     1  bot, h, sga, sgngam, sqeps, t, d1mach, dgamr, d9gmit, d9lgit,
     2  dlngam, d9lgic, dint
#endif
      external d1mach, d9gmit, d9lgic, d9lgit, dgamr, dint,
     1  dlngam

      data alneps, sqeps, bot / 3*0.d0 /

      it1=1
      it3=3
      it4=4
      if (alneps.ne.0.d0) go to 10
      alneps = -log (d1mach(it3))
      sqeps = sqrt (d1mach(it4))
      bot = log (d1mach(it1))
c
cc 10   if (x.lt.0.d0) call seteru (21hdgamit  x is negative, 21, 2, 2)

      it1=1
      it2=2
      it21=21
      it32=32

 10   if (x.lt.0.d0)call seteru ("dgamit  x is negative               ",
     &     it21, it2, it2)

      if (x.ne.0.d0) alx = log (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = sign (1.0d0, a)
      ainta = aint (a + 0.5d0*sga)
      aeps = a - ainta
c
      if (x.gt.0.d0) go to 20
      dgamit = 0.0d0
      if (ainta.gt.0.d0 .or. aeps.ne.0.d0) dgamit = dgamr(a+1.0d0)
      return
c
 20   if (x.gt.1.d0) go to 30
      if (a.ge.(-0.5d0) .or. aeps.ne.0.d0) call dlgams (a+1.0d0, algap1,
     1     sgngam)
      dgamit = d9gmit (a, x, algap1, sgngam, alx)
      return
c
 30   if (a.lt.x) go to 40
      t = d9lgit (a, x, dlngam(a+1.0d0))
      if (t.lt.bot) call erroff
      dgamit = exp (t)
      return
c
 40   alng = d9lgic (a, x, alx)
c
c evaluate dgamit in terms of log (dgamic (a, x))
c
      h = 1.0d0
      if (aeps.eq.0.d0 .and. ainta.le.0.d0) go to 50
c
      call dlgams (a+1.0d0, algap1, sgngam)
      t = log (abs(a)) + alng - algap1
      if (t.gt.alneps) go to 60
c
      if (t.gt.(-alneps)) h = 1.0d0 - sga * sgngam * exp(t)
      if (abs(h).gt.sqeps) go to 50
c
      call erroff
cc      call seteru (32hdgamit  result lt half precision, 32, 1, 1)
      call seteru ("dgamit  result lt half precision    ",it32,it1,it1)
c
 50   t = -a*alx + log(abs(h))
      if (t.lt.bot) call erroff
      dgamit = sign (exp(t), h)
      return
c
 60   t = t - a*alx
      if (t.lt.bot) call erroff
      dgamit = -sga * sgngam * exp(t)
      return
c
      end

#ifdef REAL16
      real*16 function csevl (x, cs, n)
#else
      real*8 function csevl (x, cs, n)
#endif
      implicit integer*4 (i-n)

c     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c     evaluate the n-term chebyshev series cs at x.  adapted from
c     r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
c     and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c     input arguments --

c        x      value at which the series is to be evaluated.
c        cs     array of n terms of a chebyshev series.  in eval-
c               uating cs, only half the first coef is summed.
c        n      number of terms in array cs.
c
      dimension cs(*)
c
      it1=1
      it2=2
      it25=25
      it28=28
      it31=31


cc      if (n.lt.1) call seteru (28hcsevl   number of terms le 0, 28, 2,2)
      if (n.lt.1) call seteru ("csevl   number of terms le 0        ",
     1     it28, it2,it2)
cc      if (n.gt.1000) call seteru (31hcsevl   number of terms gt 1000,
      if (n.gt.1000) call seteru
c     1     ("csevl   number of terms gt 1000     ", 1  31, 3, 2)
     1     ("csevl   number of terms gt 1000     ", it31,it3,it2)

      if (x.lt.(-1.1) .or. x.gt.1.1) call seteru (
cc     1  25hcsevl   x outside (-1,+1), 25, 1, 1)
     1  "csevl   x outside (-1,+1)           ",it25,it1,it1)
c
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
c
      csevl = 0.5 * (b0-b2)
c
      return
      end
      subroutine d9gaml (xmin, xmax)
      implicit integer*4 (i-n)

c     june 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     calculate the minimum and maximum legal bounds for x in gamma(x).
c     xmin and xmax are not the only bounds, but they are the only non-
c     trivial ones to calculate.
c

c     output arguments --

c     xmin   dble prec minimum legal value of x in gamma(x).  any smaller
c            value of x might result in underflow.
c     xmax   dble prec maximum legal value of x in gamma(x).  any larger
c             value of x might cause overflow.
c
#ifdef REAL16
      real*16 xmin, xmax, alnbig, alnsml, xln, xold, d1mach
#else
      real*8 xmin, xmax, alnbig, alnsml, xln, xold, d1mach
#endif

      external d1mach

      it1=1
      it2=2
      it27=27

      alnsml = log(d1mach(it1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)
     1    / (xmin*xln+0.5d0)
        if (abs(xmin-xold).lt.0.005d0) go to 20
 10   continue
cc      call seteru (27hd9gaml  unable to find xmin, 27, 1, 2)
      call seteru ("9gaml  unable to find xmin          ",it27,it1,it2)
c
 20   xmin = -xmin + 0.01d0
c
      alnbig = log (d1mach(it2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)
     1    / (xmax*xln-0.5d0)
        if (abs(xmax-xold).lt.0.005d0) go to 40
 30   continue
cc      call seteru (27hd9gaml  unable to find xmax, 27, 2, 2)
      call seteru ("9gaml  unable to find xmax          ",it27,it2,it2)
c
 40   xmax = xmax - 0.01d0
      xmin = dmax1 (xmin, -xmax+1.d0)
c
      return
      end

#ifdef REAL16
      real*16 function d9gmit (a, x, algap1, sgngam, alx)
#else
      real*8 function d9gmit (a, x, algap1, sgngam, alx)
#endif
      implicit integer*4 (i-n)

c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute tricomi-s incomplete gamma function for small x.
c
#ifdef REAL16
      real*16 a, x, algap1, sgngam, alx, ae, aeps, algs, alg2,
#else
      real*8 a, x, algap1, sgngam, alx, ae, aeps, algs, alg2,
#endif
     1  bot, eps, fk, s, sgng2, t, te, d1mach, dlngam
      data eps, bot / 2*0.d0 /

      it1=1
      it2=2
      it3=3
      it24=24
      it54=54

      if (eps.ne.0.d0) go to 10
      eps = 0.5d0*d1mach(it3)
      bot = log (d1mach(it1))
c
cc 10   if (x.le.0.d0) call seteru (24hd9gmit  x should be gt 0, 24, 1, 2)
 10   if (x.le.0.d0) call seteru
     1     ("d9gmit  x should be gt 0            ",it24,it1,it2)
c
      ma = a + 0.5d0
      if (a.lt.0.d0) ma = a - 0.5d0
      aeps = a - dble(float(ma))
c
      ae = a
      if (a.lt.(-0.5d0)) ae = aeps
c
      t = 1.d0
      te = ae
      s = t
      do 20 k=1,200
        fk = k
        te = -x*te/fk
        t = te/(ae+fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 30
 20   continue
cc      call seteru (54hd9gmit  no convergence in 200 terms of taylor-s se
cc     1ries, 54, 2, 2)
      call seteru ("d9gmit  no conv. in 200 terms of taylors series",
     1    it54, it2, it2)

c
 30   if (a.ge.(-0.5d0)) algs = -algap1 + log(s)
      if (a.ge.(-0.5d0)) go to 60
c
      algs = -dlngam(1.d0+aeps) + log(s)
      s = 1.0d0
      m = -ma - 1
      if (m.eq.0) go to 50
      t = 1.0d0
      do 40 k=1,m
        t = x*t/(aeps-dble(float(m+1-k)))
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 50
 40   continue
c
 50   d9gmit = 0.0d0
      algs = -dble(float(ma))*log(x) + algs
      if (s.eq.0.d0 .or. aeps.eq.0.d0) go to 60
c
      sgng2 = sgngam * sign (1.0d0, s)
      alg2 = -x - algap1 + log(abs(s))
c
      if (alg2.gt.bot) d9gmit = sgng2 * exp(alg2)
      if (algs.gt.bot) d9gmit = d9gmit + exp(algs)
      return
c
 60   d9gmit = exp (algs)
      return
      end

#ifdef REAL16
      real*16 function d9lgic (a, x, alx)
#else
      real*8 function d9lgic (a, x, alx)
#endif

      implicit integer*4 (i-n)

c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute the log complementary incomplete gamma function for large
c     x and for a .le. x.
c
#ifdef REAL16
      real*16 a, x, alx, eps, fk, p, r, s, t, xma, xpa,
#else
      real*8 a, x, alx, eps, fk, p, r, s, t, xma, xpa,
#endif
     1  d1mach
      external d1mach

#ifdef REAL16
      data eps / 0.d0 /
#else
      data eps / 0.0 /
#endif

      it1=1
      it2=2
      it3=3
      it57=57

      if (eps.eq.0.d0) eps = 0.5d0*d1mach(it3)
c
      xpa = x + 1.0d0 - a
      xma = x - 1.d0 - a
c
      r = 0.d0
      p = 1.d0
      s = p
      do 10 k=1,300
        fk = k
        t = fk*(a-fk)*(1.d0+r)
        r = -t/((xma+2.d0*fk)*(xpa+2.d0*fk)+t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 20
 10   continue
c$$$      call seteru (57hd9lgic  no convergence in 300 terms of continued f
c$$$     1raction, 57, 1, 2)
      call seteru ("d9lgic  no conv in 300 terms of contd fraction",
     1 it57, it1, it2)
c
 20   d9lgic = a*alx - x + log(s/xpa)
c
      return
      end

#ifdef REAL16
      real*16 function d9lgit (a, x, algap1)
#else
      real*8 function d9lgit (a, x, algap1)
#endif
      implicit integer*4 (i-n)

c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute the log of tricomi-s incomplete gamma function with
c     perron-s continued fraction for large x and for a .ge. x.
c
#ifdef REAL16
      real*16 a, x, algap1, ax, a1x, eps, fk, hstar, p, r, s,
#else
      real*8 a, x, algap1, ax, a1x, eps, fk, hstar, p, r, s,
#endif
     1  sqeps, t, d1mach
      external d1mach

#ifdef REAL16
      data eps, sqeps / 2*0.d0 /
#else
      data eps, sqeps / 2*0.0 /
#endif

      it1=1
      it2=2
      it3=3
      it4=4
      it7=7
      it35=35
      it39=39

      if (eps.ne.0.d0) go to 10
      eps = 0.5d0*d1mach(it3)
      sqeps = sqrt (d1mach(it4))
c
c$$$ 10   if (x.le.0.d0 .or. a.lt.x) call seteru (
c$$$     1  35hd9lgit  x should be gt 0.0 and le a, 35, 2, 2)

 10   if (x.le.0.d0 .or. a.lt.x) call seteru (
     1  "d9lgit  x should be gt 0.0 and le a ",it35,it2,it2)

c
      ax = a + x
      a1x = ax + 1.0d0
      r = 0.d0
      p = 1.d0
      s = p
      do 20 k=1,200
        fk = k
        t = (a+fk)*x*(1.d0+r)
        r = t/((ax+fk)*(a1x+fk)-t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 30
 20   continue
c$$$      call seteru (57hd9lgit  no convergence in 200 terms of continued f
c$$$     1raction, 57, 3, 2)
      call seteru ("d9lgit  no conv in 200 terms of contd fraction",
     5     it7, it3, it2)
c
 30   hstar = 1.0d0 - x*s/a1x
      if (hstar.lt.sqeps) call seteru (
CC     1  39hd9lgit  result less than half precision, 39, 1, 1)
     1  "d9lgit  result less than half precision",it39,it1,it1)
c
      d9lgit = -x - algap1 - log(hstar)
      return
c
      end


#ifdef REAL16
      real*16 function d9lgmc (x)
#else
      real*8 function d9lgmc (x)
#endif
      implicit integer*4 (i-n)

c     august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute the log gamma correction factor for x .ge. 10. so that
c     log (dgamma(x)) = log(dsqrt(2*pi)) + (x-.5)*log(x) - x +
c     d9lgmc(x)
c
#ifdef REAL16
      real*16 x, algmcs(15), xbig, xmax, dcsevl, d1mach
#else
      real*8 x, algmcs(15), xbig, xmax, dcsevl, d1mach
#endif
      external d1mach, dcsevl, initds

c     series for algm       on the interval  0.          to  1.00000e-02

c                                        with weighted error   1.28e-31
c                                         log weighted error  30.89
c                               significant figures required  29.81
c                                    decimal places required  31.48

      data algmcs(  1) / +.1666389480 4518632472 0572965082 2 d+0      /
      data algmcs(  2) / -.1384948176 0675638407 3298605913 5 d-4      /
      data algmcs(  3) / +.9810825646 9247294261 5717154748 7 d-8      /
      data algmcs(  4) / -.1809129475 5724941942 6330626671 9 d-10     /
      data algmcs(  5) / +.6221098041 8926052271 2601554341 6 d-13     /
      data algmcs(  6) / -.3399615005 4177219443 0333059966 6 d-15     /
      data algmcs(  7) / +.2683181998 4826987489 5753884666 6 d-17     /
      data algmcs(  8) / -.2868042435 3346432841 4462239999 9 d-19     /
      data algmcs(  9) / +.3962837061 0464348036 7930666666 6 d-21     /
      data algmcs( 10) / -.6831888753 9857668701 1199999999 9 d-23     /
      data algmcs( 11) / +.1429227355 9424981475 7333333333 3 d-24     /
      data algmcs( 12) / -.3547598158 1010705471 9999999999 9 d-26     /
      data algmcs( 13) / +.1025680058 0104709120 0000000000 0 d-27     /
      data algmcs( 14) / -.3401102254 3167487999 9999999999 9 d-29     /
      data algmcs( 15) / +.1276642195 6300629333 3333333333 3 d-30     /

#ifdef REAL16
      data nalgm, xbig, xmax / 0, 2*0.d0 /
#else
      data nalgm, xbig, xmax / 0, 2*0.d0 /
#endif

      it0=0
      it1=1
      it2=2
      it3=3
      it15=15
      it23=23
      it34=34
      if (nalgm.ne.0) go to 10
      nalgm = initds (algmcs, it15, sngl(d1mach(it3)) )
      xbig = 1.0d0/sqrt(d1mach(it3))
      xmax =exp(dmin1(log(d1mach(it2)/12.d0), -log(12.d0*d1mach(it1))))

C 10   if (x.lt.10.d0) call seteru (23hd9lgmc  x must be ge 10, 23, 1, 2)
 10   if (x.lt.10.d0) call seteru
     1     ("d9lgmc  x must be ge 10             ",it23, it1, it2)
      if (x.ge.xmax) go to 20
c
      d9lgmc = 1.d0/(12.d0*x)
      if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs,
     1  nalgm) / x
      return
c
 20   d9lgmc = 0.d0
C      call seteru (34hd9lgmc  x so big d9lgmc underflows, 34, 2, 0)
      call seteru ("d9lgmc  x so big d9lgmc underflows  ",it34,it2,it0)
      return
      end

#ifdef REAL16
      real*16 function d9pak (y, n)
#else
      real*8 function d9pak (y, n)
#endif

      implicit integer*4 (i-n)

c     december 1979 edition. w. fullerton, c3, los alamos scientific lab.
c
c     pack a base 2 exponent into floating point number x.  this routine is
c     almost the inverse of d9upak.  it is not exactly the inverse, because
c     abs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
c     were known to be in range we could compute
c     d9pak = x * 2.0d0**n

#ifdef REAL16
      real*16 y, aln2b, aln210, d1mach
#else
      real*8 y, aln2b, aln210, d1mach
#endif

      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094 8873623478 7031942948 9 d0 /

      it0=0
      it1=1
      it2=2
      it5=5
      it10=10
      it15=15
      it16=16
      it31=31
      it32=32
      if (nmin.ne.0) go to 10
      aln2b = 1.0d0
      print *,'1'
      if (i1mach(it10).ne.2) aln2b = d1mach(it5)*aln210
      nmin = aln2b*dble(float(i1mach(it15)))
      nmax = aln2b*dble(float(i1mach(it16)))
c
 10   call d9upak (y, d9pak, ny)
c
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
      if (nsum.gt.nmax) call seteru (
c     1  31hd9pak   packed number overflows, 31, 1, 2)
     1  "d9pak   packed number overflows     ",it31,it1,it2)
c
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
c
 20   d9pak = 0.5d0*d9pak
      nsum = nsum + 1
      if (nsum.ne.0) go to 20
      return
c
 30   d9pak = 2.0d0*d9pak
      nsum = nsum - 1
      if (nsum.ne.0) go to 30
      return
c
c 40   call seteru (32hd9pak   packed number underflows, 32, 1, 0)
 40   call seteru ("d9pak   packed number underflows    ",it32,it1,it0)
      d9pak = 0.0d0
      return
c
      end
      subroutine d9upak (x, y, n)
      implicit integer*4 (i-n)
c     august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c     unpack floating point number x so that x = y * 2.0**n, where
c     0.5 .le. abs(y) .lt. 1.0 .
c
#ifdef REAL16
      real*16 x, y, absx
#else
      real*8 x, y, absx
#endif

      absx = abs(x)
      n = 0
      y = 0.0d0
      if (x.eq.0.0d0) return

 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10

 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20

 30   y = sign (absx, x)
      return
      end

#ifdef REAL16
      real*16 function dcsevl (x, a, n)
#else
      real*8 function dcsevl (x, a, n)
#endif
      implicit integer*4 (i-n)

c     evaluate the n-term chebyshev series a at x.  adapted from
c     r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c     input arguments --
c       x      dble prec value at which the series is to be evaluated.
c       a      dble prec array of n terms of a chebyshev series.  in eval-
c              uating a, only half the first coef is summed.
c       n      number of terms in array a.
c
#ifdef REAL16
      real*16 a(n), x, twox, b0, b1, b2
#else
      real*8 a(n), x, twox, b0, b1, b2
#endif

c      if (n.lt.1) call seteru (28hdcsevl  number of terms le 0, 28, 2,2)

      it1=1
      it2=2
      it3=3
      it25=25
      it28=28
      it31=31

      if (n.lt.1) call seteru ("dcsevl  number of terms le 0        ",
     1     it28, it2,it2)
c      if (n.gt.1000) call seteru (31hdcsevl  number of terms gt 1000,
      if (n.gt.1000) call seteru
     1   ("dcsevl  number of terms gt 1000       ",it31,it3,it2)
c     1   ("dcsevl  number of terms gt 1000       ", 1 31,it3,it2)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru (
c     1  25hdcsevl  x outside (-1,+1), 25, 1, 1)
     1  "dcsevl  x outside (-1,+1)           ", it25, it1, it1)
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end

#ifdef REAL16
      real*16 function dgamma (x)
#else
      real*8 function dgamma (x)
#endif

      implicit integer*4 (i-n)
c     jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
c     jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml

#ifdef REAL16
      real*16 x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1    xmin, y, d9lgmc, dcsevl, d1mach, dint, xsml
#else
      real*8 x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1    xmin, y, d9lgmc, dcsevl, d1mach, dint, xsml      
#endif

      external d1mach, d9lgmc, dcsevl, dint,  initds

c     series for gam        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.00
c                                    decimal places required  32.05
c
      data gam cs(  1) / +.8571195590 9893314219 2006239994 2 d-2      /
      data gam cs(  2) / +.4415381324 8410067571 9131577165 2 d-2      /
      data gam cs(  3) / +.5685043681 5993633786 3266458878 9 d-1      /
      data gam cs(  4) / -.4219835396 4185605010 1250018662 4 d-2      /
      data gam cs(  5) / +.1326808181 2124602205 8400679635 2 d-2      /
      data gam cs(  6) / -.1893024529 7988804325 2394702388 6 d-3      /
      data gam cs(  7) / +.3606925327 4412452565 7808221722 5 d-4      /
      data gam cs(  8) / -.6056761904 4608642184 8554829036 5 d-5      /
      data gam cs(  9) / +.1055829546 3022833447 3182350909 3 d-5      /
      data gam cs( 10) / -.1811967365 5423840482 9185589116 6 d-6      /
      data gam cs( 11) / +.3117724964 7153222777 9025459316 9 d-7      /
      data gam cs( 12) / -.5354219639 0196871408 7408102434 7 d-8      /
      data gam cs( 13) / +.9193275519 8595889468 8778682594 0 d-9      /
      data gam cs( 14) / -.1577941280 2883397617 6742327395 3 d-9      /
      data gam cs( 15) / +.2707980622 9349545432 6654043308 9 d-10     /
      data gam cs( 16) / -.4646818653 8257301440 8166105893 3 d-11     /
      data gam cs( 17) / +.7973350192 0074196564 6076717535 9 d-12     /
      data gam cs( 18) / -.1368078209 8309160257 9949917230 9 d-12     /
      data gam cs( 19) / +.2347319486 5638006572 3347177168 8 d-13     /
      data gam cs( 20) / -.4027432614 9490669327 6657053469 9 d-14     /
      data gam cs( 21) / +.6910051747 3721009121 3833697525 7 d-15     /
      data gam cs( 22) / -.1185584500 2219929070 5238712619 2 d-15     /
      data gam cs( 23) / +.2034148542 4963739552 0102605193 2 d-16     /
      data gam cs( 24) / -.3490054341 7174058492 7401294910 8 d-17     /
      data gam cs( 25) / +.5987993856 4853055671 3505106602 6 d-18     /
      data gam cs( 26) / -.1027378057 8722280744 9006977843 1 d-18     /
      data gam cs( 27) / +.1762702816 0605298249 4275966074 8 d-19     /
      data gam cs( 28) / -.3024320653 7353062609 5877211204 2 d-20     /
      data gam cs( 29) / +.5188914660 2183978397 1783355050 6 d-21     /
      data gam cs( 30) / -.8902770842 4565766924 4925160106 6 d-22     /
      data gam cs( 31) / +.1527474068 4933426022 7459689130 6 d-22     /
      data gam cs( 32) / -.2620731256 1873629002 5732833279 9 d-23     /
      data gam cs( 33) / +.4496464047 8305386703 3104657066 6 d-24     /
      data gam cs( 34) / -.7714712731 3368779117 0390152533 3 d-25     /
      data gam cs( 35) / +.1323635453 1260440364 8657271466 6 d-25     /
      data gam cs( 36) / -.2270999412 9429288167 0231381333 3 d-26     /
      data gam cs( 37) / +.3896418998 0039914493 2081663999 9 d-27     /
      data gam cs( 38) / -.6685198115 1259533277 9212799999 9 d-28     /
      data gam cs( 39) / +.1146998663 1400243843 4761386666 6 d-28     /
      data gam cs( 40) / -.1967938586 3451346772 9510399999 9 d-29     /
      data gam cs( 41) / +.3376448816 5853380903 3489066666 6 d-30     /
      data gam cs( 42) / -.5793070335 7821357846 2549333333 3 d-31     /

      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /

      it0=0
      it1=1
      it2=2
      it3=3
      it4=4
      it5=5
      it14=14
      it31=31
      it32=32
      it35=0
      it42=42
      it54=54
      it61=61
      it68=68

      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, it42, 0.1*sngl(d1mach(it3)) )
c
      it4=4
      call d9gaml (xmin, xmax)
      xsml = exp (dmax1 (log(d1mach(it1)), -log(d1mach(it2)))+0.01d0)
      dxrel = sqrt (d1mach(it4))
c
 10   y = abs(x)
      if (y.gt.10.d0) go to 50
c
c     compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
c     gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c     compute gamma(x) for x .lt. 1.0
c
      n = -n
c      if (x.eq.0.d0) call seteru (14hdgamma  x is 0, 14, 4, 2)
      if (x.eq.0.d0) call seteru
c     1     ("dgamma  x is 0                      ",  1     14, it4, it2)
     1     ("dgamma  x is 0                      ", it14, it4, it2)

      if (x.lt.0.0d0 .and. x+dble(float(n-2)).eq.0.d0) call seteru (
c     1  31hdgamma  x is a negative integer*4, 31, 4, 2)
     1  "dgamma  x is a negative integer*4     ", it31, it4, it2)
      if (x.lt.(-0.5d0) .and. abs((x-dint(x-0.5d0))/x).lt.dxrel) call
cc     1  seteru (68hdgamma  answer lt half precision because x too near n
cc     2egative integer*4, 68, 1, 1)
     1  seteru ( "dgamma  answer lt half precision because x too near n
     2egative integer*4", it68, it1, it1)

      if (y.lt.xsml) call seteru (
c     1  54hdgamma  x is so close to 0.0 that the result overflows,
     1  "dgamma  x is so close to 0.0 that the result overflows",
     2  it54, it5, it2)
c
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
c
c     gamma(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
c
c     gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).
c
c 50   if (x.gt.xmax) call seteru (32hdgamma  x so big gamma overflows,
 50   if (x.gt.xmax) call seteru
c     1     ("dgamma  x so big gamma overflows    ", 1  32, 3, 2)
     1     ("dgamma  x so big gamma overflows    ",it32, it3, it2)
c
      dgamma = 0.d0
c      if (x.lt.xmin) call seteru (35hdgamma  x so small gamma underflows
      if (x.lt.xmin) call seteru ("dgamma  x so small gamma underflows "
     1  , it35, it2, it0)
      if (x.lt.xmin) return
c
      dgamma = exp ((y-0.5d0)*log(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
c
      if (abs((x-dint(x-0.5d0))/x).lt.dxrel) call seteru (
c     1  61hdgamma  answer lt half precision, x too near negative integer*4
     1  "dgamma  answer lt half precision, x too near negative integer"
     2  , it61, it1, it1)
c
      sinpiy = sin (pi*y)
      if (sinpiy.eq.0.d0) call seteru (
c     1  31hdgamma  x is a negative integer*4, 31, 4, 2)
     1  "dgamma  x is a negative integer     ", it31, it4, it2)
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end

#ifdef REAL16
      real*16 function dgamr (x)
#else
      real*8 function dgamr (x)
#endif

      implicit integer*4 (i-n)
c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c     this routine, not dgamma(x), should be the fundamental one.

#ifdef REAL16
      real*16 x, alngx, sgngx, dgamma, dint
#else
      real*8 x, alngx, sgngx, dgamma, dint
#endif
      external  dgamma, dint

      dgamr = 0.0d0
      if (x.le.0.0d0 .and. dint(x).eq.x) return

      it1=1
      call entsrc (irold, it1)
      if (abs(x).gt.10.0d0) go to 10
      dgamr = 1.0d0/dgamma(x)
      call erroff
      call entsrc (ir, irold)
      return

 10   call dlgams (x, alngx, sgngx)
      call erroff
      call entsrc (ir, irold)
      dgamr = sgngx * exp(-alngx)
      return

      end

#ifdef REAL16
      real*16 function dint (x)
#else
      real*8 function dint (x)
#endif
      implicit integer*4 (i-n)

c     december 1983 edition. w. fullerton, c3, los alamos scientific lab.
c
c     dint is the real*16 equivalent of aint.  this portable
c     version is quite efficient when the argument is reasonably small (a
c     common case), and so no faster machine-dependent version is needed.

#ifdef REAL16
      real*16 x, xscl, scale, xbig, xmax, part, d1mach, r1mach
#else
      real*8 x, xscl, scale, xbig, xmax, part, d1mach, r1mach
#endif
      external d1mach,  i1mach, r1mach
      data npart, scale, xbig, xmax / 0, 3*0.0d0 /

      if (npart.ne.0) go to 10

      it4=4
      it9=9
      it10=10

      ibase = i1mach(it10)
      xmax = 1.0d0/d1mach(it4)
      xbig = dmin1 (dble (i1mach(it9)), 1.d0/d1mach(it4))
      scale = ibase**int(log(xbig)/log(dble(float(ibase)))-0.5d0)
      npart = log(xmax)/log(scale) + 1.0d0

 10   if (x.lt.(-xbig) .or. x.gt.xbig) go to 20

      dint = int(sngl(x))
      return

 20   xscl = abs(x)
      if (xscl.gt.xmax) go to 50

      do 30 i=1,npart
        xscl = xscl/scale
 30   continue

      dint = 0.0d0
      do 40 i=1,npart
        xscl = xscl*scale
        ipart = xscl
        part = ipart
        xscl = xscl - part
        dint = dint*scale + part
 40   continue

      if (x.lt.0.0d0) dint = -dint
      return

      it1=1
 1    it68=68
 50   call seteru ("dint abs(x) may be too big to be represented a
     1s an exact integer*4", it68, it1, it1)
      dint = x
      return
      end

      subroutine dlgams (x, dlgam, sgngam)
      implicit integer*4 (i-n)
c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.

c     evaluate log abs (gamma(x)) and return the sign of gamma(x) in sgngam.
c     sgngam is either +1.0 or -1.0.

#ifdef REAL16
      real*16 x, dlgam, sgngam, dint, dlngam
#else
      real*8 x, dlgam, sgngam, dint, dlngam
#endif
      external dint, dlngam

      dlgam = dlngam(x)
      sgngam = 1.0d0
      if (x.gt.0.d0) return

      int = mod (-aint(x), 2.0d0) + 0.1d0
      if (int.eq.0) sgngam = -1.0d0

      return
      end

#ifdef REAL16
      real*16 function dlngam (x)
#else
      real*8 function dlngam (x)
#endif
      implicit integer*4 (i-n)

c     august 1980 edition.   w. fullerton, c3, los alamos scientific lab.
#ifdef REAL16
      real*16 x, dxrel, pi, sinpiy, sqpi2l, sq2pil,
#else
      real*8 x, dxrel, pi, sinpiy, sqpi2l, sq2pil,
#endif
     1  y, xmax, dint, dgamma, d9lgmc, d1mach
      external d1mach, d9lgmc, dgamma, dint

      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
c     sq2pil = alog (sqrt(2*pi)),  sqpi2l = alog(sqrt(pi/2))
      data sqpi2l / +.2257913526 4472743236 3097614947 441 e+0    /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /

      data xmax, dxrel / 2*0.d0 /

      it1=1
      it2=2
      it3=3
      it4=4
      it31=31
      it39=39
      it68=68
      if (xmax.ne.0.d0) go to 10
      xmax = d1mach(it2)/log(d1mach(it2))
      dxrel = sqrt (d1mach(it4))
c
 10   y = abs (x)
      if (y.gt.10.d0) go to 20
c
c     log (abs (dgamma(x)) ) for abs(x) .le. 10.0
c
      dlngam = log (abs (dgamma(x)) )
      return
c
c     log ( abs (dgamma(x)) ) for abs(x) .gt. 10.0
c
 20   if (y.gt.xmax) call seteru (
c     1  39hdlngam  abs(x) so big dlngam overflows, 39, 2, 2)
     1 "dlngam  abs(x) so big dlngam overflows",it39,it2,it2)
c
      if (x.gt.0.d0) dlngam = sq2pil + (x-0.5d0)*log(x) - x + d9lgmc(y)
      if (x.gt.0.d0) return
c
      sinpiy = abs (sin(pi*y))
      if (sinpiy.eq.0.d0) call seteru (
c     1  31hdlngam  x is a negative integer*4, 31, 3, 2)
     1  "dlngam  x is a negative integer*4     ",it31,it3,it2)
c
      dlngam = sqpi2l + (x-0.5d0)*log(y) - x - log(sinpiy) - d9lgmc(y)
c
      if (abs((x-dint(x-0.5d0))*dlngam/x).lt.dxrel) call seteru (
c     1  68hdlngam  answer lt half precision because x too near negative
     1 "dlngam  answer lt half precision because x too near negative
     2integer*4", it68, it1, it1)
      return
      end

      subroutine e9rint(messg,nw,nerr,save)
      implicit integer*4 (i-n)
c
c     this routine stores the current error message or prints the old one,
c     if any, depending on whether or not save = .true. .
c
      character*1 messg(nw)
      logical save
      external i1mach, i8save
c
c     messgp stores at least the first 72 characters of the previous
c     message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer*4 word).
c
c     integer*4 messgp(36),fmt(14),ccplus
      character*1 messgp(36),fmt(14),ccplus
c
c     start with no previous message.
c
      data messgp(1)/'1'/, nwp/0/, nerrp/0/
c
c     set up the format for printing the error message.  the format is
c     simply (a1,14x,72axx) where xx=i1mach(6) is the number of
c     characters stored per integer*4 word.
c
      data ccplus  /'+'/
c
      data fmt( 1) /'('/
      data fmt( 2) /'a'/
      data fmt( 3) /'1'/
      data fmt( 4) /','/
      data fmt( 5) /'1'/
      data fmt( 6) /'4'/
      data fmt( 7) /'x'/
      data fmt( 8) /','/
      data fmt( 9) /'7'/
      data fmt(10) /'2'/
      data fmt(11) /'a'/
      data fmt(12) /'x'/
      data fmt(13) /'x'/
      data fmt(14) /')'/
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do i=1,nw
           messgp(i)=messg(i)
        enddo

        go to 30

        it0=0
        it1=1
 20   if (i8save(it1,it0,.false.).eq.0) go to 30

c  print the message.

      it4=4
      iwunit=i1mach(it4)
      write(iwunit,9000) nerrp
 9000 format(7h error ,i4,4h in )
     
      it2=2
      it6=6
      call s88fmt(it2,i1mach(it6),fmt(12))
c     write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
      print *, ccplus,(messgp(i),i=1,nwp)

 30   return

      end

      subroutine entsrc(irold,irnew)
      implicit integer*4 (i-n)

c     this routine returns irold = lrecov and sets lrecov = irnew.
c
c     if there is an active error state, the message is printed and
c     execution stops.
c
c     irnew = 0 leaves lrecov unchanged, while
c     irnew = 1 gives recovery and
c     irnew = 2 turns recovery off.
c
c     error states -
c
c     1 - illegal value of irnew.
c     2 - called while in an error state.
c
      external i8save
c
      it0=0
      it1=1
      it2=2
      it31=31
      it39=39

      if (irnew.lt.it0 .or. irnew.gt.it2)


c     1   call seterr(31hentsrc - illegal value of irnew,31,1,2)
     1   call seterr("entsrc - illegal value of irnew     ",
     2                it31,it1,it2)

      irold=i8save(it2,irnew,irnew.ne.0)

c  if have an error state, stop execution.

      if (i8save(it1,it0,.false.) .ne. 0) call seterr
c     1   (39hentsrc - called while in an error state,39,2,2)
     1   ("entsrc - called while in an error state",
     2    it39,it2,it2)

      return
      end

      subroutine eprint
      implicit integer*4 (i-n)

c     this subroutine prints the last error message, if any.

c      parameter (messglength=100)
      parameter (messglength=36)
      character*1 messg(messglength)
      it1=1
      call e9rint(messg,it1,it1,.false.)
      return
      end

      subroutine erroff
      implicit integer*4 (i-n)

c     turns off the error state off by setting lerror=0.

      external i8save

      it0=0
      it1=1
      i=i8save(it1,it0,.true.)
      return
      end

      integer*4 function i8save(isw,ivalue,set)
      implicit integer*4 (i-n)

c     if (isw = 1) i8save returns the current error number and
c     sets it to ivalue if set = .true. .
c
c     if (isw = 2) i8save returns the current recovery switch and
c     sets it to ivalue if set = .true. .
c
      logical set
c
      integer*4 iparam(2)
c     iparam(1) is the error number and iparam(2) is the recovery switch.

c     start execution error free and with recovery turned off.

      data iparam(1) /0/,  iparam(2) /2/

      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
      return
      end

      integer*4 function initds (dos, nos, eta)
      implicit integer*4 (i-n)

c     june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c     initialize the real*16 orthogonal series dos so that
c     initds is the number of terms needed to insure the error is no
c     larger than eta.  ordinarily eta will be chosen to be one-tenth
c     machine precision.
c
c     input arguments --
c     dos    dble prec array of nos coefficients in an orthogonal series.
c     nos    number of coefficients in dos.
c     eta    requested accuracy of series.
c
#ifdef REAL16
      real*16 dos(nos)
#else
      real*8 dos(nos)
#endif

      it1=1
      it2=2
      it28=28
      it35=35

      if (nos.lt.1) call seteru (
c     1  35hinitds  number of coefficients lt 1, 35, 2, 2)
     1  "initds  number of coefficients lt 1 ",it35,it2,it2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue

c 20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28,
 20   if (i.eq.nos) call seteru ("initds  eta may be too small        ",
     1     it28, it1, it2)
      initds = i

      return
      end

      integer*4 function inits (os, nos, eta)
      implicit integer*4 (i-n)

c     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c     initialize the orthogonal series so that inits is the number of terms
c     needed to insure the error is no larger than eta.  ordinarily, eta
c     will be chosen to be one-tenth machine precision.
c
c     input arguments --
c     os     array of nos coefficients in an orthogonal series.
c     nos    number of coefficients in os.
c     eta    requested accuracy of series.
c
      dimension os(nos)
c
      it1=1
      it2=2
      it28=28
      it35=35

      if (nos.lt.1) call seteru (
c     1  35hinits   number of coefficients lt 1, 35, 2, 2)
     1  "inits   number of coefficients lt 1 ",it35,it2,it2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue

c 20   if (i.eq.nos) call seteru (28hinits   eta may be too small, 28,
 20   if (i.eq.nos) call seteru ("inits   eta may be too small        ",
     1     it28,  it1 , it2)
c     1     28,  1  1, 2)
      inits = i
c
      return
      end

      subroutine r9upak (x, y, n)
      implicit integer*4 (i-n)
c     august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c     unpack floating point number x so that x = y * 2.0**n, where
c     0.5 .le. abs(y) .lt. 1.0 .
c
      absx = abs(x)
      n = 0
      y = 0.0
      if (x.eq.0.0) return
c
 10   if (absx.ge.0.5) go to 20
      n = n - 1
      absx = absx*2.0
      go to 10
c
 20   if (absx.lt.1.0) go to 30
      n = n + 1
      absx = absx*0.5
      go to 20
c
 30   y = sign (absx, x)
      return
c
      end

      subroutine s88fmt( n, w, ifmt )
      implicit integer*4 (i-n)
c
c     s88fmt  replaces ifmt(1), ... , ifmt(n) with
c     the characters corresponding to the n least significant
c     digits of w.
c
      integer*4 n,nt,w,wt
      character*1 ifmt(n),digits(10)
c
      data digits( 1) /'0'/
      data digits( 2) /'1'/
      data digits( 3) /'2'/
      data digits( 4) /'3'/
      data digits( 5) /'4'/
      data digits( 6) /'5'/
      data digits( 7) /'6'/
      data digits( 8) /'7'/
      data digits( 9) /'8'/
      data digits(10) /'9'/
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
      end

      subroutine seterr (messg, nmessg, nerr, iopt)
      implicit integer*4 (i-n)
c
c     this version modified by w. fullerton to dump if iopt = 1 and
c     not recovering.
c     seterr sets lerror = nerr, optionally prints the message and dumps
c     according to the following rules...
c
c     if iopt = 1 and recovering      - just remember the error.
c     if iopt = 1 and not recovering  - print, dump and stop.
c     if iopt = 2                     - print, dump and stop.
c
c     input
c
c     messg  - the error message.
c     nmessg - the length of the message, in characters.
c     nerr   - the error number. must have nerr non-zero.
c     iopt   - the option. must have iopt=1 or 2.
c
c     error states -
c
c     1 - message length not positive.
c     2 - cannot have nerr=0.
c     3 - an unrecovered error followed by another error.
c     4 - bad value for iopt.
c
c     only the first 72 characters of the message are printed.
c
c     the error handler calls a subroutine named fdump to produce a
c     symbolic dump. to complete the package, a dummy version of fdump
c     is supplied, but it should be replaced by a locally written version
c     which at least gives a trace-back.
c
c      parameter (messglength=100)
      parameter (messglength=36)
      character*1 messg(messglength)
      external i1mach, i8save
c
c  the unit for error messages.
c
      i4=4
      iwunit=i1mach(i4)

      if (nmessg.ge.1) go to 10

c  a message of non-positive length is fatal.
      
      write(iwunit,9000)
 9000 format(52h1error    1 in seterr - message length not positive.)
      go to 60
c       
c  nw is the number of words the message occupies.
c
      i6=6
      i72=72
      
 10   continue
      i6=i1mach(i6)
      nw=(min0(nmessg,i72)-1)/i6+1
      
      it1=1
      if (nerr.ne.0) go to 20
      
c  cannot turn the error state off using seterr.

      write(iwunit,9001)
 9001 format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
      call e9rint(messg,nw,nerr,.true.)
      it1=1
      itemp=i8save(it1,it1,.true.)
      go to 50
     
c  set lerror and test for a previous unrecovered error.
      
 20   continue

      if (i8save(it1,nerr,.true.).eq.0) go to 30

      write(iwunit,9002)
 9002 format(23h1error    3 in seterr -,
     1     48h an unrecovered error followed by another error.//
     2     48h the previous and current error messages follow.///)
      call eprint
      call e9rint(messg,nw,nerr,.true.)
      go to 50
      
c  save this message in case it is not recovered from properly.

 30   call e9rint(messg,nw,nerr,.true.)

      if (iopt.eq.1 .or. iopt.eq.2) go to 40

c  must have iopt = 1 or 2.

      write(iwunit,9003)
 9003 format(42h1error    4 in seterr - bad value for iopt//
     1     34h the current error message follows///)
      go to 50
      
c     test for recovery.
      
 40   if (iopt.eq.2) go to 50
      
      it0=0
      it2=2
      if (i8save(it2,it0,.false.).eq.1) return
      
c     call eprint
c     stop

 50   call eprint

 60   call fdump
      stop 'seterr'
      end

      subroutine seteru (messg, nmessg, nerr, iopt)
      implicit integer*4 (i-n)
c      parameter (messglength=100)
      parameter (messglength=36)
      common /cseter/ iunflo
      character*1 messg(messglength)
      data iunflo / 0 /

      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
      if (iunflo.le.0) return
      it1=1
      call seterr (messg, nmessg, nerr, it1)
      return
      end

      subroutine fdump
      return
      end


