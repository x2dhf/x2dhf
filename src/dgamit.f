      real*8 function dgamit (a, x) 
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
c     close to a negative integer (but not a negative integer), there is
c     a loss of accuracy, which is reported if the result is less than
c     half machine precision.

c
c     ref. -- w. gautschi, an evaluation procedure for incomplete gamma
c     functions, acm trans. math. software.
c

      real*8 a, x, aeps, ainta, algap1, alneps, alng, alx,
     1  bot, h, sga, sgngam, sqeps, t, d1mach, dgamr, d9gmit, d9lgit,
     2  dlngam, d9lgic, dint
      external d1mach, d9gmit, d9lgic, d9lgit, dgamr, dint,
     1  dlngam
c
      data alneps, sqeps, bot / 3*0.d0 /
c
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

c
      if (x.ne.0.d0) alx = log (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = dsign (1.0d0, a)
      ainta = dint (a + 0.5d0*sga)
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
      dgamit = dsign (exp(t), h)
      return
c
 60   t = t - a*alx
      if (t.lt.bot) call erroff
      dgamit = -sga * sgngam * exp(t)
      return
c
      end

      real*8 function csevl (x, cs, n)
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
      real*8 xmin, xmax, alnbig, alnsml, xln, xold, d1mach
      external d1mach
c
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
      real*8 function d9gmit (a, x, algap1, sgngam, alx)
      implicit integer*4 (i-n)


c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute tricomi-s incomplete gamma function for small x.
c
c
      real*8 a, x, algap1, sgngam, alx, ae, aeps, algs, alg2,
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
      sgng2 = sgngam * dsign (1.0d0, s)
      alg2 = -x - algap1 + log(abs(s))
c
      if (alg2.gt.bot) d9gmit = sgng2 * exp(alg2)
      if (algs.gt.bot) d9gmit = d9gmit + exp(algs)
      return
c
 60   d9gmit = exp (algs)
      return
c
      end
      real*8 function d9lgic (a, x, alx)
      implicit integer*4 (i-n)

c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute the log complementary incomplete gamma function for large
c     x and for a .le. x.
c
      real*8 a, x, alx, eps, fk, p, r, s, t, xma, xpa,
     1  d1mach
      external d1mach
      data eps / 0.d0 /
c
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
      real*8 function d9lgit (a, x, algap1)
      implicit integer*4 (i-n)

c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute the log of tricomi-s incomplete gamma function with
c     perron-s continued fraction for large x and for a .ge. x.
c
      real*8 a, x, algap1, ax, a1x, eps, fk, hstar, p, r, s,
     1  sqeps, t, d1mach
      external d1mach
      data eps, sqeps / 2*0.d0 /
c
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

      real*8 function d9lgmc (x)
      implicit integer*4 (i-n)

c     august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     compute the log gamma correction factor for x .ge. 10. so that
c     log (dgamma(x)) = log(dsqrt(2*pi)) + (x-.5)*log(x) - x +
c     d9lgmc(x)
c
      real*8 x, algmcs(15), xbig, xmax, dcsevl, d1mach
      external d1mach, dcsevl,   initds
c
c     series for algm       on the interval  0.          to  1.00000e-02

c                                        with weighted error   1.28e-31
c                                         log weighted error  30.89
c                               significant figures required  29.81
c                                    decimal places required  31.48
c
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
c
      data nalgm, xbig, xmax / 0, 2*0.d0 /
c
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
c
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
c
      end
      real*8 function d9pak (y, n)
      implicit integer*4 (i-n)

c     december 1979 edition. w. fullerton, c3, los alamos scientific lab.
c
c     pack a base 2 exponent into floating point number x.  this routine is
c     almost the inverse of d9upak.  it is not exactly the inverse, because
c     abs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
c     were known to be in range we could compute
c     d9pak = x * 2.0d0**n
c
      real*8 y, aln2b, aln210, d1mach
      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094 8873623478 7031942948 9 d0 /
c
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
      real*8 x, y, absx
c
      absx = abs(x)
      n = 0
      y = 0.0d0
      if (x.eq.0.0d0) return
c
 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10
c
 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20
c
 30   y = dsign (absx, x)
      return
c
      end
      real*8 function dcsevl (x, a, n)
      implicit integer*4 (i-n)
c
c     evaluate the n-term chebyshev series a at x.  adapted from
c     r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c     input arguments --
c       x      dble prec value at which the series is to be evaluated.
c       a      dble prec array of n terms of a chebyshev series.  in eval-
c              uating a, only half the first coef is summed.
c       n      number of terms in array a.
c
      real*8 a(n), x, twox, b0, b1, b2
c
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

      real*8 function dgamma (x)
      implicit integer*4 (i-n)
c     jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
c     jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml

      real*8 x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1  xmin, y, d9lgmc, dcsevl, d1mach, dint, xsml
      external d1mach, d9lgmc, dcsevl, dint,  initds
c
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
c
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /
c
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
c     1  31hdgamma  x is a negative integer, 31, 4, 2)
     1  "dgamma  x is a negative integer     ", it31, it4, it2)
      if (x.lt.(-0.5d0) .and. abs((x-dint(x-0.5d0))/x).lt.dxrel) call
cc     1  seteru (68hdgamma  answer lt half precision because x too near n
cc     2egative integer, 68, 1, 1)
     1  seteru ( "dgamma  answer lt half precision because x too near n
     2egative integer", it68, it1, it1)

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
c     1  61hdgamma  answer lt half precision, x too near negative integer
     1  "dgamma  answer lt half precision, x too near negative integer"
     2  , it61, it1, it1)
c
      sinpiy = sin (pi*y)
      if (sinpiy.eq.0.d0) call seteru (
c     1  31hdgamma  x is a negative integer, 31, 4, 2)
     1  "dgamma  x is a negative integer     ", it31, it4, it2)
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end

      real*8 function dgamr (x)
      implicit integer*4 (i-n)
c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c     this routine, not dgamma(x), should be the fundamental one.
c
      real*8 x, alngx, sgngx, dgamma, dint
      external  dgamma, dint
c
      dgamr = 0.0d0
      if (x.le.0.0d0 .and. dint(x).eq.x) return
c
      it1=1
      call entsrc (irold, it1)
      if (abs(x).gt.10.0d0) go to 10
      dgamr = 1.0d0/dgamma(x)
      call erroff
      call entsrc (ir, irold)
      return
c
 10   call dlgams (x, alngx, sgngx)
      call erroff
      call entsrc (ir, irold)
      dgamr = sgngx * exp(-alngx)
      return
c
      end

      real*8 function dint (x)
      implicit integer*4 (i-n)

c     december 1983 edition. w. fullerton, c3, los alamos scientific lab.
c
c     dint is the real*8 equivalent of aint.  this portable
c     version is quite efficient when the argument is reasonably small (a
c     common case), and so no faster machine-dependent version is needed.
c     
      real*8 x, xscl, scale, xbig, xmax, part, d1mach, r1mach
      external d1mach,  i1mach, r1mach
      data npart, scale, xbig, xmax / 0, 3*0.0d0 /
c
      if (npart.ne.0) go to 10

      it4=4
      it9=9
      it10=10

      ibase = i1mach(it10)
      xmax = 1.0d0/d1mach(it4)
      xbig = dmin1 (dble (i1mach(it9)), 1.d0/r1mach(it4))
      scale = ibase**int(log(xbig)/log(dble(float(ibase)))-0.5d0)
      npart = log(xmax)/log(scale) + 1.0d0
c
 10   if (x.lt.(-xbig) .or. x.gt.xbig) go to 20
c
      dint = int(sngl(x))
      return
c
 20   xscl = abs(x)
      if (xscl.gt.xmax) go to 50
c
      do 30 i=1,npart
        xscl = xscl/scale
 30   continue
c
      dint = 0.0d0
      do 40 i=1,npart
        xscl = xscl*scale
        ipart = xscl
        part = ipart
        xscl = xscl - part
        dint = dint*scale + part
 40   continue
c
      if (x.lt.0.0d0) dint = -dint
      return
c
c$$$ 50   call seteru (68hdint    abs(x) may be too big to be represented a
c$$$     1s an exact integer, 68, 1, 1)
      it1=1
 1    it68=68
 50   call seteru ("dint abs(x) may be too big to be represented a
     1s an exact integer", it68, it1, it1)
      dint = x
      return
c
      end
      subroutine dlgams (x, dlgam, sgngam)
      implicit integer*4 (i-n)
c     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c     evaluate log abs (gamma(x)) and return the sign of gamma(x) in sgngam.
c     sgngam is either +1.0 or -1.0.
c     
      real*8 x, dlgam, sgngam, dint, dlngam
      external dint, dlngam
c
      dlgam = dlngam(x)
      sgngam = 1.0d0
      if (x.gt.0.d0) return
c
      int = dmod (-dint(x), 2.0d0) + 0.1d0
      if (int.eq.0) sgngam = -1.0d0
c
      return
      end

      real*8 function dlngam (x)
      implicit integer*4 (i-n)

     
c     august 1980 edition.   w. fullerton, c3, los alamos scientific lab.
      real*8 x, dxrel, pi, sinpiy, sqpi2l, sq2pil,
     1  y, xmax, dint, dgamma, d9lgmc, d1mach
      external d1mach, d9lgmc, dgamma, dint
c
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
c     sq2pil = alog (sqrt(2*pi)),  sqpi2l = alog(sqrt(pi/2))
      data sqpi2l / +.2257913526 4472743236 3097614947 441 d+0    /
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c
      data xmax, dxrel / 2*0.d0 /
c
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
c     1  31hdlngam  x is a negative integer, 31, 3, 2)
     1  "dlngam  x is a negative integer     ",it31,it3,it2)
c
      dlngam = sqpi2l + (x-0.5d0)*log(y) - x - log(sinpiy) - d9lgmc(y)
c
      if (abs((x-dint(x-0.5d0))*dlngam/x).lt.dxrel) call seteru (
c     1  68hdlngam  answer lt half precision because x too near negative
     1 "dlngam  answer lt half precision because x too near negative
     2integer", it68, it1, it1)
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
c       1 + 71/(the number of characters stored per integer word).
c
c     integer messgp(36),fmt(14),ccplus
      character*1 messgp(36),fmt(14),ccplus
c
c     start with no previous message.
c
      data messgp(1)/'1'/, nwp/0/, nerrp/0/
c
c     set up the format for printing the error message.  the format is
c     simply (a1,14x,72axx) where xx=i1mach(6) is the number of
c     characters stored per integer word.
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
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
        it0=0
        it1=1
 20   if (i8save(it1,it0,.false.).eq.0) go to 30
c
c  print the message.
c
      it4=4
        iwunit=i1mach(it4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
c
        it2=2
        it6=6
        call s88fmt(it2,i1mach(it6),fmt(12))
c        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
        print *, ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end
      subroutine entsrc(irold,irnew)
      implicit integer*4 (i-n)
c
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
c
      irold=i8save(it2,irnew,irnew.ne.0)
c
c  if have an error state, stop execution.
c
      if (i8save(it1,it0,.false.) .ne. 0) call seterr
c     1   (39hentsrc - called while in an error state,39,2,2)
     1   ("entsrc - called while in an error state",
     2    it39,it2,it2)
c
      return
c
      end
      subroutine eprint
      implicit integer*4 (i-n)
c
c     this subroutine prints the last error message, if any.
c
c      parameter (messglength=100)
      parameter (messglength=36)
      character*1 messg(messglength)
      it1=1
      call e9rint(messg,it1,it1,.false.)
      return
c
      end
      subroutine erroff
      implicit integer*4 (i-n)
c
c     turns off the error state off by setting lerror=0.
c
      external i8save
c
      it0=0
      it1=1
      i=i8save(it1,it0,.true.)
      return
c
      end
      integer*4 function i8save(isw,ivalue,set)
      implicit integer*4 (i-n)
c
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
c     
c     start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end
      integer*4 function initds (dos, nos, eta)
      implicit integer*4 (i-n)

c     june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c     initialize the real*8 orthogonal series dos so that
c     initds is the number of terms needed to insure the error is no
c     larger than eta.  ordinarily eta will be chosen to be one-tenth
c     machine precision.
c
c     input arguments --
c     dos    dble prec array of nos coefficients in an orthogonal series.
c     nos    number of coefficients in dos.
c     eta    requested accuracy of series.
c
      real*8 dos(nos)
c
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
c
c 20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28,
 20   if (i.eq.nos) call seteru ("initds  eta may be too small        ",
     1     it28, it1, it2)
      initds = i
c
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
c
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
c
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
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
        i6=6
        i72=72
 10     nw=(min0(nmessg,i72)-1)/i1mach(i6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        it1=1
        itemp=i8save(it1,it1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(it1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      it0=0
      it2=2
      if (i8save(it2,it0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
 60   call fdump
      stop 'seterr'
c
      end
      subroutine seteru (messg, nmessg, nerr, iopt)
      implicit integer*4 (i-n)
c      parameter (messglength=100)
      parameter (messglength=36)
      common /cseter/ iunflo
      character*1 messg(messglength)
      data iunflo / 0 /
c
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      it1=1
      call seterr (messg, nmessg, nerr, it1)
c
      return
      end

      integer*4 FUNCTION I1MACH(I)
      implicit integer*4 (i-n)

C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
C    INTEGERS HAVE FORM SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C    I1MACH( 7) = A, THE BASE.
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C    FLOATS HAVE FORM  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C               WHERE  EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, THE BASE.
C  SINGLE-PRECISION
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C  DOUBLE-PRECISION
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
      INTEGER*4 CRAY1, IMACH(16), OUTPUT, SANITY, SMALL(2)
      COMMON /D8MACH/ CRAY1
      SAVE IMACH, SANITY
      REAL*8 RMACH
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
      IF (SANITY .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532
     *           .AND. SMALL(2) .EQ. -448790528)
     *       .OR.     (SMALL(2) .EQ. 1117925532
     *           .AND. SMALL(1) .EQ. -448790528)) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -128
               IMACH(13) = 127
               IMACH(14) = 53
               IMACH(15) = -1024
               IMACH(16) = 1023
               SANITY = 987
            ELSE
cc               CRAY1 = 4617762693716115456.D0
               CRAY1 = 16115456.D0
               IF (SMALL(1) .NE. CRAY1) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
               IMACH(1) = 5
               IMACH(2) = 6
               IMACH(3) = 102
               IMACH(4) = 6
               IMACH(5) = 64
               IMACH(6) = 8
               IMACH(7) = 2
               IMACH(8) = 63
cc               IMACH(9) = 372036854775807.D0
               IMACH(9) = 54775807.D0
               IMACH(10) = 2
               IMACH(11) = 47
               IMACH(12) = -8189
               IMACH(13) = 8190
               IMACH(14) = 94
               IMACH(15) = -8099
               IMACH(16) = 8190
               SANITY = 987
               GO TO 10
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
         SANITY = 987
         END IF
 9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/
     * ' statements appropriate for your machine and setting'/
     * ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
 9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/
     * ' appropriate for your machine.')
 10   IF (I .LT. 1  .OR.  I .GT. 16) GO TO 30
      I1MACH = IMACH(I)
C REMOVE THE FOLLOWING LINE IF FORTRAN66 IS PREFERRED TO FORTRAN77.
      IF (I .EQ. 6) I1MACH = 1
      RETURN
 30   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
      STOP
      end


      subroutine fdump
      return
      end

      REAL*8 FUNCTION D1MACH(I)
      implicit integer*4 (i-n)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.  IF YOU DO NOT
C  KNOW WHICH SET TO USE, TRY BOTH AND SEE WHICH GIVES PLAUSIBLE
C  VALUES.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR D1MACH.
C
      INTEGER*4 SMALL(2)
      INTEGER*4 LARGE(2)
      INTEGER*4 RIGHT(2)
      INTEGER*4 DIVER(2)
      INTEGER*4 LOG10(2)
      INTEGER*4 SC
C/6S
C/7S
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
C/
      REAL*8 DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR BIG-ENDIAN IEEE ARITHMETIC (BINARY FORMAT)
C     MACHINES IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST,
C     SUCH AS THE AT&T 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G.
C     SUN 3), AND MACHINES THAT USE SPARC, HP, OR IBM RISC CHIPS.
C
c      DATA SMALL(1),SMALL(2) /    1048576,          0 /
c      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
c      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
c      DATA DIVER(1),DIVER(2) / 1018167296,          0 /
c      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR LITTLE-ENDIAN (BINARY) IEEE ARITHMETIC
C     MACHINES IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST,
C     E.G. IBM PCS AND OTHER MACHINES THAT USE INTEL 80X87 OR DEC
C     ALPHA CHIPS.
C
      DATA SMALL(1),SMALL(2) /          0,    1048576 /
      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /, SC/987/
C

C  ***  ISSUE STOP 779 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE
               WRITE(*,*)'Adjust D1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine.'
            STOP 779
            END IF
         SC = 987
         END IF
C
C  ***  ISSUE STOP 778 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...
      IF (DMACH(4) .GE. 1.0D0) STOP 778
*C/6S
*C     IF (I .LT. 1  .OR.  I .GT. 5)
*C    1   CALL SETERR(24HD1MACH - I OUT OF BOUNDS,24,1,2)
*C/7S
*      IF (I .LT. 1  .OR.  I .GT. 5)
*     1   CALL SETERR('D1MACH - I OUT OF BOUNDS',24,1,2)
*C/
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
C
* /* C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return DBL_MIN;
*	  case 2: return DBL_MAX;
*	  case 3: return DBL_EPSILON/FLT_RADIX;
*	  case 4: return DBL_EPSILON;
*	  case 5: return log10(FLT_RADIX);
*	  }
*
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1);
*	return 0; /* for compilers that complain of missing return values */
*	}
      END

      REAL*8 FUNCTION R1MACH(I)
      implicit integer*4 (i-n)
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR R1MACH.
C
      INTEGER*4 SMALL(2)
      INTEGER*4 LARGE(2)
      INTEGER*4 RIGHT(2)
      INTEGER*4 DIVER(2)
      INTEGER*4 LOG10(2)
      INTEGER*4 SC
C/6S
C/7S
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
C/
      REAL*8 RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C      DATA SMALL(1) /     8388608 /
C      DATA LARGE(1) /  2139095039 /
C      DATA RIGHT(1) /   864026624 /
C      DATA DIVER(1) /   872415232 /
C      DATA LOG10(1) /  1050288283 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00564000000000000000B /
C      DATA RMACH(2) / 37767777777777777776B /
C      DATA RMACH(3) / 16414000000000000000B /
C      DATA RMACH(4) / 16424000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / O"00564000000000000000" /
C      DATA RMACH(2) / O"37767777777777777776" /
C      DATA RMACH(3) / O"16414000000000000000" /
C      DATA RMACH(4) / O"16424000000000000000" /
C      DATA RMACH(5) / O"17164642023241175720" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /, SC/987/
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /, SC/987/
C
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH(1) = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      SMALL(1) .EQ. 1117925532
     *          .AND. SMALL(2) .EQ. -448790528) THEN
*              *** IEEE BIG ENDIAN ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2146435071
               LARGE(2) = -1
               RIGHT(1) = 1017118720
               RIGHT(2) = 0
               DIVER(1) = 1018167296
               DIVER(2) = 0
               LOG10(1) = 1070810131
               LOG10(2) = 1352628735
            ELSE IF ( SMALL(2) .EQ. 1117925532
     *          .AND. SMALL(1) .EQ. -448790528) THEN
*              *** IEEE LITTLE ENDIAN ***
               SMALL(2) = 1048576
               SMALL(1) = 0
               LARGE(2) = 2146435071
               LARGE(1) = -1
               RIGHT(2) = 1017118720
               RIGHT(1) = 0
               DIVER(2) = 1018167296
               DIVER(1) = 0
               LOG10(2) = 1070810131
               LOG10(1) = 1352628735
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*              *** VAX WITH D_FLOATING ***
               SMALL(1) = 128
               SMALL(2) = 0
               LARGE(1) = -32769
               LARGE(2) = -1
               RIGHT(1) = 9344
               RIGHT(2) = 0
               DIVER(1) = 9472
               DIVER(2) = 0
               LOG10(1) = 546979738
               LOG10(2) = -805796613
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2147483647
               LARGE(2) = -1
               RIGHT(1) = 856686592
               RIGHT(2) = 0
               DIVER(1) = 873463808
               DIVER(2) = 0
               LOG10(1) = 1091781651
               LOG10(2) = 1352628735
            ELSE
               WRITE(*,*)'Adjust autodoubled R1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine.'
               STOP 777
               END IF
         ELSE
            RMACH(1) = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*              *** IEEE ***
               SMALL(1) = 8388608
               LARGE(1) = 2139095039
               RIGHT(1) = 864026624
               DIVER(1) = 872415232
               LOG10(1) = 1050288283
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*              *** VAX ***
               SMALL(1) = 128
               LARGE(1) = -32769
               RIGHT(1) = 13440
               DIVER(1) = 13568
               LOG10(1) = 547045274
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*              *** IBM ***
               SMALL(1) = 1048576
               LARGE(1) = 2147483647
               RIGHT(1) = 990904320
               DIVER(1) = 1007681536
               LOG10(1) = 1091781651
            ELSE
               WRITE(*,*)'Adjust R1MACH by uncommenting'
               WRITE(*,*)'data statements appropriate for your machine.'
               STOP 777
               END IF
            END IF
         SC = 987
         END IF
C
C  ***  ISSUE STOP 776 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...
      IF (RMACH(4) .GE. 1.0) STOP 776
*C/6S
*C     IF (I .LT. 1  .OR.  I .GT. 5)
*C    1   CALL SETERR(24HR1MACH - I OUT OF BOUNDS,24,1,2)
*C/7S
*      IF (I .LT. 1  .OR.  I .GT. 5)
*     1   CALL SETERR('R1MACH - I OUT OF BOUNDS',24,1,2)
*C/
C
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'R1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      R1MACH = RMACH(I)

      RETURN
C
* /* C source for R1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*
*float r1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return FLT_MIN;
*	  case 2: return FLT_MAX;
*	  case 3: return FLT_EPSILON/FLT_RADIX;
*	  case 4: return FLT_EPSILON;
*	  case 5: return log10(FLT_RADIX);
*	  }
*
*	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
*	exit(1);
*	return 0; /* for compilers that complain of missing return values */
*	}
      END
