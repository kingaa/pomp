      REAL FUNCTION gengam(a,r)
C**********************************************************************
C
C     REAL FUNCTION GENGAM( A, R )
C           GENerates random deviates from GAMma distribution
C
C
C                              Function
C
C
C     Generates random deviates from the gamma distribution whose
C     density is
C          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
C
C
C                              Arguments
C
C
C     JJV added the argument ranges supported
C     A --> Location parameter of Gamma distribution
C                              REAL A ( A > 0 )
C
C     R --> Shape parameter of Gamma distribution
C                              REAL R ( R > 0 )
C
C
C                              Method
C
C
C     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C               (Case R >= 1.0)
C               Ahrens, J.H. and Dieter, U.
C               Generating Gamma Variates by a
C               Modified Rejection Technique.
C               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
C     Algorithm GD
C
C     JJV altered the following to reflect sgamma argument ranges
C               (Case 0.0 < R < 1.0)
C               Ahrens, J.H. and Dieter, U.
C               Computer Methods for Sampling from Gamma,
C               Beta, Poisson and Binomial Distributions.
C               Computing, 12 (1974), 223-246/
C     Adapted algorithm GS.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL a,r
C     ..
C     .. External Functions ..
      REAL sgamma
      EXTERNAL sgamma
C     ..
C     .. Executable Statements ..

C     JJV added argument value checker
      IF ( a.GT.0.0 .AND. r.GT.0.0 ) GO TO 10
      WRITE (*,*) 'In GENGAM - Either (1) Location param A <= 0.0 or'
      WRITE (*,*) '(2) Shape param R <= 0.0 - ABORT!'
      WRITE (*,*) 'A value: ',a,'R value: ',r
      STOP 'Location or shape param out of range in GENGAM - ABORT!'
C     JJV end addition

 10   gengam = sgamma(r)/a
C      gengam = gengam/a
      RETURN

      END
      
      REAL FUNCTION sgamma(a)
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  G A M M A  DISTRIBUTION                             C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C               PARAMETER  A >= 1.0  !                                 C
C                                                                      C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               GENERATING GAMMA VARIATES BY A                         C
C               MODIFIED REJECTION TECHNIQUE.                          C
C               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  C
C                                                                      C
C     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     C
C                                 (STRAIGHTFORWARD IMPLEMENTATION)     C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C                                                                      C
C               PARAMETER  0.0 < A < 1.0  !                            C
C                                                                      C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              C
C               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              C
C               COMPUTING, 12 (1974), 223 - 246.                       C
C                                                                      C
C     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    C
C                                                                      C
C**********************************************************************C
C
C
C     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
C     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
C
C     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
C     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
C     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
C
C     .. Scalar Arguments ..
      REAL a
C     ..
C     .. Local Scalars .. (JJV added B0 to fix rare and subtle bug)
      REAL a1,a2,a3,a4,a5,a6,a7,aa,aaa,b,b0,c,d,e,e1,e2,e3,e4,e5,p,q,q0,
     +     q1,q2,q3,q4,q5,q6,q7,r,s,s2,si,sqrt32,t,u,v,w,x
C     ..
C     .. External Functions ..
      REAL ranf,sexpo,snorm
      EXTERNAL ranf,sexpo,snorm
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,alog,exp,sign,sqrt
C     ..
C     .. Save statement ..
C     JJV added Save statement for vars in Data satatements
      SAVE aa,aaa,s2,s,d,q0,b,si,c,q1,q2,q3,q4,q5,q6,q7,a1,a2,a3,a4,a5,
     +     a6,a7,e1,e2,e3,e4,e5,sqrt32
C     ..
C     .. Data statements ..
C
C     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
C     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
C
      DATA q1,q2,q3,q4,q5,q6,q7/.04166669,.02083148,.00801191,.00144121,
     +     -.00007388,.00024511,.00024240/
      DATA a1,a2,a3,a4,a5,a6,a7/.3333333,-.2500030,.2000062,-.1662921,
     +     .1423657,-.1367177,.1233795/
      DATA e1,e2,e3,e4,e5/1.,.4999897,.1668290,.0407753,.0102930/
      DATA aa/0.0/,aaa/0.0/,sqrt32/5.656854/
C     ..
C     .. Executable Statements ..
C
      IF (a.EQ.aa) GO TO 10
      IF (a.LT.1.0) GO TO 130
C
C     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
C
      aa = a
      s2 = a - 0.5
      s = sqrt(s2)
      d = sqrt32 - 12.0*s
C
C     STEP  2:  T=STANDARD NORMAL DEVIATE,
C               X=(S,1/2)-NORMAL DEVIATE.
C               IMMEDIATE ACCEPTANCE (I)
C
   10 t = snorm()
      x = s + 0.5*t
      sgamma = x*x
      IF (t.GE.0.0) RETURN
C
C     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
C
      u = ranf()
      IF (d*u.LE.t*t*t) RETURN
C
C     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
C
      IF (a.EQ.aaa) GO TO 40
      aaa = a
      r = 1.0/a
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
C
C               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
C               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
C               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
C
      IF (a.LE.3.686) GO TO 30
      IF (a.LE.13.022) GO TO 20
C
C               CASE 3:  A .GT. 13.022
C
      b = 1.77
      si = .75
      c = .1515/s
      GO TO 40
C
C               CASE 2:  3.686 .LT. A .LE. 13.022
C
   20 b = 1.654 + .0076*s2
      si = 1.68/s + .275
      c = .062/s + .024
      GO TO 40
C
C               CASE 1:  A .LE. 3.686
C
   30 b = .463 + s + .178*s2
      si = 1.235
      c = .195/s - .079 + .16*s
C
C     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
C
   40 IF (x.LE.0.0) GO TO 70
C
C     STEP  6:  CALCULATION OF V AND QUOTIENT Q
C
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 50
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 60

   50 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
C
C     STEP  7:  QUOTIENT ACCEPTANCE (Q)
C
   60 IF (alog(1.0-u).LE.q) RETURN
C
C     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
C               U= 0,1 -UNIFORM DEVIATE
C               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
C
   70 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = b + sign(si*e,u)
C
C     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
C
   80 IF (t.LT. (-.7187449)) GO TO 70
C
C     STEP 10:  CALCULATION OF V AND QUOTIENT Q
C
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 90
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 100

   90 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
C
C     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
C
  100 IF (q.LE.0.0) GO TO 70
      IF (q.LE.0.5) GO TO 110
C
C     JJV modified the code through line 125 to handle large Q case
C
      IF (q.LT.15.0) GO TO 105
C
C     JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
C     JJV so reformulate test at 120 in terms of one EXP, if not too big
C     JJV 87.49823 is close to the largest real which can be
C     JJV exponentiated (87.49823 = log(1.0E38))
C
      IF ((q+e-0.5*t*t).GT.87.49823) GO TO 125
      IF (c*abs(u).GT.exp(q+e-0.5*t*t)) GO TO 70
      GO TO 125

 105  w = exp(q) - 1.0
      GO TO 120

  110 w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
C
C               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
C
  120 IF (c*abs(u).GT.w*exp(e-0.5*t*t)) GO TO 70
 125  x = s + 0.5*t
      sgamma = x*x
      RETURN
C
C     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
C
C     JJV changed B to B0 (which was added to declarations for this)
C     JJV in 130 to END to fix rare and subtle bug.
C     JJV Line: '130 aa = 0.0' was removed (unnecessary, wasteful).
C     JJV Reasons: the state of AA only serves to tell the A .GE. 1.0
C     JJV case if certain A-dependant constants need to be recalculated.
C     JJV The A .LT. 1.0 case (here) no longer changes any of these, and
C     JJV the recalculation of B (which used to change with an
C     JJV A .LT. 1.0 call) is governed by the state of AAA anyway.
C
 130  b0 = 1.0 + .3678794*a
  140 p = b0*ranf()
      IF (p.GE.1.0) GO TO 150
      sgamma = exp(alog(p)/a)
      IF (sexpo().LT.sgamma) GO TO 140
      RETURN

  150 sgamma = -alog((b0-p)/a)
      IF (sexpo().LT. (1.0-a)*alog(sgamma)) GO TO 140
      RETURN

      END
     
      REAL FUNCTION sexpo()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM THE                 C
C               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  C
C               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       C
C     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
C
C     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
C     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
C
C     JJV added a Save statement for q (in Data statement)
C     .. Local Scalars ..
      REAL a,q1,u,umin,ustar
      INTEGER i
C     ..
C     .. Local Arrays ..
      REAL q(8)
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Equivalences ..
      EQUIVALENCE (q(1),q1)
C     ..
C     .. Save statement ..
      SAVE q
C     ..
C     .. Data statements ..
      DATA q/.6931472,.9333737,.9888778,.9984959,.9998293,.9999833,
     +     .9999986,.9999999/
C     ..
C
   10 a = 0.0
      u = ranf()
      GO TO 30

   20 a = a + q1
   30 u = u + u
C     JJV changed the following to reflect the true algorithm and
C     JJV prevent unpredictable behavior if U is initially 0.5.
C      IF (u.LE.1.0) GO TO 20
      IF (u.LT.1.0) GO TO 20
   40 u = u - 1.0
      IF (u.GT.q1) GO TO 60
   50 sexpo = a + u
      RETURN

   60 i = 1
      ustar = ranf()
      umin = ustar
   70 ustar = ranf()
      IF (ustar.LT.umin) umin = ustar
   80 i = i + 1
      IF (u.GT.q(i)) GO TO 70
   90 sexpo = a + umin*q1
      RETURN

      END

      REAL FUNCTION snorm()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
C
C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
C
C     .. Local Scalars ..
      REAL aa,s,tt,u,ustar,w,y
      INTEGER i
C     ..
C     .. Local Arrays ..
      REAL a(32),d(31),h(31),t(31)
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC float,int
C     ..
C     .. Save statement ..
C     JJV added a Save statement for arrays initialized in Data statmts
      SAVE a,d,t,h
C     ..
C     .. Data statements ..
      DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991,
     +     .2372021,.2776904,.3186394,.3601299,.4022501,.4450965,
     +     .4887764,.5334097,.5791322,.6260990,.6744898,.7245144,
     +     .7764218,.8305109,.8871466,.9467818,1.009990,1.077516,
     +     1.150349,1.229859,1.318011,1.417797,1.534121,1.675940,
     +     1.862732,2.153875/
      DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     +     .1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     +     .1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     +     .1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     +     .1134023,.1114027,.1095039/
      DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,
     +     .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,
     +     .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,
     +     .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,
     +     .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,
     +     .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     +     .5847031/
      DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,
     +     .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,
     +     .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,
     +     .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,
     +     .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,
     +     .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,
     +     .7010474/
C     ..
C     .. Executable Statements ..
C
   10 u = ranf()
      s = 0.0
      IF (u.GT.0.5) s = 1.0
      u = u + u - s
   20 u = 32.0*u
      i = int(u)
      IF (i.EQ.32) i = 31
      IF (i.EQ.0) GO TO 100
C
C                                START CENTER
C
   30 ustar = u - float(i)
      aa = a(i)
   40 IF (ustar.LE.t(i)) GO TO 60
      w = (ustar-t(i))*h(i)
C
C                                EXIT   (BOTH CASES)
C
   50 y = aa + w
      snorm = y
      IF (s.EQ.1.0) snorm = -y
      RETURN
C
C                                CENTER CONTINUED
C
   60 u = ranf()
      w = u* (a(i+1)-aa)
      tt = (0.5*w+aa)*w
      GO TO 80

   70 tt = u
      ustar = ranf()
   80 IF (ustar.GT.tt) GO TO 50
   90 u = ranf()
      IF (ustar.GE.u) GO TO 70
      ustar = ranf()
      GO TO 40
C
C                                START TAIL
C
  100 i = 6
      aa = a(32)
      GO TO 120

  110 aa = aa + d(i)
      i = i + 1
  120 u = u + u
      IF (u.LT.1.0) GO TO 110
  130 u = u - 1.0
  140 w = u*d(i)
      tt = (0.5*w+aa)*w
      GO TO 160

  150 tt = u
  160 ustar = ranf()
      IF (ustar.GT.tt) GO TO 50
  170 u = ranf()
      IF (ustar.GE.u) GO TO 150
      u = ranf()
      GO TO 140

      END
