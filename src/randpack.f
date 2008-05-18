      REAL FUNCTION ranf()
C**********************************************************************
C
C     REAL FUNCTION RANF()
C                RANDom number generator as a Function
C
C     Returns a random floating point number from a uniform distribution
C     over 0 - 1 (endpoints of this interval are not returned) using the
C     current generator
C
C     This is a transcription from Pascal to Fortran of routine
C     Uniform_01 from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C**********************************************************************
C     .. External Functions ..
      INTEGER ignlgi
      EXTERNAL ignlgi
C     ..
C     .. Executable Statements ..
C
C     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
C      and is currently 2147483563. If M1 changes, change this also.
C
      ranf = ignlgi()*4.656613057E-10
      RETURN

      END

      INTEGER FUNCTION ignlgi()
C**********************************************************************
C
C     INTEGER FUNCTION IGNLGI()
C               GeNerate LarGe Integer
C
C     Returns a random integer following a uniform distribution over
C     (1, 2147483562) using the current generator.
C
C     This is a transcription from Pascal to Fortran of routine
C     Random from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER curntg,k,s1,s2,z
      LOGICAL qqssd
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn,inrgcm,rgnqsd,setall
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C
C     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
C     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
C     THIS ROUTINE  2) A CALL TO SETALL.
C
      IF (.NOT. (qrgnin())) CALL inrgcm()
      CALL rgnqsd(qqssd)
      IF (.NOT. (qqssd)) CALL setall(1234567890,123456789)
C
C     Get Current Generator
C
      CALL getcgn(curntg)
      s1 = cg1(curntg)
      s2 = cg2(curntg)
      k = s1/53668
      s1 = a1* (s1-k*53668) - k*12211
      IF (s1.LT.0) s1 = s1 + m1
      k = s2/52774
      s2 = a2* (s2-k*52774) - k*3791
      IF (s2.LT.0) s2 = s2 + m2
      cg1(curntg) = s1
      cg2(curntg) = s2
      z = s1 - s2
      IF (z.LT.1) z = z + m1 - 1
      IF (qanti(curntg)) z = m1 - z
      ignlgi = z
      RETURN

      END

      LOGICAL FUNCTION qrgnin()
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNIN()
C               Q Random GeNerators INitialized?
C
C     A trivial routine to determine whether or not the random
C     number generator has been initialized.  Returns .TRUE. if
C     it has, else .FALSE.
C
C**********************************************************************
C     .. Scalar Arguments ..
      LOGICAL qvalue
C     ..
C     .. Local Scalars ..
      LOGICAL qinit
C     ..
C     .. Entry Points ..
      LOGICAL qrgnsn
C     ..
C     .. Save statement ..
      SAVE qinit
C     ..
C     .. Data statements ..
      DATA qinit/.FALSE./
C     ..
C     .. Executable Statements ..
      qrgnin = qinit
      RETURN

      ENTRY qrgnsn(qvalue)
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNSN( QVALUE )
C               Q Random GeNerators Set whether iNitialized
C
C     Sets state of whether random number generator is initialized
C     to QVALUE.
C
C     This routine is actually an entry in QRGNIN, hence it is a
C     logical function.  It returns the (meaningless) value .TRUE.
C
C**********************************************************************
      qinit = qvalue
      qrgnsn = .TRUE.
      RETURN

      END

      SUBROUTINE getcgn(g)
      INTEGER g
C**********************************************************************
C
C      SUBROUTINE GETCGN(G)
C                         Get GeNerator
C
C     Returns in G the number of the current random number generator
C
C
C                              Arguments
C
C
C     G <-- Number of the current random number generator (1..32)
C                    INTEGER G
C
C**********************************************************************
C
      INTEGER curntg,numg
      SAVE curntg
      PARAMETER (numg=32)
      DATA curntg/1/
C
      g = curntg
      RETURN

      ENTRY setcgn(g)
C**********************************************************************
C
C     SUBROUTINE SETCGN( G )
C                      Set GeNerator
C
C     Sets  the  current  generator to G.    All references to a generat
C     are to the current generator.
C
C
C                              Arguments
C
C
C     G --> Number of the current random number generator (1..32)
C                    INTEGER G
C
C**********************************************************************
C
C     Abort if generator number out of range
C
      IF (.NOT. (g.LT.0.OR.g.GT.numg)) GO TO 10
      WRITE (*,*) ' Generator number out of range in SETCGN:',
     +  ' Legal range is 1 to ',numg,' -- ABORT!'
      STOP ' Generator number out of range in SETCGN'

   10 curntg = g
      RETURN

      END

      SUBROUTINE inrgcm()
C**********************************************************************
C
C     SUBROUTINE INRGCM()
C          INitialize Random number Generator CoMmon
C
C
C                              Function
C
C
C     Initializes common area  for random number  generator.  This saves
C     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
C     assuring that the routine is loaded with the other routines.
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER i
      LOGICAL qdum
C     ..
C     .. External Functions ..
      LOGICAL qrgnsn
      EXTERNAL qrgnsn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     V=20;                            W=30;
C
C     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
C     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
C
C   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
C    An efficient way to precompute a**(2*j) MOD m is to start with
C    a and square it j times modulo m using the function MLTMOD.
C
      m1 = 2147483563
      m2 = 2147483399
      a1 = 40014
      a2 = 40692
      a1w = 1033780774
      a2w = 1494757890
      a1vw = 2082007225
      a2vw = 784306273
      DO 10,i = 1,numg
          qanti(i) = .FALSE.
   10 CONTINUE
C
C     Tell the world that common has been initialized
C
      qdum = qrgnsn(.TRUE.)
      RETURN

      END
      
      SUBROUTINE setall(iseed1,iseed2)
C**********************************************************************
C
C      SUBROUTINE SETALL(ISEED1,ISEED2)
C               SET ALL random number generators
C
C     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
C     initial seeds of the other generators are set accordingly, and
C     all generators states are set to these seeds.
C
C     This is a transcription from Pascal to Fortran of routine
C     Set_Initial_Seed from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISEED1 -> First of two integer seeds
C                                   INTEGER ISEED1
C
C     ISEED2 -> Second of two integer seeds
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
      LOGICAL qssd
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g,ocgn
      LOGICAL qqssd
C     ..
C     .. External Functions ..
      INTEGER mltmod
      LOGICAL qrgnin
      EXTERNAL mltmod,qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn,initgn,inrgcm,setcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/,qqssd
C     ..
C     .. Data statements ..
      DATA qqssd/.FALSE./
C     ..
C     .. Executable Statements ..
C
C     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
C      HAS BEEN CALLED.
C
      qqssd = .TRUE.
      CALL getcgn(ocgn)
C
C     Initialize Common Block if Necessary
C
      IF (.NOT. (qrgnin())) CALL inrgcm()
      ig1(1) = iseed1
      ig2(1) = iseed2
      CALL initgn(-1)
      DO 10,g = 2,numg
          ig1(g) = mltmod(a1vw,ig1(g-1),m1)
          ig2(g) = mltmod(a2vw,ig2(g-1),m2)
          CALL setcgn(g)
          CALL initgn(-1)
   10 CONTINUE
      CALL setcgn(ocgn)
      RETURN

      ENTRY rgnqsd(qssd)
C**********************************************************************
C
C     SUBROUTINE RGNQSD
C                    Random Number Generator Query SeeD set?
C
C     Returns (LOGICAL) QSSD as .TRUE. if SETALL has been invoked,
C     otherwise returns .FALSE.
C
C**********************************************************************
      qssd = qqssd
      RETURN

      END
 
      INTEGER FUNCTION mltmod(a,s,m)
C**********************************************************************
C
C     INTEGER FUNCTION MLTMOD(A,S,M)
C
C                    Returns (A*S) MOD M
C
C     This is a transcription from Pascal to Fortran of routine
C     MULtMod_Decompos from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     A, S, M  -->
C                         INTEGER A,S,M
C
C**********************************************************************
C     .. Parameters ..
      INTEGER h
      PARAMETER (h=32768)
C     ..
C     .. Scalar Arguments ..
      INTEGER a,m,s
C     ..
C     .. Local Scalars ..
      INTEGER a0,a1,k,p,q,qh,rh
C     ..
C     .. Executable Statements ..
C
C     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
C      machine. On a different machine recompute H
C
      IF (.NOT. (a.LE.0.OR.a.GE.m.OR.s.LE.0.OR.s.GE.m)) GO TO 10
      WRITE (*,*) ' A, M, S out of order in MLTMOD - ABORT!'
      WRITE (*,*) ' A = ',a,' S = ',s,' M = ',m
      WRITE (*,*) ' MLTMOD requires: 0 < A < M; 0 < S < M'
      STOP ' A, M, S out of order in MLTMOD - ABORT!'

   10 IF (.NOT. (a.LT.h)) GO TO 20
      a0 = a
      p = 0
      GO TO 120

   20 a1 = a/h
      a0 = a - h*a1
      qh = m/h
      rh = m - h*qh
      IF (.NOT. (a1.GE.h)) GO TO 50
      a1 = a1 - h
      k = s/qh
      p = h* (s-k*qh) - k*rh
   30 IF (.NOT. (p.LT.0)) GO TO 40
      p = p + m
      GO TO 30

   40 GO TO 60

   50 p = 0
C
C     P = (A2*S*H)MOD M
C
   60 IF (.NOT. (a1.NE.0)) GO TO 90
      q = m/a1
      k = s/q
      p = p - k* (m-a1*q)
      IF (p.GT.0) p = p - m
      p = p + a1* (s-k*q)
   70 IF (.NOT. (p.LT.0)) GO TO 80
      p = p + m
      GO TO 70

   80 CONTINUE
   90 k = p/qh
C
C     P = ((A2*H + A1)*S)MOD M
C
      p = h* (p-k*qh) - k*rh
  100 IF (.NOT. (p.LT.0)) GO TO 110
      p = p + m
      GO TO 100

  110 CONTINUE
  120 IF (.NOT. (a0.NE.0)) GO TO 150
C
C     P = ((A2*H + A1)*H*S)MOD M
C
      q = m/a0
      k = s/q
      p = p - k* (m-a0*q)
      IF (p.GT.0) p = p - m
      p = p + a0* (s-k*q)
  130 IF (.NOT. (p.LT.0)) GO TO 140
      p = p + m
      GO TO 130

  140 CONTINUE
  150 mltmod = p
C
      RETURN

      END
      
      SUBROUTINE initgn(isdtyp)
C**********************************************************************
C
C     SUBROUTINE INITGN(ISDTYP)
C          INIT-ialize current G-e-N-erator
C
C     Reinitializes the state of the current generator
C
C     This is a transcription from Pascal to Fortran of routine
C     Init_Generator from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISDTYP -> The state to which the generator is to be set
C
C          ISDTYP = -1  => sets the seeds to their initial value
C          ISDTYP =  0  => sets the seeds to the first value of
C                          the current block
C          ISDTYP =  1  => sets the seeds to the first value of
C                          the next block
C
C                                   INTEGER ISDTYP
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER isdtyp
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      INTEGER mltmod
      EXTERNAL qrgnin,mltmod
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' INITGN called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' INITGN called before random number generator initialized'

   10 CALL getcgn(g)
      IF ((-1).NE. (isdtyp)) GO TO 20
      lg1(g) = ig1(g)
      lg2(g) = ig2(g)
      GO TO 50

   20 IF ((0).NE. (isdtyp)) GO TO 30
      CONTINUE
      GO TO 50
C     do nothing
   30 IF ((1).NE. (isdtyp)) GO TO 40
      lg1(g) = mltmod(a1w,lg1(g),m1)
      lg2(g) = mltmod(a2w,lg2(g),m2)
      GO TO 50

   40 STOP 'ISDTYP NOT IN RANGE'

   50 cg1(g) = lg1(g)
      cg2(g) = lg2(g)
      RETURN

      END
     