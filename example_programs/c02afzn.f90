    Subroutine c02afzn(a,ndeg,scal,z,du,deflat,ierr)
!     MARK 14 RELEASE. NAG COPYRIGHT 1989.
!     MARK 15 REVISED. IER-891 (APR 1991).
!     MARK 15 REVISED. IER-941 (APR 1991).
!     MARK 17 REVISED. IER-1629 (JUN 1995).
!     MARK 20 REVISED. IER-3293 (JUL 2001).
!     THREAD-SAFETY (APRIL 2001)
!     REMOVED UNINITIALISED ACCESS TO OVFLOW AND UNFLOW ORIGINALLY
!     ASSUMED TO COME FROM COMMON BLOCK In MAIN PROGRAM.
!     BASED ON THE ROUTINE  ZERPOL, WRITTEN BY BRIAN T. SMITH

!     THIS SUBROUTINE COMPUTES THE N ZEROS OF THE COMPLEX POLYNOMIAL

!          N
!         SUM [A(1,k)+A(2,k)*I] * Z**(N-k) = 0
!         k=0

!     GIVEN BY CMPLX(Z(1,j),Z(2,j)), WHERE j = 1,2,...,N.

!     GAMA AND THETA ARE ARBITRARY PARAMETERS WHICH FOR THIS
!     IMPLEMENTATION HAVE BEEN SET TO 1.0 AND 2.0 RESPECTIVELY.
!     THERE IS NO INHERENT LIMITATION ON THE DEGREE OTHER THAN
!     AS THE DEGREE OF A POLYNOMIAL INCREASES, ITS ROOTS BECOME
!     ILL-CONDITIONED AS A FUNCTION OF THE COEFFICIENTS.

!     THIS PROGRAM IS DESIGNED TO TAKE ADVANTAGE OF SYSTEMS THAT REPORT
!     OVERFLOW CONDITIONS IN AN EFFICIENT WAY.  THAT IS, IF, WHENEVER
!     AN OVERFLOW OCCURS, A FLAG IS SET (THAT IS, THE LOGICAL VARIABLE
!     OVFLOW), C02AFXN CAN USE THIS INDICATOR TO INDICATE THAT THE ROOTS
!     OVERFLOW AND CANNOT BE REPRESENTED.

!     X02BJFN -- NUMBER OF DIGITS IN THE MANTISSA OF MODEL NUMBERS OF
!                TYPE DOUBLE PRECISION.
!     X02AJFN -- RELATIVE MACHINE PRECISION FOR ENTITIES OF TYPE
!                DOUBLE PRECISION.
!     X02ALFN -- LARGEST POSITIVE MACHINE REPRESENTABLE NUMBER OF TYPE
!                DOUBLE PRECISION.
!     X02BLFN -- MAXIMUM EXPONENT OF ENTITIES OF TYPE DOUBLE PRECISION.
!     X02AMFN -- SAFE RANGE PARAMETER FOR ENTITIES OF TYPE DOUBLE
!                PRECISION.
!     C02AGYN -- SCALE THE FIRST ARGUMENT BY A VALUE WITH AN EXPONENT
!                EQUAL TO THE SECOND ARGUMENT.
!     C02AGXN -- DETERMINE THE EXPONENT OF A NUMBER IN TERMS OF THE
!                MODEL.
!     C02AGSN -- RETURNS TRUE IF THE ARGUMENT IS UNNORMALIZED OR ZERO.
!     C02AFYN -- COMPUTE THE POLYNOMIAL VALUE, FIRST DERIVATIVE, AND
!                SECOND DERIVATIVE AT A COMPLEX POINT.
!     C02AFXN -- DETERMINE THE ROOTS OF A QUADRATIC EQUATION WITH
!                COMPLEX COEFFICIENTS.
!     C02AFWN -- PERFORMS THE COMPLEX DIVISION (CR,CI) = (AR,AI)/(BR,BI)
!                USING A MODIFIED VERSION OF NAG ROUTINE F06CLFN.

!     .. Use Statements ..
      Use nag_a02_ib_aux, Only: a02abfn, a02acfn
      Use nag_c02_ib_aux, Only: c02afwn, c02afxn, c02afyn, c02agsn, c02agxn,   &
                                c02agyn
      Use nag_precisions, Only: wp
!     .. Implicit None Statement ..
      Implicit None
!     .. Parameters ..
      Real (Kind=wp), Parameter        :: bigone = 1.0001E0_wp
      Real (Kind=wp), Parameter        :: gama = 1.0E0_wp
      Real (Kind=wp), Parameter        :: half = 0.5E0_wp
      Real (Kind=wp), Parameter        :: one = 1.0E0_wp
      Real (Kind=wp), Parameter        :: onepqt = 1.25E0_wp
      Real (Kind=wp), Parameter        :: rconst = 1.445E0_wp
      Real (Kind=wp), Parameter        :: small = 1.0E-3_wp
      Real (Kind=wp), Parameter        :: smlone = 0.99999E0_wp
      Real (Kind=wp), Parameter        :: theta = 2.0E0_wp
      Real (Kind=wp), Parameter        :: two = 2.0E0_wp
      Real (Kind=wp), Parameter        :: zero = 0.0E0_wp
!     .. Scalar Arguments ..
      Integer                          :: ierr
      Integer, Intent (In)             :: ndeg
      Logical                          :: scal
!     .. Array Arguments ..
      Real (Kind=wp), Intent (In)      :: a(2,0:ndeg)
      Real (Kind=wp)                   :: deflat(2,0:ndeg), du(2,0:ndeg),      &
                                          z(2,ndeg)
!     .. Local Scalars ..
      Real (Kind=wp)                   :: abdir, abdiro, abscl, deps, dpnewl,  &
                                          dpnewu, dx, dz0i, dz0r, dzni, dznr,  &
                                          e, f0, fact, fejer, finity, fn, g,   &
                                          lowerb, mxcoef, r, ratio, rtn, s, t, &
                                          temp, tiny, upperb, x2n, x2n1, xn,   &
                                          xn1, xn2, xn2n
      Integer                          :: dbase, emax, emin, expdep, i, ihalf, &
                                          ispir, iter, k, lrgexp, mnexp,       &
                                          mxcfex, mxexp, n, newl, newu, scbyex
      Logical                          :: cauchy, contin, ovf, ovflow, spiral, &
                                          startd, unf, unflow
!     .. Local Arrays ..
      Real (Kind=wp)                   :: c(2), cdir(2), cdiro(2), cf(2),      &
                                          cf1(2), cf2(2), cl(2), cr(2),        &
                                          cspir(2), ctemp(2)
!     .. Intrinsic Procedures ..
      Intrinsic                        :: abs, digits, epsilon, exp, log, max, &
                                          maxexponent, min, minexponent, real, &
                                          sqrt
!     .. Executable Statements ..
      Continue

      tiny = 0.22250738585072018772E-307_wp
      finity = one/tiny
      expdep = digits(one) + 1
      emin = minexponent(one) - 1
      emax = maxexponent(one) - 1
      lrgexp = emax + 1 - expdep
      deps = epsilon(one)/two
      ierr = 0
      iter = 0
      ihalf = 0
      ispir = 0
      n = ndeg

!     SET OVERFLOW/UNDERFLOW INDICATORS.

      ovf = .False.
      unf = .False.
      ovflow = .False.
      unflow = .False.

!     MOVE THE REAL AND IMAGINARY PARTS OF THE COEFFICIENTS TO DU(1,I)
!     AND DU(2,I) RESPECTIVELY AND DETERMINE THE LARGEST COEFFICIENT
!     IN MAGNITUDE.

      mxcoef = zero
      Do i = 0, n
        du(1,i) = a(1,i)
        du(2,i) = a(2,i)
        mxcoef = max(mxcoef,a02abfn(xr=a(1,i),xi=a(2,i)))
      End Do
      If (mxcoef==zero) Then
        Do i = 1, n
          z(1,i) = finity
          z(2,i) = zero
        End Do
        n = 0
        ovf = .True.
      Else

!       DETERMINE A SCALING FOR THE COEFFICIENTS SO THAT THE LARGEST
!       COEFFICIENT IN MAGNITUDE IS LARGE IN MAGNITUDE -- THAT IS, NEAR
!       DEPS * BASE ** X02BLFN(), UNLESS THE LARGEST COEFFICIENT IN
!       MAGNITUDE IS LARGER THAN THIS QUANTITY.  IN THIS CASE, SET
!       SCAL TO .FALSE. AND DO NOT SCAL THE COEFFICIENTS.

        mxcfex = c02agxn(dx=mxcoef,dpnewl=dpnewl,dpnewu=dpnewu,deps=deps,      &
          temp=temp,fact=fact,dbase=dbase,mnexp=mnexp,mxexp=mxexp,newl=newl,   &
          newu=newu)
        If (mxcfex>lrgexp) Then
          scbyex = 0
          scal = .False.
        Else
          scbyex = lrgexp - mxcfex
        End If
      End If

!     INDICATE THAT THE CAUCHY REGION CONTAINING THE SMALLEST ZEROS
!     OF THE CURRENT POLYNOMIAL HAS NOT BEEN COMPUTED.

      cauchy = .False.
      abdiro = zero
      lowerb = zero
      rtn = zero
      xn = zero
      xn1 = zero
      xn2n = zero
      x2n = zero
      x2n1 = zero
100   Continue

!     DO LOOP WHILE N>2.

      If (n>2) Then

!       IF SCAL = .TRUE., SCALE THE COEFFICIENTS SO THAT THE LARGEST
!       COEFFICIENT IN MAGNITUDE IS LARGE.

        If (scal) Then
          If (scbyex/=0) Then
            Do i = 0, n
              du(1,i) = c02agyn(dx=du(1,i),iexp=scbyex,dpnewl=dpnewl,          &
                fact=fact,dbase=dbase,newl=newl,newu=newu)
              du(2,i) = c02agyn(dx=du(2,i),iexp=scbyex,dpnewl=dpnewl,          &
                fact=fact,dbase=dbase,newl=newl,newu=newu)
            End Do
            scbyex = 0
          End If
        End If
        unf = unflow .Or. unf

!       FIND THE NUMBER I OF CONSECUTIVE LEADING COEFFICIENTS EQUAL TO
!       ZERO.

        Do i = 0, n - 1
          If (c02agsn(x=du(1,i)) .And. c02agsn(x=du(2,i))) Then

!           EACH VANISHED LEADING COEFFICIENT YIELDS AN INFINITE
!           ZERO.

            z(1,n-i) = finity
            z(2,n-i) = zero
          Else

!           EXIT THE LOOP ON I -- THE FIRST NON-ZERO COEFFICIENT IS
!           THE I-TH COEFFICIENT.

            Go To 110
          End If
        End Do
110     Continue
        If (i/=0) Then

!         SLIDE BACK THE COEFFICIENTS AND DECLARE OVERFLOW.

          Do k = i, n
            du(1,k-i) = du(1,k)
            du(2,k-i) = du(2,k)
          End Do
          n = n - i

!         GIVE AN ERROR RETURN IF THE VANISHING LEADING COEFFICIENTS
!         HAVE OCCURRED BECAUSE THE POLYNOMIAL HAS BEEN SCALED DOWN.
!         THIS CAN ONLY HAPPEN IF OVERFLOW WAS DETECTED DURING
!         THE COMPUTATION WITH THE ORIGINAL INPUT COEFFICIENTS.

          If (scbyex==-expdep) Then
            ierr = 3
            Go To 150
          End If

!         SIGNAL OVERFLOW AND CYCLE THE DO LOOP ON N.

          ovf = .True.
          Go To 100
        End If

!       FIND THE NUMBER I OF CONSECUTIVE TRAILING COEFFICIENTS EQUAL TO
!       ZERO.

        Do i = n, 1, -1
          If (c02agsn(x=du(1,i)) .And. c02agsn(x=du(2,i))) Then

!           EXTRACT ROOTS (IF ANY) FROM THE ORIGIN = (0., 0.)

            z(1,i) = zero
            z(2,i) = zero
          Else

!           EXIT THE LOOP ON I -- THE FIRST NON-ZERO COEFFICIENT IS
!           THE I-TH COEFFICIENT.

            Go To 120
          End If
        End Do
120     Continue
        If (i/=n) Then

!         REDUCE THE DEGREE BY THE NUMBER OF ZERO ROOTS AND
!         THEN CYCLE THE DO LOOP ON N.

          n = i
          Go To 100
        End If

!       INITIALIZE LOGICAL UNDERFLOW/OVERFLOW CONDITION STATUS
!       VARIABLES.

        ovflow = .False.
        unflow = .False.

!       HENCEFORTH  N .GT. 2,  Q(0) .NE. 0. AND  Q(N) .NE. 0.,
!       WHERE Q(K) = CMPLX(DU(1,K),DU(2,K)) FOR K = 0 AND K = N.
!       CHECK TO SEE WHETHER THE CAUCHY BOUNDS NEED TO BE COMPUTED.

        If (.Not. cauchy) Then

!         INITIALIZE SOME USEFUL CONSTANTS.

          xn = real(n,kind=wp)
          xn1 = real(n-1,kind=wp)
          xn2 = real(n-2,kind=wp)
          x2n = two/xn
          x2n1 = x2n/xn1
          xn2n = xn2/xn
          rtn = sqrt(xn)

!         CALCULATE  G, AN UPPER BOUND FOR THE SMALLEST ZERO.
!         START WITH  G = ABS( GEOMETRIC MEAN OF THE ZEROS).

          g = exp((log(a02abfn(xr=du(1,n),xi=du(2,n)))-log(                    &
            a02abfn(xr=du(1,0),xi=du(2,0))))/xn+small)

!         CALCULATE LAGUERRE-STEP  CDIR  AND  FEJER-BOUND, WHICH IS
!         AN UPPER BOUND FOR THE SMALLEST ZERO.
!         CALCULATION OF THE LAGUERRE STEP INVOLVES THE SQUARE OF
!         RECIPROCAL OF NEWTON'S STEP.  SINCE IT CAN EASILY OVERFLOW,
!         THE FEJER BOUND IS CALCULATED WITH NO SUCH OVERFLOWS AND THE
!         LAGUERRE STEP IS CALCULATED FROM IT.

          ovflow = .False.
          Call c02afwn(ari=du(1,n-1),aii=du(2,n-1),bri=du(1,n),bii=du(2,n),    &
            cr=cr(1),ci=cr(2),fail=ovflow)

!         IF OVFLOW, A ROOT OF POLYNOMIAL IS WITHIN
!         N * BASE ** (X02BKFN()-1) OF  ZERO.

          If (ovflow) Then

!           THUS, ASSUME A ROOT IS ZERO, BY ASSUMING CMPLX(DU(1,N),
!           DU(2,N)) IS ZERO.

            z(1,n) = zero
            z(2,n) = zero
            n = n - 1

!           CYCLE THE DO LOOP ON N.

            Go To 100
          End If

!         THE LAGUERRE STEP AND FEJER BOUNDS ARE COMPUTED FROM THE
!         SMALLER ROOT OF A QUADRATIC POLYNOMIAL.

          ctemp(1) = x2n1*du(1,n-2)
          ctemp(2) = x2n1*du(2,n-2)
          cf2(1) = x2n*du(1,n-1)
          cf2(2) = x2n*du(2,n-1)
          Call c02afxn(ar=ctemp(1),ai=ctemp(2),br=cf2(1),bi=cf2(2),cr=du(1,n), &
            ci=du(2,n),zsm=c,zlg=cf1,ovflow=ovflow,finity=finity,emax=emax,    &
            emin=emin,expdep=expdep)
          cr(1) = xn2n*cr(1)
          cr(2) = xn2n*cr(2)
          ctemp(1) = (c(1)*cr(1)-c(2)*cr(2)) + xn1

!         Multi-thread testing: cause the current thread to yield
!         to another thread. Please don't move this comment unless you
!         know why it's here.
!-MTTST  CALL NAGYLD()

          ctemp(2) = c(2)*cr(1) + c(1)*cr(2)
          Call a02acfn(xr=c(1),xi=c(2),yr=ctemp(1),yi=ctemp(2),zr=cdiro(1),    &
            zi=cdiro(2))
          abdiro = a02abfn(xr=cdiro(1),xi=cdiro(2))
          g = min(g,bigone*min(a02abfn(xr=c(1),xi=c(2)),rtn*abdiro))

!         CALCULATE THE CAUCHY-LOWER BOUND  R  FOR THE SMALLEST ZERO
!         BY SOLVING FOR THE ROOT R OF THE POLYNOMIAL EQUATION
!           ABS(Q(N)) = SUM( ABS(Q(I))*R**(N-I), I = 0, N-1 )
!         USING NEWTON'S METHOD, WHERE Q(I) = CMPLX(DU(1,I),DU(2,I)).

          r = g
          s = bigone*g
          unflow = .False.

          Do i = 0, n
            deflat(1,i) = a02abfn(xr=du(1,i),xi=du(2,i))
          End Do

!         NEWTON ITERATION LOOP FOR THE CAUCHY LOWER BOUND R.
!         TEST ON T INTRODUCED TO AVOID AN INFINITE LOOP THAT
!         OCCURRED WHEN FLOATING-POINT ROUNDING CONTROL WAS SET
!         TO 64 BITS ON A PC.

          t = one
130       Continue
          If (r<s .And. t>deps) Then
            t = deflat(1,0)
            s = zero
            ovflow = .False.
            Do i = 1, n - 1
              s = r*s + t
              t = r*t + deflat(1,i)
            End Do
            s = r*s + t

!           IT CAN BE PROVED THAT S CANNOT UNDERFLOW.

            t = (r*t-deflat(1,n))/s
            s = r
            r = r - t
            Go To 130
          End If

          If (ovflow) Then

!           THE COEFFICIENTS ARE TOO LARGE;  SCALE THEM DOWN AND
!           THEN CYCLE THE DO LOOP ON N.

            scbyex = -expdep
            Go To 100
          End If

!         ABS( SMALLEST ROOT ) < R/(2**(1/N) - 1 ) <  1.445*N*R.
!         THUS, 1.445*N*R IS ANOTHER UPPER BOUND AND THE CAUCHY BOUND
!         HAS BEEN COMPUTED, SO SET

          cauchy = .True.
          upperb = min(rconst*xn*r,g)
          lowerb = smlone*s
          unf = unflow .Or. unf
        End If

!       NOW   LOWERB < ABS( SMALLEST ZERO ) < UPPERB
!       INITIALIZE THE ITERATION TO BEGIN AT THE ORIGIN.
!       (IN THE CODE BELOW, F0 IS INITIALIZED BUT ITS VALUE NEVER
!       USEFULLY REFERENCED -- IT AVOIDS REFERENCE TO AN UNINITIALIZED
!       VARIABLE IN THE TEST TO ACCEPT THE NEXT ITERATE WHEN THE
!       ITERATION IS NOT STARTED.)

        fejer = upperb
        g = upperb
        cdir(1) = cdiro(1)
        cdir(2) = cdiro(2)
        abdir = abdiro
        ratio = abdir/g
        dznr = zero
        dzni = zero
        dz0r = zero
        dz0i = zero
        fn = a02abfn(xr=du(1,n),xi=du(2,n))
        f0 = fn
        spiral = .False.
        startd = .False.
        abscl = zero
        contin = .True.
140     Continue

!       DO WHILE (CONTIN) LOOP, SEARCHING FOR A REAL ROOT,
!       OR PAIR OF COMPLEX ROOTS.  THE NEXT ITERATE IS
!            ZN=CMPLX(DZNR , DZNI).

        If (contin) Then
          iter = iter + 1

!         RE-ENTRY POINT TO ACCEPT, MODIFY, OR REJECT THE LAGUERRE
!         STEP. REJECT  CDIR  IF  ABS(CDIR) > THETA*G .

          If (ratio>theta) Then

!           CURRENT LAGUERRE STEP IS NOT ACCEPTABLE.
!           IF STARTD, REDUCE PREVIOUS LAGUERRE STEP BY HALF.

            If (startd) Then
              ihalf = ihalf + 1
              abscl = half*abscl
              cl(1) = half*cl(1)
              cl(2) = half*cl(2)

!             HAS THE STEP BECOME NEGLIGIBLE?

              dx = abs(dznr) + abs(dzni)
              If (dx+abscl/=dx) Then
                dznr = dz0r + cl(1)
                dzni = dz0i + cl(2)
              Else

!               OTHERWISE, C02AFFN HAS HUNG-UP.

                If (fn>=e*xn**2) Then
                  ierr = 2
                  Go To 150
                End If

!               EXIT THE ITERATION LOOP  DO WHILE(CONTIN).

                contin = .False.
                Go To 140
              End If
            Else

!             IF .NOT. STARTD, HAS ZN BEEN ON THE INNER CAUCHY
!             RADIUS?

              ispir = ispir + 1
              If (spiral) Then
                c(1) = cspir(1)*dznr - cspir(2)*dzni
                c(2) = cspir(2)*dznr + cspir(1)*dzni
              Else

!               SET SPIRAL TO  .TRUE..  PUT  ZN  ON THE INNER
!               CIRCLE OF THE ANNULUS CONTAINING THE SMALLEST ZERO
!               IN THE DIRECTION OF THE LAGUERRE STEP.

                spiral = .True.
                cspir(1) = -onepqt/xn
                cspir(2) = one
                abscl = lowerb/xn**2
                ctemp(1) = cdir(1)/abdir
                ctemp(2) = cdir(2)/abdir
                c(1) = ctemp(1)*lowerb
                c(2) = ctemp(2)*lowerb
              End If

!             SET  ZN  TO THE NEXT POINT ON THE SPIRAL.

              dznr = c(1)
              dzni = c(2)
            End If
          Else

!           CDIR  AT THE ORIGIN IS IN THE DIRECTION OF DECREASING
!           FUNCTION VALUE, SO

            startd = .True.

!           ACCEPT  CDIR  IF  ABS(CDIR) <= GAMA*G.

            If (ratio>gama .And. (startd .Or. spiral .Or. lowerb<=gama*g))     &
              Then
              ratio = gama/ratio
              cdir(1) = cdir(1)*ratio
              cdir(2) = cdir(2)*ratio
              abdir = abdir*ratio
            End If

!           ACCEPT THE PREVIOUS ITERATE.  SAVE THE DATA ASSOCIATED
!           WITH THE CURRENT ITERATE.

            g = fejer
            cl(1) = cdir(1)
            cl(2) = cdir(2)
            abscl = abdir
            f0 = fn
            dz0r = dznr
            dz0i = dzni
            dznr = dz0r + cl(1)
            dzni = dz0i + cl(2)
          End If

!         BE SURE THAT THE OVERFLOW INDICATOR IS TURNED OFF.

          ovflow = .False.
          unflow = .False.

!         EVALUATE THE COMPLEX POLYNOMIAL AT A COMPLEX POINT.

          Call c02afyn(dx=dznr,dy=dzni,ndeg=n,a=du,p=cf,pprime=cf1,pdprim=cf2, &
            perr=e,deflat=deflat,ovflow=ovflow,deps=deps)
          fn = a02abfn(xr=cf(1),xi=cf(2))

!         CHECK FOR OVERFLOW.

          If (ovflow) Then

!           INDICATE THAT THE POLYNOMIAL NEEDS TO BE SCALED DOWN
!           AND CYCLE THE DO LOOP ON N.  NOTE: CAUCHY IS NOT RESET
!           AS THE CAUCHY BOUNDS NEED NOT BE RECOMPUTED.

            scbyex = -expdep
            Go To 100
          End If

!         CHECK TO SEE IF  ZN  IS A ZERO OR IF UNDERFLOW HAS OCCURRED.

          If (fn<=e .Or. unflow) Then
            If (unflow) Then
              ierr = 4
              unf = .True.
              Go To 150
            End If

!           A ROOT HAS BEEN FOUND -- EXIT THE ITERATION
!           LOOP  DO WHILE(CONTIN).

            contin = .False.
            Go To 140
          End If
          unf = unflow .Or. unf

!         HAS THE FUNCTION VALUE DECREASED?

          If (fn>=f0 .And. startd) Then

!           NO, IT HAS NOT.  INDICATE THAT THE LAGUERRE STEP IS
!           UNACCEPTABLE.  (A RATIO LARGER THAN THETA INDICATES THAT
!           THE LAGUERRE STEP SHOULD BE SHORTENED.)

            ratio = bigone*theta

!           CYCLE ITERATION LOOP  DO WHILE(CONTIN).

            Go To 140
          End If
          ovflow = .False.

!         FIND THE LAGUERRE STEP AT  ZN.

          Call c02afwn(ari=cf1(1),aii=cf1(2),bri=cf(1),bii=cf(2),cr=cr(1),     &
            ci=cr(2),fail=ovflow)

!         IF OVFLOW,  A ROOT OF POLYNOMIAL IS WITHIN
!         4 * N * BASE ** (X02BKFN()-1)  OF  ZN.

          If (ovflow) Then
            unf = .True.

!           A ROOT HAS BEEN FOUND -- EXIT THE ITERATION
!           LOOP  DO WHILE(CONTIN).

            contin = .False.
            Go To 140
          End If

!         COMPUTE THE LAGUERRE STEP  CDIR  AND THE BOUND  FEJER
!         AT  ZN .  THE LAGUERRE STEP AND FEJER BOUNDS ARE COMPUTED
!         FROM THE SMALLER ROOT OF A QUADRATIC POLYNOMIAL.

          cf2(1) = x2n1*cf2(1)
          cf2(2) = x2n1*cf2(2)
          ctemp(1) = x2n*cf1(1)
          ctemp(2) = x2n*cf1(2)
          Call c02afxn(ar=cf2(1),ai=cf2(2),br=ctemp(1),bi=ctemp(2),cr=cf(1),   &
            ci=cf(2),zsm=c,zlg=cf1,ovflow=ovflow,finity=finity,emax=emax,      &
            emin=emin,expdep=expdep)
          fejer = a02abfn(xr=c(1),xi=c(2))
          cr(1) = xn2n*cr(1)
          cr(2) = xn2n*cr(2)
          ctemp(1) = (c(1)*cr(1)-c(2)*cr(2)) + xn1
          ctemp(2) = c(2)*cr(1) + c(1)*cr(2)
          Call a02acfn(xr=c(1),xi=c(2),yr=ctemp(1),yi=ctemp(2),zr=cdir(1),     &
            zi=cdir(2))
          abdir = a02abfn(xr=cdir(1),xi=cdir(2))
          ratio = abdir/g
          fejer = min(rtn*abdir,fejer)

!         IS THE STEP SIZE NEGLIGIBLE?

          dx = abs(dznr) + abs(dzni)
          If (dx+abdir==dx) Then

!           THE STEP IS NEGLIGIBLE. ASSUME  ZN=(DZNR,DZNI) IS A ROOT.
!           EXIT THE ITERATION LOOP  DO WHILE(CONTIN).

            contin = .False.
            Go To 140
          End If

!         REPEAT THE ITERATION LOOP DO WHILE(CONTIN).

          Go To 140
        End If

!       A ROOT HAS BEEN COMPUTED.  DEFLATE THE POLYNOMIAL.

!       ACCEPT ZN AS A COMPLEX ROOT AND DEFLATE FOR A COMPLEX ROOT.
!       PUT COEFFICIENTS OF THE QUOTIENT POLYNOMIAL IN THE  DU  ARRAY.
!       DU(1,0) AND DU(2,0) ARE UNCHANGED FOR THE DEFLATED POLYNOMIAL.

        Do i = 1, n - 1
          du(1,i) = deflat(1,i)
          du(2,i) = deflat(2,i)
        End Do
        z(1,n) = dznr
        z(2,n) = dzni
        n = n - 1

!       INDICATE THAT THE CAUCHY REGION CONTAINING THE SMALLEST ZEROS
!       OF THE CURRENT POLYNOMIAL HAS NOT BEEN COMPUTED.

        cauchy = .False.

!       REPEAT THE LOOP WHILE N>2 FOR DECREASING N.

        Go To 100
      End If

!     THE POLYNOMIAL IS NOW OF DEGREE 2 OR LESS.  DETERMINE THE
!     REMAINING ROOTS DIRECTLY RATHER THAN ITERATIVELY.

      ovflow = .False.
      unflow = .False.
      If (n==2) Then
        Call c02afxn(ar=du(1,0),ai=du(2,0),br=du(1,1),bi=du(2,1),cr=du(1,2),   &
          ci=du(2,2),zsm=ctemp,zlg=c,ovflow=ovflow,finity=finity,emax=emax,    &
          emin=emin,expdep=expdep)
        z(1,1) = c(1)
        z(2,1) = c(2)
        z(1,2) = ctemp(1)
        z(2,2) = ctemp(2)
      Else If (n==1) Then
        Call a02acfn(xr=-du(1,1),xi=-du(2,1),yr=du(1,0),yi=du(2,0),zr=z(1,1),  &
          zi=z(2,1))
      Else
        ovf = ovf .Or. ovflow
        unf = unf .Or. unflow

!       PROVIDE ONLY THE RELEVANT OVER/UNDERFLOW MESSAGES.

        If (ovf) Then
          r = finity*finity
        End If
        If (unf) Then
          r = tiny*tiny
        End If
      End If
150   Continue
      Return
    End Subroutine c02afzn
