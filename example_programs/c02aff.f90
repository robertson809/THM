    Subroutine c02aff(a,n,scal,z,w,ifail)
!     MARK 14 RELEASE. NAG COPYRIGHT 1989.
!     MARK 15B REVISED. IER-944 (NOV 1991).

!     C02AFF ATTEMPTS TO FIND ALL THE ROOTS OF THE NTH ORDER COMPLEX
!     POLYNOMIAL EQUATION

!          N
!         SUM [A(1,k)+A(2,k)*I] * Z**(N-k) = 0.
!         k=0

!     THE ZEROS OF POLYNOMIALS OF DEGREE 1 AND 2 ARE CALCULATED BY
!     CAREFULLY EVALUATING THE "STANDARD" CLOSED FORMULAS
!            Z = -B/A AND
!            Z = (-B +/- SQRT(B*B-4*A*C))/(2*A) RESPECTIVELY, WHERE
!         A = CMPLX(A(1,0),A(2,0))
!         B = CMPLX(A(1,1),A(2,1)) AND
!         C = CMPLX(A(1,2),A(2,2)).
!     FOR N >= 3, THE ROOTS ARE LOCATED ITERATIVELY USING A VARIANT OF
!     LAGUERRE'S METHOD, WHICH IS CUBICALLY CONVERGENT FOR ISOLATED
!     ZEROS AND LINEARLY CONVERGENT FOR MULTIPLE ZEROS.

!     C02AFF ITSELF IS ESSENTIALLY A DUMMY ROUTINE WHOSE FUNCTION IS TO
!     PARTITION THE WORK ARRAY W FOR USE BY C02AFZN.
!     W IS PARTITIONED INTO 2 ARRAYS EACH OF SIZE 2*(N + 1).

!     .. Use Statements ..
      Use nag_c02_ib_aux, Only: c02afzn
      Use nag_precisions, Only: wp
!     .. Implicit None Statement ..
      Implicit None
!     .. Parameters ..
      Real (Kind=wp), Parameter        :: zero = 0.0E0_wp
      Character (*), Parameter         :: routine_name = 'C02AFFN'
!     .. Scalar Arguments ..
      Integer, Intent (Out)            :: ifail
      Integer, Intent (In)             :: n
      Logical, Intent (In)             :: scal
!     .. Array Arguments ..
      Real (Kind=wp), Intent (In)      :: a(2,0:n)
      Real (Kind=wp), Intent (Out)     :: w(4*(n+1)), z(2,n)
!     .. Local Scalars ..
      Real (Kind=wp)                   :: big
      Integer                          :: i, ierr, ndeg
      Logical                          :: sc
!     .. Intrinsic Procedures ..
      Intrinsic                        :: char, ichar, sqrt
!     .. Executable Statements ..
      Continue

      ierr = 0

      If (n<1) Then
        ierr = 1
      Else If ((a(1,0)==zero) .And. (a(2,0)==zero)) Then
        ierr = 1
      End If
      If (ierr/=0) Then
        Go To 100
      End If

!     End of argument checking.

!     Initialize Z to be -infinity.

      sc = scal
      ndeg = n
      big = 1.0E0_wp/(sqrt(2.0E0_wp)*0.22250738585072018772E-307_wp)
      Do i = 1, n
        z(1,i) = -big
        z(2,i) = -big
      End Do

      Call c02afzn(a=a,ndeg=ndeg,scal=sc,z=z,du=w(1),deflat=w(2*n+3),          &
        ierr=ierr)

      ifail = ierr

100   Continue
      Return

    End Subroutine c02aff
