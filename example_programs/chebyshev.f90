    Program chebyshev

!     ==========================================================================
!     Compare the performance of C02AFF with Thomas R. Cameron's 'FPML' when
!     finding the roots of the Chebyshev polynomials.
!     ==========================================================================

!     .. Use Statements ..
      Use compare_utils, Only: c02aff_helper, fpml_helper,                     &
        compute_coefficients, get_max_rel_err, quick_c02aff_call
      use nag_library, Only: nag_wp
!     .. Implicit None Statement ..
      Implicit None
!     .. Parameters ..
      Integer, Parameter :: n1 = 5, n2 = 80, nstep = 1, reps = 100,            &
        zero = 0.0_nag_wp, one = 1.0_nag_wp, two = 2.0_nag_wp,                 &
        nag_csv = 7, fpml_csv = 8, fpmlp_csv = 9, fpmlcp_csv = 10
!     .. Local Scalars ..
      Integer :: i, n, polish
      Real (Kind=nag_wp) :: fpmlerr, fpmltime, nagerr, nagtime, pi, tmp, eps
!     .. Local Arrays ..
      Complex (Kind=nag_wp), Allocatable :: a(:), as(:,:), b(:), c(:), z(:), zout(:)
      Real (Kind=nag_wp), Allocatable :: cond(:)
!     .. Intrinsic Procedures ..
      Intrinsic :: abs, cmplx, real, sum
!     .. Executable Statements ..

      Open (unit=nag_csv, file='data_files/chebyshev/c02aff.csv')
      Open (unit=fpml_csv, file='data_files/chebyshev/fpml.csv')
      Open (unit=fpmlp_csv, file='data_files/chebyshev/fpml_p.csv')
      Open (unit=fpmlcp_csv, file='data_files/chebyshev/fpml_cp.csv')

      Do i = nag_csv, fpmlcp_csv
        Write (i,*) 'Degree, Maximum Relative Error, Average Run Time'
      End Do

      eps = epsilon(1.0_nag_wp)/2.0_nag_wp
      pi = 3.1415926535897932384626433832795028841971693993751058209750E0_nag_wp

      Allocate (as(0:n2,0:n2))
      as(0:n2,0:n2) = (zero, zero)
      as(0,n2) = (one, zero)
      as(1,n2-1) = (one, zero)

!     Get Chebysev coefficients up to n2: as(n) = 2*shift(as(n-1)) - as(n-2) ==> a = b - c
      Allocate (a(0:n2), b(0:n2), c(0:n2))
      Do n = 2, n2
        c(0:n2) = as(n-2, 0:n2)
        b(0:n2-1) = two * as(n-1, 1:n2)
        b(n2) = (zero, zero)
        a = b - c
        as(n, 0:n2) = a(0:n2)
      End Do
      Deallocate (a, b, c)

      Call quick_c02aff_call()

      Do n = n1, n2, nstep

        Allocate (a(0:n), cond(n), z(n), zout(n))

!       Get coefficients from as
        a(0:n) = as(n, (n2-n):n2)

!       Set roots
        Do i = 1, n
          tmp = cos(pi*(two*i - one)/(two*n))
          If (abs(tmp) < eps) tmp = zero
          z(i) = cmplx(tmp, zero, kind=nag_wp)
        End Do

!       NAG call
        Call c02aff_helper(a, n, zout, reps, nagtime)
        nagerr = get_max_rel_err(n, z, zout)
        Write (nag_csv, 99998) n, nagerr, nagtime

!       FPML calls
        Do polish = 0, 2
          Call fpml_helper(a, n, zout, cond, reps, fpmltime, polish)
          fpmlerr = get_max_rel_err(n, z, zout)
          Write (fpml_csv + polish, 99998) n, fpmlerr, fpmltime
        End Do

        Deallocate (a, cond, z, zout)

      End Do

      Do i = nag_csv, fpmlcp_csv
        Close (i)
      End Do

99999 Format (1X, A, I5)
99998 Format (1X, I5, 1P, 2(', ', E12.4))

    End Program chebyshev
