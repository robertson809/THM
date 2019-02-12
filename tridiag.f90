! This is the tridiagonal matrix polynomial program
! use lapack for any actually linear algebra
! lapak written in the 70s, core basis of matlab

! A bit tautologic here, define dp as the kind that a double is
! this just takes whatever this system calls a double and then
! makes dp that
integer, parameter  :: dp = kind(1.d0)
module tridiag
	implicit none
	Complex (kind = dp), Allocatable :: diag(:), upper(:), lower(:)
	

end module tridiag