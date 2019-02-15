! This is the tridiagonal matrix polynomial program
! use lapack for any actually linear algebra
! lapak written in the 70s, core basis of matlab

! A bit tautologic here, define dp as the kind that a double is
! this just takes whatever this system calls a double and then
! makes dp that
integer, parameter  :: dp = kind(1.d0)
module tridiag
	implicit none
	Complex (kind = dp), Allocatable :: m_diag(:), upper(:), lower(:)
	
	!allocate memory 
	allocate(m_diag)
	
	do i = 1, 5
		m_diag(i) = i + 1
	end do
	
	allocate(upper)
	allocate(lower)
	do i = 1, 4
		upper(i) = i + 1
		lower(i) = i + 1
		
	end do
	
	do i = size(lower)
		print *, "The entry for the upper diagonal is"
		print *, upper(i)
		print *, "The entry for the lower diagonal is"
		print *, lower(i)
		print *, "The entry for the middle diagonal is"	
		print *, m_diag(i)
		if i = size(lower) - 1 then
			print *, "The last entry in the main diagonal is"
			print *, m_diag(i + 1)
		end if

end module tridiag
