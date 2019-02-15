program addNumbers

   use Circle
   use tridiag
! This simple program adds two numbers
   implicit none
   
! Type declarations
   real :: a, b, result
   character(len = 80) :: message
   integer :: test
   print *, "I made it to line 24"
   print*, Pi
   print*, "ello"
   print *, "I made it to line 27"

! Executable statements
   a = 12.0
   b = 15.0
   result = a + b + 1
   message = "A kindred soul discovered in black sand and steel"
   print *, message
   
   
   !do loop is really a while loop
   !stop is an exit statement, literally ends the program
   !cycle is a continue
   !exit exits the loop, equivalent to breakgfo
   do test = 1, 10
   		print *, "Iteration"
		if (test == 5) then
			exit
		end if
		print *, test
   end do
   print *, "I made it to line 25"
end program addNumbers



! This is the tridiagonal matrix polynomial program

! A bit tautologic here, define dp as the kind that a double is
! this just takes whatever this system calls a double and then
! makes dp that

integer, parameter  :: dp = kind(1.d0)
module tridiag
	implicit none
	Complex (kind = dp), Allocatable :: m_diag(:), upper(:), lower(:)
	
contains      
   subroutine pretty_print()          
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
   end subroutine pretty_print
   
    
   subroutine randomize(size)          
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
   end subroutine randomize	

end module tridiag
