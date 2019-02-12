program addNumbers

! This simple program adds two numbers
   implicit none

! Type declarations
   real :: a, b, result
   character(len = 80) :: message
   integer :: test

! Executable statements
   a = 12.0
   b = 15.0
   result = a + b + 1
   message = "A kindred soul discovered in black sand and steel"
   print *, message
   
   
   !do loop is really a while loop
   !stop is an exit statement, literally ends the program
   !cycle is a continue
   !exit exits the loop, equivalent to break
   do test = 1, 10
   		print *, "Iteration"
		if (test == 5) then
			exit
		end if
		print *, test
   end do
   print *, "I made it to line 25"
end program addNumbers