program tridiag_driver
    use tridiag
    !loop counter
    integer :: i 
    !tridiagonal matrix polynomial
    type(trid) :: a
	
	!allocating arrays
    allocate(a%diag(10), a%upper(9), a%lower(9))
	
	!putting numbers in the arrays
    do i = 1, 9
        a%diag(i) = i
        a%upper(i) = i
        a%lower(i) = i
    end do
        
    a%diag(10) = 10
    call print_trid(a, 10)
    
    a = s_mult(a, 5) 
    call print_trid(a, 10)
    
    deallocate(a%diag, a%upper, a%lower)
    
end program tridiag_driver