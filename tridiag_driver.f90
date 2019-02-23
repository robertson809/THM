program tridiag_driver
    use tridiag
    !loop counter
    integer :: i 
    !tridiagonal matrix polynomial
    type(trid) :: a, b
	
	!allocating arrays
    allocate(a%diag(5), a%upper(4), a%lower(4))
	
	!putting numbers in the arrays
    do i = 1, 4
        a%diag(i) = i
        a%upper(i) = i
        a%lower(i) = i
    end do
        
    a%diag(5) = 5
    call print_trid(a)
    
    write(*,'(A)') 'That was A.'
    
    a = s_mult(a, 5) 
    call print_trid(a)
    
    write(*,'(A)') 'That was A * 5.'

    b = m_add(a, a)
    call print_trid(b)
    
    write(*,'(A)') 'That was 10A.'
    
    deallocate(a%diag, a%upper, a%lower)
    
end program tridiag_driver