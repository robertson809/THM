program tridiag_driver
    use tridiag
    !loop counter
    integer :: i 
    !tridiagonal matrix polynomial
    type(trid) :: a, b, c
    type(trid_mp)  ::mp
	
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
    
    b = s_mult(a, 5 + 0i) 
    call print_trid(a)
    
    write(*,'(A)') 'That was A * 5.'

    c = m_add(a, b)
    call print_trid(b)
    
    write(*,'(A)') 'That was 10A.'
    
    
    
    !create a matrix polynomial
    
    allocate(mp%main(3))
    
    
    deallocate(a%diag, a%upper, a%lower)
end program tridiag_driver