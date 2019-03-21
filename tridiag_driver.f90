program tridiag_driver
    use tridiag
    !loop counter
    integer :: i
    !tridiagonal matrix polynomial
    type(trid) :: a, b, c, d
    type(trid_mp)  ::mp
    complex(kind=dp)        :: k
	
	!allocating arrays
    allocate(a%diag(5), a%upper(4), a%lower(4))
	
	!putting numbers in the arrays
    do i = 1, 4
        a%diag(i) = i
        a%upper(i) = i
        a%lower(i) = i
    end do
    a%diag(5) = 5 
    k = cmplx(5.0D+00,0.0D+00, dp)
    b%upper = a%upper * k
    b%diag = a%diag * k
    b%lower = a%lower * k
    
    c%upper = a%upper + b%upper
    c%diag = a%diag + b%diag
    c%lower = a%lower + a%lower
    !create a matrix polynomial
    allocate(mp%coef(3))
    !put coefficients in the mp
    mp%coef(1) = a
    mp%coef(2) = b
    mp%coef(3) = c
    mp%size = 5
    mp%degree = 2
    print *, mp%degree
    
    write(*, '(A)') 'Printing A'
    call print_trid(a)
    write(*, '(A)') new_line('A')
    
    write(*, '(A)') 'Printing B'
    call print_trid(b)
    write(*, '(A)') new_line('A')
    
    write(*, '(A)') 'Printing C'
    call print_trid(c)
    write(*, '(A)') new_line('A')
    
    write(*, '(A)') 'Printing MP'
    call print_trid_mp(MP)
    write(*, '(A)') new_line('A')
    
    
    !test Horner's method
    write(*, '(A)') 'Printing MP(3)'
    d = horner((3.0D+00,0.0D+00), mp)
    call print_trid(d)
    write(*, '(A)') 'That was mp(3)'
        
    deallocate(a%diag, a%upper, a%lower)
    deallocate(mp%coef)
end program tridiag_driver