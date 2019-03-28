program tridiag_driver
    use tridiag
    !loop counter
    integer :: i
    !tridiagonal matrix polynomial
    type(trid) :: a, b, c, d, e, f, g
    type(trid) :: output(3)
    type(trid_mp)  ::mp, mp2
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
    !testing scalar multiplication
    b%upper = a%upper * k
    b%diag = a%diag * k
    b%lower = a%lower * k
    
    !testing vector addition
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
    
    
    !Derivative tester
    !Create the MP version of (x-1)^3
	!allocating arrays
    allocate(d%diag(3), d%upper(2), d%lower(2))
    
	!putting numbers in the arrays
    do i = 1, 2
        d%diag(i) = 1
        d%upper(i) = 1
        d%lower(i) = 1
    end do
    d%diag(3) = 1
    
    e%upper = d%upper * -3
    e%diag = d%diag * -3
    e%lower = d%lower * -3
    
    f%upper = e%upper * -1
    f%diag = e%diag * -1
    f%lower = e%lower * -1
    
    g%upper = d%upper * -1
    g%diag = d%diag * -1
    g%lower = d%lower * -1
    
    allocate(mp2%coef(4))
    
    !put coefficients in the mp
    mp%coef(1) = g
    mp%coef(2) = f
    mp%coef(3) = e
    mp%coef(4) = d
    mp%size = 3
    mp%degree = 3
    
    !test Horner's method
    output = triHorner((1.0D+00,0.0D+00), mp2)
    call print_trid(output(1))
    call print_trid(output(2))
    call print_trid(output(3))
    write(*, '(A)') 'That was mp 2'
        
    deallocate(a%diag, a%upper, a%lower)
    deallocate(mp%coef)
end program tridiag_driver