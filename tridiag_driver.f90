program tridiag_driver
    use tridiag
    !loop counter
    integer :: i
    !tridiagonal matrix polynomial
    type(trid) :: a_3, a_2, a_1, a_0
    type(trid) :: output(3)
    type(trid_mp)  ::mp
    complex(kind=dp)        :: k
    
    
    !Derivative tester
    !Create the MP version of (x-1)^3
	!allocating arrays
    allocate(a_3%diag(3), a_3%upper(2), a_3%lower(2))
    allocate(a_2%diag(3), a_2%upper(2), a_2%lower(2))
    allocate(a_1%diag(3), a_1%upper(2), a_1%lower(2))
    allocate(a_0%diag(3), a_0%upper(2), a_0%lower(2))
    allocate(mp%coef(4))

    
	!putting numbers in the arrays
    write(*, '(A)') 'line65'
    do i = 1, 2
        !(x-1)^3 = x^3 -3x^2 +3x - 1
        a_0%diag(i) = -1
        a_1%diag(i) = 3
        a_2%diag(i) = -3
        a_3%diag(i) = 1
        a_0%upper(i) = -1
        a_1%upper(i) = 3
        a_2%upper(i) = -3
        a_3%upper(i) = 1
        a_0%lower(i) = -1
        a_1%lower(i) = 3
        a_2%lower(i) = -3
        a_3%lower(i) = 1
    end do
    a_0%diag(3) = -1
    a_1%diag(3) = 3
    a_2%diag(3) = -3
    a_3%diag(3) = 1
    
    !put coefficients in the mp
    mp%coef(1) = a_0
    mp%coef(2) = a_1
    mp%coef(3) = a_2
    mp%coef(4) = a_3
    mp%size = 3
    mp%degree = 3
    
    
    write(*, '(A)') 'Printing A_0'
    call print_trid(a)
    write(*, '(A)') new_line('A_0')
    
    write(*, '(A)') 'Printing A_1'
    call print_trid(b)
    write(*, '(A)') new_line('A_1')
    
    write(*, '(A)') 'Printing A_2'
    call print_trid(c)
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'Printing MP'
    call print_trid_mp(mp)
    write(*, '(A)') new_line('A')
    
    
    !test Horner's method
    write(*, '(A)') 'That was mp 2'
    output = triHorner((1.0D+00,0.0D+00), mp2)
    call print_trid(output(1))
    call print_trid(output(2))
    call print_trid(output(3))
    write(*, '(A)') 'That was mp 2'
        
    deallocate(a%diag, a%upper, a%lower)
    deallocate(mp%coef)
end program tridiag_driver