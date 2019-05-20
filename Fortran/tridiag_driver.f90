program tridiag_driver
    use tridiag
    !loop counter
    integer :: i
    !tridiagonal matrix polynomial
    type(trid) :: a_3, a_2, a_1, a_0, eval
    complex(kind=dp), allocatable :: y(:)
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
    !(x-1)^3 = x^3 -3x^2 +3x - 1
    do i = 1, 2
        
        a_0%lower(i) = -1
        a_0%diag(i) = -1
        a_0%upper(i) = -1
       
        a_1%lower(i) = 3 
        a_1%diag(i) = 3
        a_1%upper(i) = 3
        
        a_2%lower(i) = -3
        a_2%diag(i) = -3
        a_2%upper(i) = -3
        
        a_3%lower(i) = 1
        a_3%diag(i) = 1
        a_3%upper(i) = 1
                 
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
    call print_trid(a_0)
    write(*, '(A)') new_line('A_0')
    
    write(*, '(A)') 'Printing A_1'
    call print_trid(a_1)
    write(*, '(A)') new_line('A_1')
    
    write(*, '(A)') 'Printing A_2'
    call print_trid(a_2)
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'Printing A_3'
    call print_trid(a_3)
    write(*, '(A)') new_line('A_3')
    
    write(*, '(A)') 'Printing MP'
    call print_trid_mp(mp)
    write(*, '(A)') new_line('A')
    
    
    !test Horner's method
    write(*, '(A)') 'Calling Horners on the MP'
    output = triHorner(cmplx(3,0, kind = dp), mp)
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'MP(3)'
    call print_trid(output(1))
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') '1st Derivative(3)'
    call print_trid(output(2))
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') '2nd Derivative(3)'
    call print_trid(output(3))
    write(*, '(A)') new_line('A_2')
    
    !make my test case
    allocate(eval%diag(5), eval%lower(4), eval%upper(4))
    eval%diag(1) = 1
    eval%diag(2) = 4
    eval%diag(3) = 7
    eval%diag(4) = 10
    eval%diag(5) = 13
    eval%lower(1) = 3
    eval%lower(2) = 6
    eval%lower(3) = 9
    eval%lower(4) = 12
    eval%upper(1) = 2
    eval%upper(2) = 5
    eval%upper(3) = 8
    eval%upper(4) = 11
    
    n = size(eval%diag)
    !test back substitution
    allocate(y(n-1)) 
    y = cmplx(0,0,kind = dp)
    y(n-1) = eval%diag(n)
    y(n-2) = eval%upper(n-1) !get y from the Horner's evaluation
    
    write(*, '(A)') 'Calling back substitution on MP(3)'
    call triBack(eval%lower, eval%diag(2:n-1), eval%upper(2:n-2), y)
    write(*, '(A)') 'x, from R(u)x = y from Hymans'
    do i = 1, size(y)
        write(*,'(F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i")') y(i)
    end do
    write(*, '(A)') new_line('A_2')
        
    deallocate(a_0%lower, a_0%diag, a_0%upper, a_1%lower, a_1%diag, a_1%upper, a_2%lower, a_2%diag, a_2%upper, &
    a_3%lower, a_3%diag, a_3%upper, mp%coef, y, eval%lower, eval%diag, eval%upper)
end program tridiag_driver