program tridiag_driver
    use tridiag
    !loop counter and input
    integer :: i, j, k, deg, size
    !tridiagonal matrix polynomial
    type(trid) :: a_3, a_2, a_1, a_0
    type(trid) :: output(3)
    type(trid_mp)     ::mp
    character(len=32)    :: arg, poly_file
    integer :: test
    complex(kind=dp)        :: aux
    
    ! load polynomial
    read(arg, '(A)') poly_file
    !don't really get what's going on here
    write(*,*) 'Starting: Running on '//poly_file
    !goes to poly_file, but defaults to data_files if "problem" with poly_file
    !unit = 1 is the number associated with the file
    open(unit=1,file="testTRID.txt", status='OLD')
    read(1,*) deg
    read(1,*) size 
    print *, deg, size
    allocate(mp%coef(deg))
    do i = 1, deg + 1
        allocate(mp%coef(i)%diag(size), mp%coef(i)%upper(size - 1), mp%coef(i)%lower(size - 1))
    end do 
        
    !allocate(p(deg+1), roots(deg), berr(deg), cond(deg), conv(deg))
    !read through every matrix
    do i = 1, deg + 1
        !upper
        print * , "ehllo"
        do k=1, size -1
            read(1,*) test
            print *, test
            !mp%coef(i)%upper(k)
        end do
        !diag
        do k=1, size
            read(1,*) mp%coef(i)%diag(k)
        end do
        !lower
        do k=1, size -1
            read(1,*) mp%coef(i)%lower(k)   
        end do
    end do
    
    
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
    call print_trid(a_0)
    write(*, '(A)') new_line('A_0')
    
    write(*, '(A)') 'Printing A_1'
    call print_trid(a_1)
    write(*, '(A)') new_line('A_1')
    
    write(*, '(A)') 'Printing A_2'
    call print_trid(a_2)
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'Printing A_2'
    call print_trid(a_3)
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'Printing MP'
    call print_trid_mp(mp)
    write(*, '(A)') new_line('A')
    
    
    !test Horner's method
    write(*, '(A)') 'Calling Horners on the MP'
    output = triHorner(cmplx(1,0, kind = dp), mp)
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'MP(1)'
    call print_trid(output(1))
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') '1st Derivative(1)'
    call print_trid(output(2))
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') '2nd Derivative(1)'
    call print_trid(output(3))
    write(*, '(A)') new_line('A_2')

        
    deallocate(a_0%lower, a_0%diag, a_0%upper, a_1%lower, a_1%diag, a_1%upper, a_2%lower, a_2%diag, a_2%upper, &
    a_3%lower, a_3%diag, a_3%upper, mp%coef)
end program tridiag_driver