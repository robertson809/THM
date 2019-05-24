!********************************************************************************
!   Tridiag_driver: Driver Program for tridiag
!   Author: Michael Robertson, Thomas R. Cameron, Davidson College
!   Last Modified: 20 May 2019
!********************************************************************************
! MIT License
!
! Copyright (c) 2019 Michael Robertson, Thomas Cameron
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!********************************************************************************
!   This module contains the following paramaters, functions, and 
!   subroutines:
!       dp: integer paramater for machine-compiler specific double 
!       precision.
!
!       eps: machine-compiler specific double precision unit roundoff.
!
!       main: main subroutine for computing roots of polynomial.
!
!       triHorner: function for evaluating a tridiagonal matrix polynomial using 
!                   an adaptation of Horner's method
!
!       triBack: subroutine for solving a special upper tridiagonal linear system
!                  using back substitution
! 
!********************************************************************************
program tridiag_driver
    use tridiag
    ! loop counter
    integer :: i
    ! tridiagonal matrix polynomial
    type(trid) :: a_3, a_2, a_1, a_0, b_0, b_1, b_2, b_3, c, c1, d, d1, eval
    complex(kind=dp), allocatable :: y(:), y2(:), v(:)
    type(trid) :: output(3)
    type(trid_mp)  ::mp, mp2, mp3
    complex(kind=dp)        :: k
    
    
    ! Derivative tester
    ! create the MP version of (x-1)^3
    allocate(a_3%diag(3), a_3%upper(2), a_3%lower(2))
    allocate(a_2%diag(3), a_2%upper(2), a_2%lower(2))
    allocate(a_1%diag(3), a_1%upper(2), a_1%lower(2))
    allocate(a_0%diag(3), a_0%upper(2), a_0%lower(2))
    allocate(mp%coef(4))

    
	! putting numbers in the arrays
    write(*, '(A)') 'line65'
    ! (x-1)^3 = x^3 -3x^2 +3x - 1
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
    
    
    ! Horner's tester
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
    
    ! Back substitution tester
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
    !N({{3,4,5,0},{0,6,7,8},{0,0,9,10},{0,0,0,12}}^-1 {{0},{0},{11},{13}})
    
    n = size(eval%diag)
    allocate(y(n-1))
    y = cmplx(0,0,kind = dp)
    y(n-1) = eval%diag(n)
    y(n-2) = eval%upper(n-1) !from Hyman's decomposition

    write(*, '(A)') 'Calling back substitution on MP(3)'
    call triBack(eval%lower, eval%diag(2:n-1), eval%upper(2:n-2), y) !from Hyman's decomposition
    write(*, '(A)') 'x, from R(u)x = y from Hymans'
    do i = 1, size(y)
        print *, y(i)
    end do
    
    
    
    ! Proper Hyman tester
    allocate(mp2%coef(3))
    a_0%lower(1) = 3
    a_0%diag(1) = 1
    a_0%upper(1) = 2
    
    a_0%lower(2) = 6
    a_0%diag(2) = 4
    a_0%upper(2) = 5
    
    
    a_1%lower(1) = 7 
    a_1%diag(1) = 8
    a_1%upper(1) = 9
    
    a_1%lower(2) = 10 
    a_1%diag(2) = 11
    a_1%upper(2) = 12
        
    a_2%lower(1) = 13
    a_2%diag(1) = 14
    a_2%upper(1) = 15
    
    a_2%lower(2) = 16
    a_2%diag(2) = 17
    a_2%upper(2) = 18
                  
    a_0%diag(3) = 19
    a_1%diag(3) = 20
    a_2%diag(3) = 21
    a_3%diag(3) = 22
    
    mp2%coef(1) = a_0
    mp2%coef(2) = a_1
    mp2%coef(3) = a_2
    mp2%size = 3
    mp2%degree = 2
    
    write(*, '(A)') 'Printing for HYMANS'
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
    
    write(*, '(A)') 'Printing MP2 for proper Hymans'
    call print_trid_mp(mp2)
    write(*, '(A)') new_line('A')
    
    !evaluate MP
    output = triHorner(cmplx(3,0, kind = dp), mp2)
    
    write(*, '(A)') 'MP2(3)'
    call print_trid(output(1))
    write(*, '(A)') new_line('A_2')
    
    write(*, '(A)') 'MP2''(3)'
    call print_trid(output(2))
    write(*, '(A)') new_line('A_2')
    
    ! Proper Hyman tester number 2
    allocate(b_3%diag(4), b_3%upper(3), b_3%lower(3))
    allocate(b_2%diag(4), b_2%upper(3), b_2%lower(3))
    allocate(b_1%diag(4), b_1%upper(3), b_1%lower(3))
    allocate(b_0%diag(4), b_0%upper(3), b_0%lower(3))
    allocate(mp3%coef(4))
    
    b_0%lower(1) = 4
    b_0%lower(2) = 15
    b_0%lower(3) = 1
    
    b_0%diag(1) = 1
    b_0%diag(2) = 7
    b_0%diag(3) = 19
    b_0%diag(4) = 6
    
    b_0%upper(1) = 6
    b_0%upper(2) = 4
    b_0%upper(3) = 23

    b_1%lower(1) = 3
    b_1%lower(2) = 14
    b_1%lower(3) = 4
    
    b_1%diag(1) = 4
    b_1%diag(2) = 6
    b_1%diag(3) = 18
    b_1%diag(4) = 2
    
    b_1%upper(1) = 2
    b_1%upper(2) = 3
    b_1%upper(3) = 22
    
    b_2%lower(1) = 2
    b_2%lower(2) = 13
    b_2%lower(3) = 1
    
    b_2%diag(1) = 1
    b_2%diag(2) = 5
    b_2%diag(3) = 17
    b_2%diag(4) = 9
    
    b_2%upper(1) = 9
    b_2%upper(2) = 2
    b_2%upper(3) = 21
    
    b_3%lower(1) = 1
    b_3%lower(2) = 12
    b_3%lower(3) = 3
    
    b_3%diag(1) = 3
    b_3%diag(2) = 5
    b_3%diag(3) = 16
    b_3%diag(4) = 5
    
    b_3%upper(1) = 5
    b_3%upper(2) = 1
    b_3%upper(3) = 20
    
    mp3%coef(1) = b_0
    mp3%coef(2) = b_1
    mp3%coef(3) = b_2
    mp3%coef(4) = b_3
    mp3%size = 4
    mp3%degree = 3
    
    write(*, '(A)') 'Printing for HYMANS'
    write(*, '(A)') 'Printing B_0'
    call print_trid(b_0)
    write(*, '(A)') new_line('B_0')
    
    write(*, '(A)') 'Printing B_1'
    call print_trid(b_1)
    write(*, '(A)') new_line('B_1')
    
    write(*, '(A)') 'Printing B_2'
    call print_trid(b_2)
    write(*, '(A)') new_line('B_2')
    
    write(*, '(A)') 'Printing B_3'
    call print_trid(b_3)
    write(*, '(A)') new_line('B_3')
    
    write(*, '(A)') 'Printing MP3 for proper Hymans'
    call print_trid_mp(mp3)
    write(*, '(A)') new_line('A')
    
    ! Proper Hyman tester number 3
    allocate(c%diag(3), c%upper(2), c%lower(2))
    allocate(c1%diag(3), c1%upper(2), c1%lower(2))
    
    c%lower(1) = 58
    c%lower(2) = 441
    c%diag(1) = 103
    c%diag(2) = 205
    c%diag(3) = 658
    c%upper(1) = 228
    c%upper(2) = 58
    
    c1%lower(1) = 42
    c1%lower(2) = 402
    c1%diag(1) = 91
    c1%diag(2) = 171
    c1%diag(3) = 552
    c1%upper(1) = 191
    c1%upper(2) = 42
    
    write(*, '(A)') 'Printing for HYMANS three dimensional test case'
    write(*, '(A)') 'Printing C'
    call print_trid(c)
    write(*, '(A)') new_line('B_0')
    
    write(*, '(A)') 'Printing c prime '
    call print_trid(c1)
    write(*, '(A)') new_line('B_1')
    
    ! Spedial 2x2 case
    allocate(d%diag(2), d%upper(1), d%lower(1))
    allocate(d1%diag(2), d1%upper(1), d1%lower(1))
    
    d%lower(1) = 58
    d%diag(1) = 103
    d%diag(2) = 205
    d%upper(1) = 228
    
    d1%lower(1) = 42
    d1%diag(1) = 91
    d1%diag(2) = 171
    d1%upper(1) = 191

    
   !evaluate MP
!     output = triHorner(cmplx(3,0, kind = dp), mp3)
    
!     write(*, '(A)') 'MP3(3)'
!     call print_trid(output(1))
!     write(*, '(A)') new_line('A_2')
!
!     write(*, '(A)') 'MP3''(3)'
!     call print_trid(output(2))
!     write(*, '(A)') new_line('A_2')
    
    !main Hyman tester, with output which is three trids, the value, and first two derivatives
!     k = PHyman(output)
!     write(*, '(A)') new_line('A_2')
!     print *, 'Step from proper Hyman is', k
!     write(*, '(A)') new_line('A_2')

    
    
    !Three dimensional test case
    write(*, '(A)') 'Printing for HYMANS three dimensional test case'
    write(*, '(A)') 'Printing C'
    call print_trid(c)
    write(*, '(A)') new_line('B_0')
    
    write(*, '(A)') 'Printing c prime '
    call print_trid(c1)
    write(*, '(A)') new_line('B_1')
    
    output(1) = c
    output(2) = c1
    k = SHyman3(output, output(1)%lower, output(1)%diag, output(1)%upper, output(2)%lower, &
    output(2)%diag, output(2)%upper)
    
    write(*, '(A)') new_line('A_2')
    print *, 'Step from Hyman 3 is', k
    print *, 'Error', abs(k) - 3.1249929625482458
    if (abs(k) - 3.1249929625482458 < 0.0000009) then
        print *, 'Error', abs(k) - 3.1249929625482458, 'is CORRECT'
    end if
    write(*, '(A)') new_line('A_2')
    
    !Two dimensional test case
 !    write(*, '(A)') 'Printing for HYMANS two dimensional test case'
!     write(*, '(A)') 'Printing D'
!     call print_trid(d)
!     write(*, '(A)') new_line('B_0')
!
!     write(*, '(A)') 'Printing D prime '
!     call print_trid(d1)
!     write(*, '(A)') new_line('B_1')
    
    output(1) = d
    output(2) = d1
    
    k = SHyman2(output(1)%lower, output(1)%diag, output(1)%upper, output(2)%lower, &
    output(2)%diag, output(2)%upper)
    
    write(*, '(A)') new_line('A_2')
    print *, 'Step from Hyman 2 is', k
    print *, 'Error:', k - 1.9787099226967431
    write(*, '(A)') new_line('A_2')
    
        
    deallocate(a_0%lower, a_0%diag, a_0%upper, a_1%lower, a_1%diag, a_1%upper, a_2%lower, a_2%diag, a_2%upper, &
    a_3%lower, a_3%diag, a_3%upper, mp%coef, y, eval%lower, eval%diag, eval%upper)
end program tridiag_driver