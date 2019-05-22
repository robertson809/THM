!********************************************************************************
!   tridiag: Tridiagonal Matrix Polynomial solver
!   Author: Michael Robertson, Thomas Cameron, Davidson College
!   Last Modified: 5 May 2019
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
module tridiag
    implicit none
    integer, parameter  :: dp = kind(1.d0)
	
    !****************************************************************
    !				           Type trid                            *
    !                      Tridiagonal Matrix                       *
    !   Can contain values as a scalar or variable coefficients     *
    !****************************************************************
    type trid
        complex(kind=dp), allocatable :: diag(:), upper(:), lower(:)
    end type trid
    
    !**************************************************************
    !				           Type trid_mp		                  *
    !                Tridiagonal Matrix Polynomial	              *
    !   An array of variable matrix polynomials contains values   *
    ! for size and degree, which could be removed to save memory  *
    !**************************************************************
    ! This is a tridiagonal matrix polynoial
    type trid_mp
        type(trid), allocatable :: coef(:)
        integer                 :: size
        integer                 :: degree
    end type trid_mp
    
contains
    !**************************************************************
    !				Subroutine print_trid                         
    !   Prints a tridiagonal matrix diagonal by diagonal          
    !**************************************************************
    subroutine print_trid(a)
        implicit none
        ! argument variables
        type(trid)      :: a
        ! local variables
        integer         :: i
        integer         :: n
        
        n = size(a%diag)
			
        ! print statements
        write(*,'(F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i")') a%diag(1), a%upper(1)
        do i=2, n - 1
            write(*,'(F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i")') a%lower(i-1), a%diag(i), a%upper(i)
        end do
        write(*,'(F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i")') a%lower(n-1), a%diag(n)
    end subroutine print_trid
    
    !**************************************************************
    !				Subroutine print_trid_mp                       
    !   Prints a tridiagonal matrix polynomial as additions of 
    !   matrix polynomials        
    !**************************************************************
    subroutine print_trid_mp(a)
        implicit none
        !arguement variables
        type(trid_mp)       :: a
        !local variables
        integer             :: i
        integer             :: degree
        degree = size(a%coef)
        
        !printing
        do i = 1, degree 
            call print_trid(a%coef(i))
            !write (*,'(A, F8.3)') 'The answer is x = ', degree
            write(*,'(A, I1)') 'x^', i - 1
            if (i /= degree) then
                write(*, '(A)') '+'
            end if
        end do
    end subroutine print_trid_mp
    
    
    !*****************************************************************
    !                 Tridiagonal Matrix Multiplication                                  
    !*****************************************************************
    !    Multiplies a upper tridiagonal matrix A by a vector x
    !                                                                
    !  @param x - vector                                              
    !  @param diag - the diagonal
	!  @param upper - the 1 diagonal
	!  @param super - the 2 diagonal
    !  @return comp - Ax                                 
    !*****************************************************************
    function triMult(diag, upper, super, x) result(comp)
        implicit none
        !arguement variables
        complex(kind = dp)  :: x(:), diag(:), upper(:), super(:)
        !local variables
        integer             :: i, n
        !return variables
        complex(kind = dp), allocatable :: comp(:)
        
        n = size(diag)
        allocate(comp(n))
        
!Testing
!         print *, 'super has entries'
!         print *, super(1)
!         print *, 'but has size'
!         print *, size(super)
!
!         print *, 'printing diag'
!         do i = 1, size(diag)
!             print *, diag(i)
!         end do
!
!         print *, 'printing upper'
!         do i = 1, size(upper)
!             print *, upper(i)
!         end do
!
!         print *, 'printing super'
!         do i = 1, size(super)
!             print *, super(i)
!         end do
!
!         print *, 'printing mulitiplication terms'
!         do i = 1, size(super)
!             print *, super(i)
!         end do
        
        if(size(super) > 0) then
            print *, 'doing regular multiplication'
            !upper triangular tridiagonal matrix multiplication
            do i = 1, size(diag)
                comp(i) = diag(i) * x(i) + upper(i) * x(i+1) + super(i) * x(i + 2)
            end do
            comp(n) = diag(n) * x(n)
        else
            print *, 'doing 2x2 multiplication'
            !upper triangular tridiagonal matrix multiplication with 2x2 matrix
            do i = 1, size(diag) - 1
                comp(i) = diag(i) * x(i) + upper(i) * x(i+1)
            end do
            comp(n) = diag(n) * x(n)
        end if
        return
    end function
    
    
    !*****************************************************************
    !                       Reversal horner's method for             *
    !                       matrix polynomials                       *
    !*****************************************************************
    !    evaluates a matrix polynomial                               *
    !    for a given scalar input                                    *
    !                                                                *
    !  @param x - scalar at which to evaluate the MP                 *
    !  @param mp - matrix polynomial, assumes a list of tridiagonal  *
    !       matrices of data type trid, in *increasing* degree       *
    !       e.g. A_0, A_1, A_2                                       *
    !  @return comp - [comp, compD1, compD1], the value of a function*
    !                   plus the values of its derivatives. In the   *
    !                   future, these will likely be passed in       *
    !                   as the running values of an outside program  *
    !                                                                *                                     
    !*****************************************************************
    function revTriHorner(x, mp) 
        implicit none
        !arguement variables
        complex(kind = dp)  :: x
        type(trid_mp)       :: mp
        !local variables
        integer             :: i, xx
        !return variables
        type(trid)          :: comp
        type(trid)          :: compD2
        type(trid)          :: compD1
        type(trid)          :: revTriHorner(3)
            
        !initial coefficients, initialize the derivatives to 0
        comp = mp%coef(mp%degree + 1)
        !get 0 arrays for the first and second derivatives
        !first derivative, three trid data types, all zero, of the right length (size)
        allocate(compD1%upper(mp%size -1), compD1%diag(mp%size), compD1%lower(mp%size -1))
        compD1%upper = cmplx(0,0,kind = dp)
        compD1%diag = cmplx(0,0, kind = dp)
        compD1%lower = cmplx(0,0, kind = dp)
        
        !second derivative
        allocate(compD2%upper(mp%size -1), compD2%diag(mp%size), compD2%lower(mp%size -1))
        compD2%upper = cmplx(0,0, kind = dp)
        compD2%diag = cmplx(0,0, kind = dp)
        compD2%lower = cmplx(0,0, kind = dp)
        
        !horner's method
        !evaluate the reverseal polynomial
        do i = mp%degree, 1, -1
            xx = 1/x
            !second derivative
            compD2%lower = compD2%lower * xx + compD1%lower
            compD2%diag = compD2%diag * xx + compD1%lower
            compD2%upper = compD2%upper * xx + compD1%lower
            !first derivative
            compD1%lower = compD1%lower * xx + comp%lower
            compD1%diag = compD1%diag * xx + comp%lower
            compD1%upper = compD1%upper * xx + comp%lower
            !full value
            comp%lower = comp%lower * xx + mp%coef(i)%lower
            comp%diag = comp%diag * xx + mp%coef(i)%diag
            comp%upper = comp%upper * xx + mp%coef(i)%upper
        end do

        revTriHorner = [comp, compD1, compD2]
        return
    end function
   
    !*****************************************************************
    !                       horner's method for                      *
    !                       matrix polynomials                       *
    !*****************************************************************
    !    evaluates a matrix polynomial                               *
    !    for a given scalar input                                    *
    !                                                                *
    !  @param x - scalar at which to evaluate the MP                 *
    !  @param mp - matrix polynomial, assumes a list of tridiagonal  *
    !       matrices of data type trid, in *increasing* degree       *
    !       e.g. A_0, A_1, A_2                                       *
    !  @return comp - [comp, compD1, compD1], the value of a function*
    !                   plus the values of its derivatives. In the   *
    !                   future, these will likely be passed in       *
    !                   as the running values of an outside program  *
    !                                                                *                             
    !*****************************************************************
    function triHorner(x, mp)
        implicit none
        !arguement variables
        complex(kind = dp)  :: x
        type(trid_mp)       :: mp
        !local variables
        integer             :: i, xx
        !return variables
        type(trid)          :: comp
        type(trid)          :: compD2
        type(trid)          :: compD1
        type(trid)          :: triHorner(3)
            
            
        !first derivative, three trid data types, all zero, of the right length (size)
        allocate(compD1%upper(mp%size -1), compD1%diag(mp%size), compD1%lower(mp%size -1))
        compD1%lower = cmplx(0,0,kind = dp) !sets everything in the array to zero
        compD1%diag = cmplx(0,0,kind = dp) 
        compD1%upper = cmplx(0,0,kind = dp) 
        

        !second derivative
        allocate(compD2%upper(mp%size -1), compD2%diag(mp%size), compD2%lower(mp%size -1))
        compD2%lower = cmplx(0,0,kind = dp)
        compD2%diag = cmplx(0,0,kind = dp)
        compD2%upper = cmplx(0,0,kind = dp)
        
        !horner's method
        comp = mp%coef(mp%degree + 1)
        do i = mp%degree, 1, -1 
            !second derivative
            compD2%lower = compD2%lower * x + compD1%lower
            compD2%diag = compD2%diag * x + compD1%diag
            compD2%upper = compD2%upper * x + compD1%upper
                   
            !first derivative
            compD1%lower = compD1%lower * x + comp%lower
            compD1%diag = compD1%diag * x + comp%diag
            compD1%upper = compD1%upper * x + comp%upper
        
            !full value
            comp%lower = comp%lower * x + mp%coef(i)%lower
            comp%diag = comp%diag * x + mp%coef(i)%diag
            comp%upper = comp%upper * x + mp%coef(i)%upper
        end do
        compD2%lower = compD2%lower * 2.0
        compD2%diag = compD2%diag * 2.0
        compD2%upper = compD2%upper * 2.0
        
        triHorner = [comp, compD1, compD2]
        return
    end function
    
    !*****************************************************************
    !                       Proper Hyman's method                
    !*****************************************************************
    !    Assumes a proper tridiagonal matrix with no zero entries 
    !      along the diagonals, excluding the top and bottom entries 
    !      of the main and upper diagonal, which can be zero
    !                        
    !  @param R(u):  A tridiagonal matrix with entries on the main,  
    !                upper and superupper diagonals                  
    !  @return x: the solution of Rx = y, where R and y come from the
    !              Hyman decomposition                                                           
    !*****************************************************************  
    function PHyman(Rs) result(res)
        implicit none
        !arguement variables
        type(trid)          :: Rs(3) !evaluation in first entry, 1st deriv in second, 2nd deriv in third
        !local variables
        complex(kind=dp), allocatable :: bottom(:), middle(:), top(:), y(:), y1(:), x1(:)
        integer          :: n,i 
        !return variable
        complex(kind=dp)   :: q, q1, res
        
        
        n = size(Rs(1)%diag)
        !get x = R^-1*y
        allocate(y(n-1))
        y = cmplx(0,0,kind = dp)
        y(n-1) = Rs(1)%diag(n)
        y(n-2) = Rs(1)%upper(n-1) !create y Hyman's decomposition
        !triBack writes x into y
        bottom = Rs(1)%lower
        middle = Rs(1)%diag(2:n-1)
        top = Rs(1)%upper(2:n-2)
        call triBack(bottom, middle, top, y) !solve Hyman's decomposition for x
        
        !at this point, y is x
        q = -(Rs(1)%diag(1) * y(1) + Rs(1)%upper(1) * y(2))
        
        
        !get x' = y' - R'(R^-1y)
        allocate(y1(n-1)) 
        y1 = cmplx(0,0,kind = dp)
        y1(n-1) = Rs(2)%diag(n)
        y1(n-2) = Rs(2)%upper(n-1) !create y' Hyman's decomposition
        
!         print *, 'printing y prime'
!         do i =1, size(y1)
!             print *, y1(i)
!         end do
!
!         print *, 'printing x'
!         do i =1, size(y)
!             print *, y(i)
!         end do
!
!         print *, 'printing R'
!         print *, 'printing lower'
!         do i = 1, size(Rs(2)%lower)
!             print *, Rs(2)%lower(i)
!         end do
!         print *, 'printing main'
!         do i = 1, size(Rs(2)%diag(2:n-1))
!             print *, Rs(2)%diag(i+1)
!         end do
!         print *, 'size of upper is', size(Rs(2)%upper(2:n-2))
        
        allocate(x1(size(y)))
        x1 = triMult(Rs(2)%lower, Rs(2)%diag(2:n-1), Rs(2)%upper(2:n-2), y)
!         print *, 'printing output of Rprimex'
!         do i =1, size(x1)
!             print *, x1(i)
!         end do
                    
        x1 = y1 - triMult(Rs(2)%lower, Rs(2)%diag(2:n-1), Rs(2)%upper(2:n-2), y) !Rx' = x1
        call triback(bottom, middle, top, x1)!solve Rx' = x1 for x' 
		
!         print *, 'printing x prime'
!         do i =1, size(x1)
!             print *, x1(i)
!         end do
        
!         print *, 'htprime1 is', Rs(2)%diag(1)
!         print *, 'htprime2 is', Rs(2)%upper(1)
!
!         print *, 'x1 is', y(1)
!         print *, 'x2 is', y(2)
!
!         print *, 'ht1 is', Rs(1)%diag(1)
!         print *, 'ht2 is', Rs(1)%upper(1)
!
!         print *, 'xprime1 is', x1(1)
!         print *, 'xprime2 is', x1(2)
!
!         print *, 'the first term in the dot product is', Rs(2)%diag(1)*y(1) + Rs(2)%upper(1)*y(2)
!         print *, 'the second term in the dot product is', Rs(1)%diag(1)*x1(1) + Rs(1)%upper(1)*x1(2)
        
        !y. is x, x1 is x'
        q1 = -((Rs(2)%diag(1)*y(1) + Rs(2)%upper(1)*y(2)) + (Rs(1)%diag(1)*x1(1) + Rs(1)%upper(1)*x1(2)))
        
        !print *, ' q prime is', q1
        !get r'/r using derivative of logarithm
        res = cmplx(0,0,kind = dp)
        do i = 1, n-1
            res = res + Rs(2)%lower(i)/Rs(1)%lower(i)
        end do
        res = res + q1/q
        return
    end function PHyman
    
    !*****************************************************************
    !                       Back substituion for                     *
    !               "upper" tridiagonal matrix polynomials           *
    !*****************************************************************
    !    Solves Ax = y, where A has three diagonal bands on the main *
    !    diagonal, the upper diagonal and the "super" upper diagonal *
    !    written for the R matrix of an upper hessenburg Matrix      *
    !    Polynomial Hyman's method decomposition. Below, n is the    *
    !    size of the pre-Hyman's method matrix polynomial.            *
    !  @param y:  n - 1 vertical vector, where n is the size of the  *
    !              original matrix polynomial                        *
    !  @param R(u):  A tridiagonal matrix with entries on the main,  *
    !                upper and superupper diagonals                  *
    !  @return y: the solution to the system will be overwritten in  *
    !              the y parameter                                   *                          
    !*****************************************************************
    subroutine triBack(bottom, middle, top, y)
        implicit none
        !arguement variables
        complex(kind=dp) :: bottom(:), middle(:), top(:), y(:)
        !local variables
        integer             :: n, i
        n = size(bottom)
        y(n) = y(n)/bottom(n) !first case 
        y(n-1) = (y(n-1) - middle(n-1) * y(n))/bottom(n-1) !second case
        do i = n-2,1,-1
            y(i) = (y(i) - top(i)*y(i+2) - middle(i)*y(i+1))/bottom(i)
        end do 
    end subroutine triBack   
    
end module tridiag