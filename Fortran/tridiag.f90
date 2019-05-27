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
    integer, parameter          :: dp = kind(1.d0)
    real(kind=dp), parameter    :: eps = 2.0_dp**(-52), mu = 2.0_dp**(-53)
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
    function revTriHorner(mp,x) 
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
    function triHorner(mp, x)
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
    
    
    !*****************************************************************
    !                 Tridiagonal Matrix Multiplication                                  
    !*****************************************************************
    !    Multiplies a upper tridiagonal matrix A by a vector x
    !                                                                
    !  @param x - vector                                              
    !  @param bottom - the diagonal
	!  @param middle - the 1 diagonal
	!  @param top - the 2 diagonal
    !  @return comp - Ax                                 
    !*****************************************************************
    function triMult(bottom, middle, top, x) result(comp)
        implicit none
        !arguement variables
        complex(kind = dp)  :: x(:), bottom(:), middle(:), top(:)
        !local variables
        integer             :: i, n
        !return variables
        complex(kind = dp), allocatable :: comp(:)
        
        n = size(bottom)
        allocate(comp(n))
        if(size(top) > 0) then
            !upper triangular tridiagonal matrix multiplication
            do i = 1, n-2
                comp(i) = bottom(i) * x(i) + middle(i) * x(i+1) + top(i) * x(i + 2)
            end do
            comp(n-1) = bottom(n-1) * x(n-1) + middle(n-1)*x(n)
            comp(n) = bottom(n) * x(n)
        end if
        return
    end function
    
    
      !*****************************************************************
      !                       Proper Hyman's method                
      !*****************************************************************
      !    Assumes a proper tridiagonal matrix with no zero entries 
      !      along the main diagonal
      !                        
      !  @param Rs(u):  Value of the tridiagonal Matrix Polynomial at 
      !                 at u in the first entry, value of the second 
      !                 derivative in the second entry, and value of
      !                 the third derivative in the third
      !  @return res: Newtown correction term, the derivative of the 
      !               of the determinant over the determinant
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
        
          allocate(y1(n-1)) 
          y1 = cmplx(0,0,kind = dp)
          y1(n-1) = Rs(2)%diag(n)
          y1(n-2) = Rs(2)%upper(n-1) !create y' Hyman's decomposition

          call triBack(Rs(1)%lower, Rs(1)%diag(2:n-1), Rs(1)%upper(2:n-2), y) !put x into y
        
          q = -(Rs(1)%diag(1) * y(1) + Rs(1)%upper(1) * y(2)) !at this point, y is x
                    
          y1 = y1 - triMult(Rs(2)%lower, Rs(2)%diag(2:n-1), Rs(2)%upper(2:n-2), y) !y1 = y'-R'x
          call triback(Rs(1)%lower, Rs(1)%diag(2:n-1), Rs(1)%upper(2:n-2), y1)!put x' from Rx' = y1 into y1 
          q1 = -((Rs(2)%diag(1)*y(1) + Rs(2)%upper(1)*y(2)) + (Rs(1)%diag(1)*y1(1) + Rs(1)%upper(1)*y1(2)))
        
          !get r'/r using derivative of logarithm
          res = cmplx(0,0,kind = dp)
          do i = 1, n-1
              res = res + Rs(2)%lower(i)/Rs(1)%lower(i)
          end do
          res = res + q1/q
          return
      end function PHyman 
      
    
      !*****************************************************************
      !                       Special Hyman's method 3             
      !*****************************************************************
      !    Handles the special three dimensional case of a 
      !    tridiagonal matrix polynomial
      !                        
      !  @param Rs(u):  Value of the tridiagonal Matrix Polynomial at 
      !                 at u in the first entry, value of the second 
      !                 derivative in the second entry, and value of
      !                 the third derivative in the third
      !  @return res: Newtown correction term, the derivative of the 
      !               of the determinant over the determinant
      !***************************************************************** 
      function SHyman3(Rs) result(res)
          implicit none
          !arguement variables
          type(trid)          :: Rs(3) !evaluation in first entry, 1st deriv in second, 2nd deriv in third
          !return variable
          complex(kind=dp)   :: res
                             
          res = (Rs(2)%diag(1)*(Rs(1)%diag(2)*Rs(1)%diag(3) - Rs(1)%lower(2)*Rs(1)%upper(2)) &
          + Rs(1)%diag(1)*(Rs(2)%diag(2)*Rs(1)%diag(3) + Rs(1)%diag(2)*Rs(2)%diag(3) - Rs(2)%lower(2)*Rs(1)%upper(2) &
          -Rs(1)%lower(2)*Rs(2)%upper(2)) - Rs(2)%upper(1)*(Rs(1)%lower(1)*Rs(1)%diag(3)) &
          - Rs(1)%upper(1)*(Rs(2)%lower(1)*Rs(1)%diag(3) + Rs(1)%lower(1)*Rs(2)%diag(3)))/ &
          (Rs(1)%diag(1)*(Rs(1)%diag(2)*Rs(1)%diag(3)-Rs(1)%lower(2)*Rs(1)%upper(2)) - &
          Rs(1)%upper(1)*Rs(1)%lower(1)*Rs(1)%diag(3)) ! good luck reading this
         
      end function SHyman3 
      
    
      !*****************************************************************
      !                       Special Hyman's method 2                
      !*****************************************************************
      !    Handles the special two dimensional case of a 
      !    triangular matrix
      !                        
      !  @param Rs(u):  Value of the tridiagonal Matrix Polynomial at 
      !                 at u in the first entry, value of the second 
      !                 derivative in the second entry, and value of
      !                 the third derivative in the third
      !  @return res: Newtown correction term, the derivative of the 
      !               of the determinant over the determinant
      !***************************************************************** 
      function SHyman2(Rs) result(res)
          implicit none
          !arguement variables
          type(trid)          :: Rs(3) !evaluation in first entry, 1st deriv in second, 2nd deriv in third
          !return variable
          complex(kind=dp)   :: res
                
          res = (Rs(2)%diag(1)*Rs(1)%diag(2) + Rs(2)%diag(2)*Rs(1)%diag(1) - &
          Rs(2)%upper(1)*Rs(1)%lower(1) - Rs(1)%upper(1)*Rs(2)%lower(1)) / &
          ((Rs(1)%diag(1)*Rs(1)%diag(2)-Rs(1)%upper(1)*Rs(1)%lower(1)))

          return
      end function SHyman2
      
    
    !*****************************************************************
    !                       Inflate            
    !*****************************************************************
    !    Handles the potential underflow problem that could arise 
    !    during backward substitution from small values in the 
    !    subdiagonal. Inflates all values on the subdiagaonl
    !    that are less than mu to mu.
    !
    !  @param A: tridiagonal matrix
    !  
    !  @return A: tridiagonal matrix with small subdiagonal 
    !              entries inflated to mu    
    !***************************************************************** 
    subroutine inflate(A)
        implicit none
        !arguement variables
        type(trid)       :: A
        !local variables
        integer             :: i
        do i = 1, size(A%lower)
            if (abs(A%lower(i)) < mu) then
                A%lower(i) = mu
            end if
        end do 
    end subroutine inflate 
    
    
    !*****************************************************************
    !                       Hymans           
    !*****************************************************************
    !    Helper function to HStep, handles the special two and three
    !   dimensional cases
    !                        
    !  @param A: tridiagonal matrix
    !  @return step: the step value  
    !***************************************************************** 
    function Hymans(Rs) result(step)
        implicit none
        !arguement variables
        type(trid)         :: Rs(3)
        !local varaibles
        integer             :: n
        !return variables
        complex(kind = dp)   :: step
        n = size(Rs(1)%diag)
        if (n > 3) then 
            step = PHyman(Rs)
        elseif (n == 3) then
            step = SHyman3(Rs)
        elseif (n == 2) then
            step = SHyman2(Rs)
        else 
            print *, 'Thats a scalar polynomial. You need to go somewhere else'
            call Exit(0)
        end if
        return 
    end function Hymans
        
    
    !*****************************************************************
    !                       HStep            
    !*****************************************************************
    !    Takes a matrix polynomial A and a complex scalar x and 
    !   returns the Newton step using Hyman's method for Matrix 
    !   polynomials
    !                        
    !  @param A: Matrix polynomial
    !  @param x: complex scalar
    !  @return step: Newtown correction term, the derivative of the 
    !               of the determinant over the determinant
    !***************************************************************** 
    function Hstep(A, x) result(step)
        implicit none
        !arguement variables
        complex(kind=dp)    :: x
        type(trid_mp)       :: A
        ! local variables 
        integer     :: n, jold, jnew, j
        type(trid)      :: eval(3)
        !return variable
        complex(kind = dp)   :: step
        
        eval = triHorner(A,x)
        call inflate(eval(1))
        step = Hymans(eval)
        return
    end function Hstep
end module tridiag