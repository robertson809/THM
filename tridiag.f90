module tridiag
    implicit none
    integer, parameter  :: dp = kind(1.d0)
	
    !****************************************
    !				type trid				*
    !****************************************
    ! This is a data type that we've defined, with the tre
    type trid
        complex(kind=dp), allocatable :: diag(:), upper(:), lower(:)
    end type trid
    
    !****************************************
    !				type trid_mp		    *
    !       a tridiagonal matrix polynomial	*
    !****************************************
    ! This is a tridiagonal matrix polynoial
    type trid_mp
        type(trid), allocatable :: coef(:)
        integer                 :: size
        integer                 :: degree
    end type trid_mp
    
contains
    !****************************************
    !				print trid				*
    !****************************************
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
    
    !****************************************
    !				print trid_mp			*
    !****************************************
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
            write(*,'(A, I1)') 'x^', i
            if (i /= degree) then
                write(*, '(A)') '+'
            end if
        end do
    end subroutine print_trid_mp
    

    
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
        real, dimension(2) :: parts
        !return variables
        type(trid)          :: comp
        type(trid)          :: compD2
        type(trid)          :: compD1
        type(trid)          :: triHorner(3)
            
            
        
        !initial coefficients, initialize the derivatives to 0
        comp = mp%coef(mp%degree + 1)
        !get 0 arrays for the first and second derivatives
        
        !first derivative, three trid data types, all zero, of the right length (size)
        allocate(compD1%upper(mp%size -1), compD1%diag(mp%size), compD1%lower(mp%size -1))
        compD1%upper = cmplx(0,0,kind = dp)
        compD1%diag = mp%coef(1)%diag * 0
        compD1%lower = mp%coef(1)%upper * 0
        
        !second derivative
        allocate(compD2%upper(mp%size -1), compD2%diag(mp%size), compD2%lower(mp%size -1))
        compD2%upper = mp%coef(1)%lower * 0
        compD2%diag = mp%coef(1)%diag * 0
        compD2%lower = mp%coef(1)%upper * 0
        
        !horner's method
        
        !evaluate the reverseal polynomial if |x| > 1
        if (dot_product(transfer(x, parts), transfer(x, parts)) > 1) then 
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
        else 
            do i = mp%degree, 1, -1 
            
            !second derivative
            compD2%lower = compD2%lower * x + compD1%lower
            compD2%diag = compD2%diag * x + compD1%lower
            compD2%upper = compD2%upper * x + compD1%lower
            
            !first derivative
            compD1%lower = compD1%lower * x + comp%lower
            compD1%diag = compD1%diag * x + comp%lower
            compD1%upper = compD1%upper * x + comp%lower
            
            !full value
            comp%lower = comp%lower * x + mp%coef(i)%lower
            comp%diag = comp%diag * x + mp%coef(i)%diag
            comp%upper = comp%upper * x + mp%coef(i)%upper
            end do
        end if
        !write(*, '(A)') 'Line 157'
        triHorner = [comp, compD1, compD2]
        return
    end function
    
end module tridiag