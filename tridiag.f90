module tridiag
    implicit none
    integer, parameter  :: dp = kind(1.d0)
	
    !****************************************
    !				type trid				*
    !****************************************
    !This is a data type that we've defined, with the tre
    type trid
        complex(kind=dp), allocatable :: diag(:), upper(:), lower(:)
    end type trid
    
    !****************************************
    !				type trid_mp		    *
    !       a tridiagonal matrix polynomial	*
    !****************************************
    !This is a tridiagonal matrix polynoial
    type trid_mp
        type(trid), allocatable :: coef(:)
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
        integer             :: n
        n = size(a%coef)
        
        !printing
        do i = 1, n
            call print_trid(a%coef(i))
            write(*, '(A)') 'x'
            write(*, '(A)') '+'
        end do
    end subroutine print_trid_mp
    
    !*************************************
    !       scalar multiplication        *
    !*************************************
    !    returns the value of a tridiagonal
    !           times a scalar           *          
    !*************************************
    function s_mult(tri, c) result(comp)
        implicit none 
        ! argument variables
        complex(kind = dp)         :: c
        type(trid)                 :: tri
        !local variables
        integer                    :: i
        !return variables 
        type(trid)                 :: comp
            
        allocate(comp%diag(size(tri%diag)), comp%upper(size(tri%upper)), comp%lower(size(tri%lower)))
        
        do i = 1, size(tri%upper)
            comp%lower(i) = tri%lower(i) * c
            comp%upper(i) = tri%upper(i) * c
            comp%diag(i) = tri%diag(i) * c
        end do
        comp%diag(size(tri%diag)) = tri%diag(size(tri%diag)) * c
        return
    end function
    
    !*************************************
    !           matrix addition          *
    !*************************************
    !    returns the value of a tridiagonal
    !           times a scalar           *          
    !*************************************
    function m_add(a, b) result(comp)
        implicit none
        !argument variables
        type(trid)      :: a, b
        !local variables
        integer        :: i
        !return variables
        type(trid)      :: comp
            
        allocate(comp%diag(size(a%diag)), comp%upper(size(a%upper)), comp%lower(size(a%lower)))
        
        do i = 1, size(a%lower)
            comp%lower(i) = a%lower(i) + b%lower(i)
            comp%upper(i) = a%upper(i) + b%upper(i)
            comp%diag(i) = a%diag(i) + b%diag(i)
        end do
        comp%diag(size(a%diag)) = a%diag(size(a%diag)) + b%diag(size(a%diag))
        
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
    !  @param a - matrix polynomial, assumes a list of tridiagonal   *
    !       matrices of data type trid, in *decreasing* degree       *
    !  @return comp - the matrix output                              *                                                 
    !*****************************************************************
    function horner(x, mp) result(comp)
        implicit none
        !arguement variables
        complex(kind = dp)  :: x
        type(trid_mp)       :: mp
        !local variables
        integer             :: i
        integer             :: degree
        !return variablescoef
        type(trid)          :: comp
            
        degree = size(mp%coef)

        !horner's method
        comp = mp%coef(1)
        do i = 2, degree
            comp = m_add(mp%coef(i), s_mult(comp, x))
        end do
        
    end function
    
end module tridiag