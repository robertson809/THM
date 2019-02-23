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
    
    !*************************************
    !       scalar multiplication        *
    !*************************************
    !    returns the value of a tridiagonal
    !           times a scalar           *          
    !*************************************
    function s_mult(a, c) result(comp)
        implicit none 
        ! argument variables
        integer         :: c
        type(trid)      :: a
        !local variables
        integer         :: i
        !return variables 
        type(trid)      :: comp
            
        allocate(comp%diag(size(a%diag)), comp%upper(size(a%upper)), comp%lower(size(a%lower)))
        
        do i = 1, size(a%upper)
            comp%lower(i) = a%lower(i) * c
            comp%upper(i) = a%upper(i) * c
            comp%diag(i) = a%diag(i) * c
        end do
        comp%diag(size(a%diag)) = a%diag(size(a%diag)) * c
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
    
end module tridiag