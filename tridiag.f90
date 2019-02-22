module tridiag
    implicit none
    integer, parameter  :: dp = kind(1.d0)
	
    !****************************************
    !				type trid				*
    !****************************************
    !This is a data type that we've defined
    type trid
        complex(kind=dp), allocatable :: diag(:), upper(:), lower(:)
    end type trid
contains
	
    !****************************************
    !				print trid				*
    !****************************************
    subroutine print_trid(a,n)
        implicit none
        ! argument variables
        integer         :: n
        type(trid)      :: a
        ! local variables
        integer         :: i
			
        ! print statements
        write(*,'(F0.0,SP,F0.0,"i", F0.0,SP,F0.0,"i")') a%diag(1), a%upper(1)
        do i=2,n-1
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
        
        do i = 1, size(a%diag) - 1
            comp%lower(i) = a%lower(i) * c
            comp%upper(i) = a%upper(i) * c
            comp%diag(i) = a%diag(i) * c
        end do
        comp%diag(i + 1) = a%diag(i + 1) * c
        return
    end function
end module tridiag