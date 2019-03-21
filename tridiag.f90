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
    
    !****************************************************************
    !                   Performs the horners step                   *
    !****************************************************************
    !                                                               *
    ! @param coef - the current coefficient in the matrix polynomial*
    ! @param x - the scalar at which the polynomial is being        *
    ! evaluated                                                     *
    ! @param comp - the running computation in horner's method      *
    !****************************************************************
    subroutine horner_step(coef, x, comp, size)
        implicit none
        ! argument variables
        type(trid)               :: coef
        complex(kind=dp)         :: x
        type(trid)               :: comp
        integer                  :: size
        ! local variables
        integer                  :: i
        
        !Horner step
        do i = 1, size - 1
                comp%lower(i) = (comp%lower(i) * x) + coef%lower(i)
                comp%diag(i) = (comp%diag(i) * x) + coef%lower(i)
                comp%upper(i) = (comp%upper(i) * x) + coef%lower(i)
        end do
        comp%diag(size) = (comp%diag(size) * x) + coef%diag(size)
        
    end subroutine
    
    !*****************************************************************
    !                       horner's method for                      *
    !                       matrix polynomials                       *
    !*****************************************************************
    !    evaluates a matrix polynomial                               *
    !    for a given scalar input                                    *
    !                                                                *
    !  @param x - scalar at which to evaluate the MP                 *
    !  @param a - matrix polynomial, assumes a list of tridiagonal   *
    !       matrices of data type trid, in *increasing* degree       *
    !  @return comp - the matrix output                              *                                                 
    !*****************************************************************
    function horner(x, mp) result(comp)
        implicit none
        !arguement variables
        complex(kind = dp)  :: x
        type(trid_mp)       :: mp
        !local variables
        integer             :: i
        !return variables
        type(trid)          :: comp

        !horner's method
        comp = mp%coef(mp%degree + 1)
        do i = mp%degree, 1, -1 
            !this function updates comp
            call horner_step(mp%coef(i), x, comp, mp%size)
        end do
        
    end function
    
end module tridiag