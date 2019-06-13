module TMP_EigSolver2
    use fpml
    implicit none
    !****************************************************************
    !				           Type tri                             *
    !****************************************************************
    !   Tridiagonal matrix type contains complex allocatable arrays
    !   for storing the diagonals dl, d, and du (lower, main and uper).
    !****************************************************************
    type trid
        complex(kind=dp), allocatable   :: dl(:), d(:), du(:)
    end type trid
    !****************************************************************
    !				           Type trid_mp		                    *
    !****************************************************************
    !   Tridagonal matrix polynomial type contains allocatable array
    !   of tridiagonal matrices (coefficients), as well as the size 
    !   and degree of the matrix polynomial.
    !****************************************************************
    type trid_mp
        type(trid), allocatable         :: coeff(:)
        integer                         :: size
        integer                         :: degree
    end type trid_mp
contains
    !****************************************************************
    !				           EigSolver                            *
    !****************************************************************
    !   Main subroutine for computing the eigenvalues and eigenvectors
    !   of a tridiagonal matrix polynomial. 
    !****************************************************************
    subroutine EigSolver(mat_poly,eigval,eigvec,berr,cond,conv,itmax)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        complex(kind=dp), intent(out)   :: eigval(:), eigvec(:,:)
        real(kind=dp), intent(out)      :: berr(:), cond(:)
        integer, intent(out)            :: conv(:)
        integer, intent(in)             :: itmax
        ! local variables
        integer                         :: i, j, num_eig, total_eig
        real(kind=dp)                   :: r
        real(kind=dp), allocatable      :: alpha(:)
        complex(kind=dp)                :: lag_term1, lag_term2, z
        ! intrinsic functions
        intrinsic                       :: abs
        
        call InitEst(mat_poly,eigval,eigvec)
    end subroutine EigSolver
    !****************************************************************
    !				           InitEst                              *
    !****************************************************************
    !   Computes the initial eigenvalue and eigenvector estimates
    !   for the tridiagonal matrix polynomial.
    !****************************************************************
    subroutine InitEst(mat_poly,eigval,eigvec)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        complex(kind=dp), intent(out)   :: eigval(:), eigvec(:,:)
        ! local variables
        integer                         :: clock, j, k
        type(trid)                      :: mat
        complex(kind=dp)                :: x(mat_poly%size), y(mat_poly%size)
        ! lapack variables
        integer                         :: iseed(4), ipiv(mat_poly%size), info
        complex(kind=dp)                :: du2(mat_poly%size)
        ! fpml variables
        integer, parameter              :: itmax = 30
        integer                         :: conv(mat_poly%degree)
        real(kind=dp)                   :: berr(mat_poly%degree), cond(mat_poly%degree)
        complex(kind=dp)                :: poly(mat_poly%degree+1)
        ! intrinsic subroutines
        intrinsic                       :: system_clock
        ! external subroutines
        external                        :: zgttrf, zgttrs, zlarnv
        ! external functions
        real(kind=dp)                   :: dznrm2
        complex(kind=dp)                :: zdotc
        external                        :: zdotc, dznrm2
        
        ! initial iseed
        call system_clock(count = clock)
        iseed = clock + 37*(/ (j-1, j=1,4) /)
        iseed = (/ (mod(iseed(j),4095), j=1,4)/)
        if(mod(iseed(4),2)==0) iseed(4)=iseed(4)+1
        ! set mat to leading coefficient and perform lu decomposition
        mat = mat_poly%coeff(mat_poly%degree+1)
        call zgttrf(mat_poly%size,mat%dl,mat%d,mat%du,du2,ipiv,info)
        ! initial eigenvalue estimates
        do j=1,mat_poly%size
            ! create random vector of unit length
            call zlarnv(5,iseed,mat_poly%size,x)
            x = x/dznrm2(mat_poly%size,x,1)
            ! build coefficients of scalar polynomial
            do k=1,mat_poly%degree+1
                ! multiply random vector by kth coefficient
                y = TriMult(mat_poly%size,mat_poly%coeff(k)%dl,mat_poly%coeff(k)%d,mat_poly%coeff(k)%du,x)
                ! solve linear system Am*x=y
                call zgttrs('N',mat_poly%size,1,mat%dl,mat%d,mat%du,du2,ipiv,y,mat_poly%size,info)
                ! compute dot product x^H*y
                poly(k) = zdotc(mat_poly%size,x,1,y,1)
            end do
            ! compute roots of polynomial
            call fpml_main(poly,mat_poly%degree,eigval((j-1)*mat_poly%degree+1:j*mat_poly%degree),berr,cond,conv,itmax)
        end do
        ! initial eigenvector estimates
        
    end subroutine InitEst
    !****************************************************************
    !				           TriMult                              *
    !****************************************************************
    !   Multiply a vector by a tridiagonal matrix.
    !****************************************************************
    function TriMult(n,dl,d,du,x) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(in)    :: dl(:), d(:), du(:), x(:)
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: res(n)
        
        ! row 1
        res(1) = d(1)*x(1) + du(1)*x(2)
        ! row 2,...,n-1
        do i=2,n-1
            res(i) = dl(i-1)*x(i-1) + d(i)*x(i) + du(i)*x(i+1)
        end do
        ! row n
        res(n) = dl(n-1)*x(n-1) + d(n)*x(n)
        return
    end function TriMult
end module TMP_EigSolver2