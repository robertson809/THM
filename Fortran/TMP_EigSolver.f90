module TMP_EigSolver
    use fpml
    implicit none
    !****************************************************************
    !				           Type tri                             *
    !****************************************************************
    !   Tridiagonal matrix type contains complex allocatable arrays
    !   for storing lower, upper, and main diagonal
    !****************************************************************
    type trid
        complex(kind=dp), allocatable :: lower(:), upper(:), main(:)
    end type trid
    !****************************************************************
    !				           Type trid_mp		                    *
    !****************************************************************
    !   Tridagonal matrix polynomial type contains allocatable array
    !   of tridiagonal matrices, as well as the size and degree of
    !   the matrix polynomial.
    !****************************************************************
    type trid_mp
        type(trid), allocatable :: coeff(:)
        integer                 :: size
        integer                 :: degree
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
        integer                         :: i, j, num_eig
        real(kind=dp)                   :: r
        real(kind=dp), allocatable      :: alpha(:)
        complex(kind=dp)                :: lag_corr1, lag_corr2, z
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! store norm of matrix coefficients
        allocate(alpha(mat_poly%degree+1))
        do i=1,mat_poly%degree+1
            call FroNorm(mat_poly%size,mat_poly%coeff(i),alpha(i))
        end do
        ! initial estimates
        
        ! main loop
        
        ! deallocate alpha memory
        deallocate(alpha)
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
    end subroutine InitEst
    !****************************************************************
    !				           FroNorm                              *
    !****************************************************************
    !   Computes the Frobenius Norm of a tridiagonal matrix.
    !****************************************************************
    subroutine FroNorm(n,mat,alpha)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        type(trid), intent(in)          :: mat
        real(kind=dp), intent(out)      :: alpha
        ! local variables
        integer                         :: j
        real(kind=dp)                   :: scale, sum
        complex(kind=dp)                :: tmp(3)
        ! external subroutines
        external                        :: zlassq
        ! intrinsic functions
        intrinsic                       :: sqrt
        
        ! initialize scale and sum
        scale = 0.0_dp
        sum = 1.0_dp
        ! column 1
        tmp(1) = mat%main(1)
        tmp(2) = mat%lower(1)
        call zlassq(2, tmp, 1, scale, sum)
        ! columns 2,...,n-1
        do j=2,n-1
            tmp(1) = mat%upper(j-1)
            tmp(2) = mat%main(j)
            tmp(3) = mat%lower(j)
            call zlassq(3, tmp, 1, scale, sum)
        end do
        ! column n
        tmp(1) = mat%upper(n-1)
        tmp(2) = mat%main(n)
        call zlassq(2, tmp, 1, scale, sum)
        ! store result
        alpha = scale*sqrt(sum)
    end subroutine
end module TMP_EigSolver