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
        call InitEst(mat_poly,eigval,eigvec)
        ! print out roots (testing purposes)
        do j=1,mat_poly%size*mat_poly%degree
            write(*,*) eigval(j)
        end do
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
        integer                         :: j, k
        type(trid)                      :: mat
        complex(kind=dp)                :: qmat(mat_poly%size-1,2), x(mat_poly%size)
        ! fpml variables
        integer, parameter              :: itmax = 30
        integer                         :: conv(mat_poly%degree)
        real(kind=dp)                   :: berr(mat_poly%degree), cond(mat_poly%degree)
        complex(kind=dp)                :: poly(mat_poly%degree+1)
        ! external function
        external                        :: zdotc
        complex(kind=dp)                :: zdotc
        
        ! set mat to leading coefficient
        mat = mat_poly%coeff(mat_poly%degree+1)
        ! perform qr decomposition of mat
        call TriQR(mat_poly%size,mat,qmat)
        ! build vectors that span column space of mat
        do j=1,mat_poly%size
            ! x is set to jth unit vector
            x = (/ (cmplx(0,0,kind=dp), k=1,mat_poly%size) /)
            x(j) = 1.0_dp
            ! x is multiplied by qmat
            call QMult(mat_poly%size,qmat,x)
            ! build coefficients of scalar polynomial
            do k=1,mat_poly%degree+1
                poly(k) = zdotc(mat_poly%size,x,1,TriMult(mat_poly%size,mat_poly%coeff(k),x),1)
            end do
            ! compute roots of polynomial
            call main(poly, mat_poly%degree, eigval((j-1)*mat_poly%degree+1:j*mat_poly%degree), berr, cond, conv, itmax)
        end do
    end subroutine InitEst
    !****************************************************************
    !				           TriMult                              *
    !****************************************************************
    !   Multiply a vector by a tridiagonal matrix.
    !****************************************************************
    function TriMult(n,mat,x) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        type(trid), intent(in)          :: mat
        complex(kind=dp), intent(in)    :: x(n)
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: res(n)
        
        ! row 1
        res(1) = mat%main(1)*x(1) + mat%upper(1)*x(2)
        ! row 2,...,n-1
        do i=2,n-1
            res(i) = mat%lower(i-1)*x(i-1) + mat%main(i)*x(i) + mat%upper(i)*x(i+1)
        end do
        ! row n
        res(n) = mat%lower(n-1)*x(n-1) + mat%main(n)*x(n)
        return
    end function TriMult
    !****************************************************************
    !				           QMult                                *
    !****************************************************************
    !   Multiply a vector by a qmat that stores the rotators used in
    !   a QR decomposition.
    !****************************************************************
    subroutine QMult(n,qmat,x)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(in)    :: qmat(n-1,2)
        complex(kind=dp), intent(out)   :: x(n)
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: tmp(2)
        ! intrinsic function
        intrinsic                       :: conjg
        
        ! rotations 1,...,n-1
        do i=1,n-1
            ! update x
            tmp(1) = conjg(qmat(i,1))*x(i) + conjg(qmat(i,2))*x(i+1)
            tmp(2) = -qmat(i,2)*x(i) + qmat(i,1)*x(i+1)
            x(i) = tmp(1)
            x(i+1) = tmp(2)
        end do
    end subroutine QMult
    !****************************************************************
    !				           TriQR                                *
    !****************************************************************
    !   Forms the QR decomposition of a tridiagonal matrix. On input
    !   mat is a tridiagonal matrix stored using type(trid). On output
    !   the three diagonals of R are stored in mat and the rotators
    !   used are stored in qmat.
    !****************************************************************
    subroutine TriQR(n,mat,qmat)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        type(trid), intent(inout)       :: mat
        complex(kind=dp), intent(out)   :: qmat(n-1,2)
        ! local variables
        integer                         :: i
        real(kind=dp)                   :: scale, sum
        complex(kind=dp)                :: tmp(2)
        ! external subroutine
        external                        :: zlassq
        ! intrinsic function
        intrinsic                       :: conjg, sqrt
        
        ! rotations 1,...,n-2
        do i=1,n-2
            ! compute r = sqrt(main(i)^2+lower(i)^2)
            tmp(1) = mat%main(i)
            tmp(2) = mat%lower(i)
            scale = 0.0_dp
            sum = 1.0_dp
            call zlassq(2, tmp, 1, scale, sum)
            sum = scale*sqrt(sum)
            ! store rotator
            qmat(i,1) = mat%main(i)/sum
            qmat(i,2) = mat%lower(i)/sum
            ! update mat
            mat%main(i) = sum
            tmp(1) = conjg(qmat(i,1))*mat%upper(i) + conjg(qmat(i,2))*mat%main(i+1)
            tmp(2) = -qmat(i,2)*mat%upper(i) + qmat(i,1)*mat%main(i+1)
            mat%upper(i) = tmp(1)
            mat%main(i+1) = tmp(2)
            tmp(1) = conjg(qmat(i,2))*mat%upper(i+1)
            tmp(2) = qmat(i,1)*mat%upper(i+1)
            mat%lower(i) = tmp(1)
            mat%upper(i+1) = tmp(2)
        end do
        ! rotation n-1
        i = n-1
        ! compute r = sqrt(main(i)^2+lower(i)^2)
        tmp(1) = mat%main(i)
        tmp(2) = mat%lower(i)
        scale = 0.0_dp
        sum = 1.0_dp
        call zlassq(2, tmp, 1, scale, sum)
        sum = scale*sqrt(sum)
        ! store rotator
        qmat(i,1) = mat%main(i)/sum
        qmat(i,2) = mat%lower(i)/sum
        ! update mat
        mat%main(i) = sum
        tmp(1) = conjg(qmat(i,1))*mat%upper(i) + conjg(qmat(i,2))*mat%main(i+1)
        tmp(2) = -qmat(i,2)*mat%upper(i) + qmat(i,1)*mat%main(i+1)
        mat%upper(i) = tmp(1)
        mat%main(i+1) = tmp(2)
        mat%lower(i) = 0.0_dp
    end subroutine
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