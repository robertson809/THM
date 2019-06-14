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
        
        ! store norm of matrix coefficients
        allocate(alpha(mat_poly%degree+1))
        write(*,*) 'forbenius norm of coefficients'
        do i=1,mat_poly%degree+1
            call FroNorm(mat_poly%size,mat_poly%coeff(i),alpha(i))
            alpha(i) = alpha(i)*(3.8_dp*(i-1)+1)
        end do
        ! initial estimates
        write(*,*) 'initial estimates'
        call InitEst(mat_poly,eigval,eigvec)
    end subroutine EigSolver
    !****************************************************************
    !				           Hyman                                *
    !****************************************************************
    !   Applies the original polynomial and Hyman's method to compute
    !   the Laguerre correction terms. In the process, the backward
    !   error is computed and stopping criterion is checked. If the
    !   stopping criterion is true, then conv is set to 1 and the
    !   eigenvalue approximation will not be updated.
    !****************************************************************
    subroutine Hyman(mat_poly,alpha,z,r,v,lag_term1,lag_term2,berr,cond,conv)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        integer, intent(out)            :: conv
        real(kind=dp), intent(in)       :: alpha(mat_poly%degree+1), r
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: z, v(mat_poly%size)
        complex(kind=dp), intent(out)   :: lag_term1, lag_term2
        ! local variables
        integer                         :: k
        type(trid)                      :: triHorner(3)
        real(kind=dp)                   :: res
        complex(kind=dp)                :: y(mat_poly%size-1,3), q(3), sum
        ! external function
        real(kind=dp)                   :: dznrm2
        external                        :: dznrm2
        
        ! polynomial evaluation (note triHorner(1) gets allocated memory seperately from mat_poly%coeff)
        triHorner(1) = mat_poly%coeff(mat_poly%degree+1)
        ! polynomial derivative evaluation
        allocate(triHorner(2)%du(mat_poly%size -1), triHorner(2)%d(mat_poly%size), triHorner(2)%dl(mat_poly%size -1))
        triHorner(2)%du = cmplx(0,0,kind=dp)
        triHorner(2)%d = cmplx(0,0,kind=dp)
        triHorner(2)%dl = cmplx(0,0,kind=dp)
        ! polynomial 2nd derivative evaluation
        allocate(triHorner(3)%du(mat_poly%size -1), triHorner(3)%d(mat_poly%size), triHorner(3)%dl(mat_poly%size -1))
        triHorner(3)%du = cmplx(0,0,kind=dp)
        triHorner(3)%d = cmplx(0,0,kind=dp)
        triHorner(3)%dl = cmplx(0,0,kind=dp)
        ! Horner's method
        do k=mat_poly%degree,1,-1
            ! update polynomial 2nd derivative evaluation
            triHorner(3)%du = triHorner(3)%du*z + triHorner(2)%du
            triHorner(3)%d = triHorner(3)%d*z + triHorner(2)%d
            triHorner(3)%dl = triHorner(3)%dl*z + triHorner(2)%dl    
            ! update polynomial derivative evaluation
            triHorner(2)%du = triHorner(2)%du*z + triHorner(1)%du
            triHorner(2)%d = triHorner(2)%d*z + triHorner(1)%d
            triHorner(2)%dl = triHorner(2)%dl*z + triHorner(1)%dl
            ! update polynomial evaluation
            triHorner(1)%du = triHorner(1)%du*z + mat_poly%coeff(k)%du
            triHorner(1)%d = triHorner(1)%d*z + mat_poly%coeff(k)%d
            triHorner(1)%dl = triHorner(1)%dl*z + mat_poly%coeff(k)%dl
        end do
        ! backward error terms
        berr = alpha(mat_poly%degree+1)
        do k=mat_poly%degree,1,-1
            berr = r*berr + alpha(k)
        end do
        res = dznrm2(mat_poly%size,TriMult(mat_poly%size,triHorner(1)%dl,triHorner(1)%d,triHorner(1)%du,v),1)
        ! compute lag terms or backward error and condition number
        if(res > eps*berr) then
            ! check for subdiagonal zeros in triHorner(1)
            do k=1,mat_poly%size-1
                if(abs(triHorner(1)%dl(k)) .le. small) triHorner(1)%dl(k) = small
            end do
            ! store initial y values
            y = cmplx(0,0,kind=dp)
            y(mat_poly%size-2,1) = triHorner(1)%du(mat_poly%size-1)
            y(mat_poly%size-1,1) = triHorner(1)%d(mat_poly%size)
            y(mat_poly%size-2,2) = triHorner(2)%du(mat_poly%size-1)
            y(mat_poly%size-1,2) = triHorner(2)%d(mat_poly%size)
            y(mat_poly%size-2,3) = triHorner(3)%du(mat_poly%size-1)
            y(mat_poly%size-1,3) = triHorner(3)%d(mat_poly%size)
            ! solve linear systems: Rx=y
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,1))
            ! solve linear system Rx'=y'-R'x
            y(:,2) = y(:,2) - TriMult(mat_poly%size-1,triHorner(2)%dl,triHorner(2)%d,triHorner(2)%du,y(:,1))
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,2))
            ! solve linear system Rx''=y''-R''x-2R'x'
            y(:,3) = y(:,3) - TriMult(mat_poly%size-1,triHorner(3)%dl,triHorner(3)%d,triHorner(3)%du,y(:,1)) &
                                - 2.0_dp*TriMult(mat_poly%size-1,triHorner(2)%dl,triHorner(2)%d,triHorner(2)%du,y(:,2))
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,3))
            ! compute q=eta-h^Tx, q'=eta'-h'^Tx-h^Tx', q''=eta''-h''^Tx-2h'^Tx'-h^Tx''
        else
            ! compute backward error and condition number
            berr = res/berr
            conv = 1
        end if
    end subroutine Hyman
    !****************************************************************
    !				           HymanLinearSolve                     *
    !****************************************************************
    !   Solves upper triangular system using backsubstitution. This
    !   particular system arrises when solving a system with the R
    !   amatrix of Hyman's method. In general, we reference R by its 
    !   three diagonals: d, du, du2. 
    !****************************************************************
    subroutine HymanLinearSolve(n,d,du,du2,v)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(in)    :: d(:), du(:), du2(:)
        complex(kind=dp), intent(inout) :: v(:)
        ! local variables
        integer                         :: i
        
        ! row n
        v(n) = v(n)/d(n)
        ! row n-1
        v(n-1) = (v(n-1) - du(n-1)*v(n))/d(n-1)
        ! row n-2,...,1
        do i=n-2,1,-1
            v(i) = (v(i) - du(i)*v(i+1) - du2(i)*v(i+2))/d(i)
        end do
    end subroutine HymanLinearSolve
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
        complex(kind=dp)                :: x(mat_poly%size), y(mat_poly%size), z
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
        ! note memory is allocated seperately from mat_poly%coeff
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
        do j=1,mat_poly%size*mat_poly%degree
            z = eigval(j)
            ! store normalized random eigenvector
            call zlarnv(5,iseed,mat_poly%size,eigvec(:,j))
            eigvec(:,j) = eigvec(:,j)/dznrm2(mat_poly%size,eigvec(:,j),1)
            ! update eigenvector
            if(abs(z)>1) then
                z = 1/z
                call RevEigvecUpd(mat_poly,z,eigvec(:,j))
            else
                call EigvecUpd(mat_poly,z,eigvec(:,j))
            end if
        end do
    end subroutine InitEst
    !****************************************************************
    !				           RevEigvecUpd                         *
    !****************************************************************
    !   Applies the reversal polynomial to update the eigenvector
    !   approximation given a new eigenvalue approximation.
    !****************************************************************
    subroutine RevEigvecUpd(mat_poly,zz,v)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        complex(kind=dp), intent(in)    :: zz
        complex(kind=dp), intent(inout) :: v(:)
        ! local variables
        integer                         :: k
        type(trid)                      :: mat
        ! lapack variables
        integer                         :: ipiv(mat_poly%size), info
        complex(kind=dp)                :: du2(mat_poly%size)
        ! external subroutines
        external                        :: zgtsv
        ! external functions
        real(kind=dp)                   :: dznrm2
        external                        :: dznrm2
        
        ! reversal Horner's method
        ! note mat memory is allocated seperately from mat_poly%coeff
        mat = mat_poly%coeff(1)
        do k=2,mat_poly%degree+1
            mat%du = mat%du*zz + mat_poly%coeff(k)%du
            mat%d = mat%d*zz + mat_poly%coeff(k)%d
            mat%dl = mat%dl*zz + mat_poly%coeff(k)%dl
        end do
        ! solve linear system
        call zgtsv(mat_poly%size,1,mat%dl,mat%d,mat%du,v,mat_poly%size,info)
        ! normalize
        v = v/dznrm2(mat_poly%size,v,1)
        ! deallocate mat memory
        deallocate(mat%dl,mat%d,mat%du)
    end subroutine RevEigvecUpd
    !****************************************************************
    !				           EigvecUpd                            *
    !****************************************************************
    !   Applies the original polynomial to update the eigenvector
    !   approximation given a new eigenvalue approximation.
    !****************************************************************
    subroutine EigvecUpd(mat_poly,z,v)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        complex(kind=dp), intent(in)    :: z
        complex(kind=dp), intent(inout) :: v(:)
        ! local variables
        integer                         :: k
        type(trid)                      :: mat
        ! lapack variables
        integer                         :: ipiv(mat_poly%size), info
        complex(kind=dp)                :: du2(mat_poly%size)
        ! external subroutines
        external                        :: zgtsv
        ! external functions
        real(kind=dp)                   :: dznrm2
        external                        :: dznrm2
        
        ! standard Horner's method
        ! note mat memory is allocated seperately from mat_poly%coeff
        mat = mat_poly%coeff(mat_poly%degree+1)
        do k=mat_poly%degree,1,-1
            mat%du = mat%du*z + mat_poly%coeff(k)%du
            mat%d = mat%d*z + mat_poly%coeff(k)%d
            mat%dl = mat%dl*z + mat_poly%coeff(k)%dl
        end do
        ! solve linear system
        call zgtsv(mat_poly%size,1,mat%dl,mat%d,mat%du,v,mat_poly%size,info)
        ! normalize
        v = v/dznrm2(mat_poly%size,v,1)
        ! deallocate mat memory
        deallocate(mat%dl,mat%d,mat%du)
    end subroutine EigvecUpd
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
        tmp(1) = mat%d(1)
        tmp(2) = mat%dl(1)
        call zlassq(2, tmp, 1, scale, sum)
        ! columns 2,...,n-1
        do j=2,n-1
            tmp(1) = mat%du(j-1)
            tmp(2) = mat%d(j)
            tmp(3) = mat%dl(j)
            call zlassq(3, tmp, 1, scale, sum)
        end do
        ! column n
        tmp(1) = mat%du(n-1)
        tmp(2) = mat%d(n)
        call zlassq(2, tmp, 1, scale, sum)
        ! store result
        alpha = scale*sqrt(sum)
    end subroutine FroNorm
end module TMP_EigSolver2