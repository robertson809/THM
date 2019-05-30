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
        integer                         :: i, j, num_eig, total_eig
        real(kind=dp)                   :: r
        real(kind=dp), allocatable      :: alpha(:)
        complex(kind=dp)                :: lag_term1, lag_term2, z
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! store norm of matrix coefficients
        allocate(alpha(mat_poly%degree+1))
        do i=1,mat_poly%degree+1
            call FroNorm(mat_poly%size,mat_poly%coeff(i),alpha(i))
        end do
        ! initial estimates
        call InitEst(mat_poly,eigval,eigvec)
        ! main loop
        num_eig = 0
        total_eig = mat_poly%size*mat_poly%degree
        conv = (/ (0, j=1,total_eig) /)
        do i=1,itmax
            do j=1,total_eig
                if(conv(j)==0) then
                    z = eigval(j)
                    r = abs(z)
                    if(r>1) then
                        ! reversal Hyman to compute lag terms
                    else
                        ! standard Hyman to compute lag terms
                        write(*,*) 'standard Hyman'
                        call Hyman(mat_poly,alpha,z,r,eigvec(:,j),lag_term1,lag_term2,berr(j),cond(j),conv(j))
                        stop
                    end if
                    if(conv(j)==0) then
                        ! modify lag terms and update eigenvalue approximation
                        ! update eigenvector approximation
                    else
                        ! don't update eigenvalue approximation
                        num_eig = num_eig + 1
                        if(num_eig==total_eig) go to 10
                    end if
                end if
            end do
        end do
        ! final steps
        10 continue
        ! deallocate alpha memory
        deallocate(alpha)
    end subroutine EigSolver
    !****************************************************************
    !				           RevHyman                             *
    !****************************************************************
    !   Applies the reversal polynomial and Hyman's method to compute
    !   the Laguerre correction terms. In the process, the backward
    !   error is computed and stopping criterion is checked. If the
    !   stopping criterion is true, then conv is set to 1 and the
    !   eigenvalue approximation will not be updated.
    !****************************************************************
    subroutine RevHyman(mat_poly,alpha,z,v,lag_term1,lag_term2,berr,cond,conv)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        integer, intent(out)            :: conv
        real(kind=dp), intent(in)       :: alpha(:), berr, cond
        complex(kind=dp), intent(in)    :: z, v(:)
        complex(kind=dp), intent(out)   :: lag_term1, lag_term2
        ! local variables
        type(trid)                      :: triHorner(3)
            
        ! main
    end subroutine RevHyman
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
        real(kind=dp), intent(in)       :: alpha(:), r
        real(kind=dp), intent(out)      :: berr, cond
        complex(kind=dp), intent(in)    :: z, v(:)
        complex(kind=dp), intent(out)   :: lag_term1, lag_term2
        ! local variables
        integer                         :: k
        type(trid)                      :: triHorner(3)
        ! external function
        external                        :: dznrm2
        real(kind=dp)                   :: dznrm2
            
        ! call Horner's method
        call Horner(mat_poly,z,triHorner)
        ! compute backward error
        berr = alpha(mat_poly%degree+1)
        do k=mat_poly%degree,1,-1
            berr = r*berr + alpha(k)
        end do
        berr = dznrm2(mat_poly%size,TriMult(mat_poly%size,triHorner(1),v),1)
        write(*,*) berr
        ! check for subdiagonal zeros in triHorner(1)
        do k=1,mat_poly%size-1
            if(abs(triHorner(1)%lower(k)) .le. small) triHorner(1)%lower(k) = small
        end do
        
        ! deallocate triHorner memory (allocated in Horner)
        do k=1,3
            deallocate(triHorner(k)%upper,triHorner(k)%main,triHorner(k)%lower)
        end do
    end subroutine Hyman
    !****************************************************************
    !				           Horner                               *
    !****************************************************************
    !   Applies the standard Horner method to evaluate a matrix 
    !   polynomial and its derivatives. 
    !****************************************************************
    subroutine Horner(mat_poly,z,triHorner)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        complex(kind=dp), intent(in)    :: z
        type(trid), intent(out)         :: triHorner(3)
        ! local variables
        integer                         :: k
        
        ! polynomial evaluation
        triHorner(1) = mat_poly%coeff(mat_poly%degree+1)
        
        ! polynomial derivative evaluation
        allocate(triHorner(2)%upper(mat_poly%size -1), triHorner(2)%main(mat_poly%size), triHorner(2)%lower(mat_poly%size -1))
        triHorner(2)%upper = cmplx(0,0,kind=dp)
        triHorner(2)%main = cmplx(0,0,kind=dp)
        triHorner(2)%lower = cmplx(0,0,kind=dp)
        
        ! polynomial 2nd derivative evaluation
        allocate(triHorner(3)%upper(mat_poly%size -1), triHorner(3)%main(mat_poly%size), triHorner(3)%lower(mat_poly%size -1))
        triHorner(3)%upper = cmplx(0,0,kind=dp)
        triHorner(3)%main = cmplx(0,0,kind=dp)
        triHorner(3)%lower = cmplx(0,0,kind=dp)
        
        ! Horner's method
        do k=mat_poly%degree,1,-1
            ! update polynomial 2nd derivative evaluation
            triHorner(3)%upper = triHorner(3)%upper*z + triHorner(2)%upper
            triHorner(3)%main = triHorner(3)%main*z + triHorner(2)%main
            triHorner(3)%lower = triHorner(3)%lower*z + triHorner(2)%lower
            
            ! update polynomial derivative evaluation
            triHorner(2)%upper = triHorner(2)%upper*z + triHorner(1)%upper
            triHorner(2)%main = triHorner(2)%main*z + triHorner(1)%main
            triHorner(2)%lower = triHorner(2)%lower*z + triHorner(1)%lower
            
            ! update polynomial evaluation
            triHorner(1)%upper = triHorner(1)%upper*z + mat_poly%coeff(k)%upper
            triHorner(1)%main = triHorner(1)%main*z + mat_poly%coeff(k)%main
            triHorner(1)%lower = triHorner(1)%lower*z + mat_poly%coeff(k)%lower
        end do
    end subroutine Horner
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
        complex(kind=dp)                :: qmat(mat_poly%size-1,2), x(mat_poly%size), z
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
        ! intitial estimate eigenvectors
        do j=1,mat_poly%size*mat_poly%degree
            z = eigval(j)
            ! evaluate matrix polynomial
            if(abs(z)>1) then
                ! reversal Horner's method
                z = 1/z
                mat = mat_poly%coeff(1)
                do k=2,mat_poly%degree+1
                    mat%upper = mat%upper*z + mat_poly%coeff(k)%upper
                    mat%main = mat%main*z + mat_poly%coeff(k)%main
                    mat%lower = mat%lower*z + mat_poly%coeff(k)%lower
                end do
            else
                ! standard Horner's method
                mat = mat_poly%coeff(mat_poly%degree+1)
                do k=mat_poly%degree,1,-1
                    mat%upper = mat%upper*z + mat_poly%coeff(k)%upper
                    mat%main = mat%main*z + mat_poly%coeff(k)%main
                    mat%lower = mat%lower*z + mat_poly%coeff(k)%lower
                end do
            end if
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
        
        ! rotations 1,...,n-1
        do i=1,n-1
            ! apply rotation to (x(i),x(i+1))
            call ApplyRot(qmat(i,1),qmat(i,2),x(i),x(i+1))
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
        real(kind=dp)                   :: r
        
        ! rotations 1,...,n-2
        do i=1,n-2
            ! build rotation for (main(i),lower(i))
            call BuildRot(mat%main(i),mat%lower(i),qmat(i,1),qmat(i,2),r)
            ! apply rotation to (main(i),lower(i))
            mat%main(i) = r
            mat%lower(i) = cmplx(0,0,kind=dp)
            ! apply rotation to (upper(i),main(i+1))
            call ApplyRot(qmat(i,1),qmat(i,2),mat%upper(i),mat%main(i+1))
            ! apply rotation to (0,upper(i+1))
            call ApplyRot(qmat(i,1),qmat(i,2),mat%lower(i),mat%upper(i+1))
        end do
        ! rotation n-1
        i = n-1
        ! build rotation for (main(i),lower(i))
        call BuildRot(mat%main(i),mat%lower(i),qmat(i,1),qmat(i,2),r)
        ! apply rotation to (main(i),lower(i))
        mat%main(i) = r
        mat%lower(i) = cmplx(0,0,kind=dp)
        ! apply rotation to (upper(i),main(i+1))
        call ApplyRot(qmat(i,1),qmat(i,2),mat%upper(i),mat%main(i+1))
    end subroutine TriQR
    !****************************************************************
    !				           BuildRot                             *
    !****************************************************************
    !   Build rotation transformation for pair (x,y).
    !****************************************************************
    subroutine BuildRot(x,y,c,s,r)
        implicit none
        ! argument variables
        complex(kind=dp), intent(in)    :: x, y
        complex(kind=dp), intent(out)   :: c, s
        real(kind=dp), intent(out)      :: r
        ! local variables
        real(kind=dp)                   :: m
        complex(kind=dp)                :: tmp_x, tmp_y
        ! intrinsic functions
        intrinsic                       :: abs, max, sqrt
        
        ! maximum abs value
        m = max(abs(x),abs(y))
        ! compute rotation (c,s) and norm r
        if(m .le. small) then
            c = cmplx(1,0,kind=dp)
            s = cmplx(0,0,kind=dp)
            r = cmplx(0,0,kind=dp)
        else
            tmp_x = x/m
            tmp_y = y/m
            r = sqrt(abs(tmp_x)**2+abs(tmp_y)**2)
            c = tmp_x/r
            s = tmp_y/r
            r = m*r
        end if
    end subroutine BuildRot
    !****************************************************************
    !				           ApplyRot                             *
    !****************************************************************
    !   Apply rotation stored in (c,s) to pair (x,y).
    !****************************************************************
    subroutine ApplyRot(c,s,x,y)
        implicit none
        ! argument variables
        complex(kind=dp), intent(in)    :: c, s
        complex(kind=dp), intent(inout) :: x, y
        ! local variables
        complex(kind=dp)                :: tmp(2)
        ! intrinsic functions
        intrinsic                       :: conjg
        
        ! apply rotation to (x,y)
        tmp(1) = conjg(c)*x + conjg(s)*y
        tmp(2) = -s*x + c*y
        ! update (x,y)
        x = tmp(1)
        y = tmp(2)
    end subroutine ApplyRot
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
    !****************************************************************
    !				           InitRandomSeed                       *
    !****************************************************************
    !   Initiate random seed using system_clock. This seed is then 
    !   available for the random number generator in random_number for
    !   the life of the program. 
    !****************************************************************
    subroutine InitRandomSeed()
        implicit none
        ! local variables
        integer                             :: i, n, clock
        integer, dimension(:), allocatable  :: seed
        ! intrinsic subroutines
        intrinsic                           :: random_seed, system_clock
        
        ! main
        call random_seed(size = n)
        allocate(seed(n))
        
        call system_clock(count = clock)
        seed = clock + 37 * (/ (i - 1, i = 1,n) /)
        call random_seed(put = seed)
        
        deallocate(seed)
    end subroutine InitRandomSeed
    !****************************************************************
    !				           RandCmplxVec                         *
    !****************************************************************
    !   Creates Random Complex vector whose entries are complex
    !   numbers with real and imaginary parts uniformly distributed
    !   in the interval (-1,1).
    !****************************************************************
    subroutine RandomCmplxVec(vec)
        implicit none
        ! argument variables
        integer                     :: n
        complex(kind=dp)            :: vec(:)
        ! local variables
        integer                     :: i
        
        ! create random vector
        do i=1,n
            call CmplxRandomNumber(-1.0_dp,1.0_dp,vec(i))
        end do
    end subroutine RandomCmplxVec
    !****************************************************************
    !				           cmplx_random_number                  *
    !****************************************************************
    !   Creates complex random number with real and imaginary part
    !   uniformly distributed in (a,b).
    !****************************************************************
    subroutine CmplxRandomNumber(a,b,res)
        implicit none
        ! argument variables
        real(kind=dp)           :: a, b
        complex(kind=dp)        :: res
        ! local variables
        real(kind=dp)           :: r1, r2
        
        ! call random_number
        call random_number(r1)
        call random_number(r2)
        ! store result
        res = cmplx(a + (b-a)*r1,a + (b-a)*r2,kind=dp)
    end subroutine CmplxRandomNumber
end module TMP_EigSolver