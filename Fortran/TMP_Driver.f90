program TMP_Driver
    use tmp_eigsolver
    implicit none
    ! local variables
    type(trid_mp)                   :: mat_poly
    complex(kind=dp), allocatable   :: eigval(:), eigvec(:,:)
    real(kind=dp), allocatable      :: berr(:), cond(:)
    integer, parameter              :: itmax = 30
    integer, allocatable            :: conv(:)
    integer                         :: k
    
    ! initiate random seed
    call InitRandomSeed()
    ! create random matrix polynomial
    mat_poly%size = 7
    mat_poly%degree = 3
    allocate(mat_poly%coeff(mat_poly%degree+1))
    do k=1,mat_poly%degree+1
        allocate(mat_poly%coeff(k)%upper(mat_poly%size-1))
        allocate(mat_poly%coeff(k)%main(mat_poly%size))
        allocate(mat_poly%coeff(k)%lower(mat_poly%size-1))
    end do
    call RandomCmplxTriMat(mat_poly)
    ! call eigenvalue solver
    allocate(eigval(mat_poly%size*mat_poly%degree))
    allocate(eigvec(mat_poly%size*mat_poly%degree,mat_poly%size))
    allocate(berr(mat_poly%size*mat_poly%degree))
    allocate(cond(mat_poly%size*mat_poly%degree))
    allocate(conv(mat_poly%size*mat_poly%degree))
    call EigSolver(mat_poly,eigval,eigvec,berr,cond,conv,itmax)
    deallocate(eigval,eigvec,berr,cond,conv)
    ! deallocate mat_poly memory
    do k=1,mat_poly%degree+1
        deallocate(mat_poly%coeff(k)%upper)
        deallocate(mat_poly%coeff(k)%main)
        deallocate(mat_poly%coeff(k)%lower)
    end do
    deallocate(mat_poly%coeff)
contains
    !****************************************************************
    !				           InitRandomSeed                       *
    !****************************************************************
    !   Initiate random seed using system_clock. This seed is then 
    !   available for the random number generator in random_number for
    !   the life of the program. 
    !****************************************************************
    !subroutine InitRandomSeed()
    !    implicit none
    !    ! local variables
    !    integer                             :: i, n, clock
    !    integer, dimension(:), allocatable  :: seed
    !    ! intrinsic subroutines
    !    intrinsic                           :: random_seed, system_clock
    !    
    !    ! main
    !    call random_seed(size = n)
    !    allocate(seed(n))
    !    
    !    call system_clock(count = clock)
    !    seed = clock + 37 * (/ (i - 1, i = 1,n) /)
    !    call random_seed(put = seed)
    !    
    !    deallocate(seed)
    !end subroutine InitRandomSeed
    !****************************************************************
    !				           RandCmplxTriMP                       *
    !****************************************************************
    !   Creates Random Complex Tridiagonal Matrix Polynomial where
    !   the scalar coefficients are complex numbers with real and
    !   imaginary parts uniformly distributed in (-1,1).
    !****************************************************************
    subroutine RandomCmplxTriMat(mat_poly)
        implicit none
        ! argument variables
        type(trid_mp)           :: mat_poly
        ! local variables
        integer                 :: i, j, k
        
        ! loop through all coefficient matrices
        do k=1,mat_poly%degree+1
            ! column 1 of coeff(k)
            call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%main(1))
            call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%lower(1))
            ! colmns 2,...,mat_poly%size-1 of coeff(k)
            do j=2,mat_poly%size-1
                call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%upper(j-1))
                call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%main(j))
                call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%lower(j))
            end do
            ! column mat_poly%size of coeff(k)
            call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%upper(mat_poly%size-1))
            call CmplxRandomNumber(-1.0_dp,1.0_dp,mat_poly%coeff(k)%main(mat_poly%size))
        end do
    end subroutine RandomCmplxTriMat
    !****************************************************************
    !				           cmplx_random_number                  *
    !****************************************************************
    !   Creates complex random number with real and imaginary part
    !   uniformly distributed in (a,b).
    !****************************************************************
    !subroutine CmplxRandomNumber(a,b,res)
    !    implicit none
    !    ! argument variables
    !    real(kind=dp)           :: a, b
    !    complex(kind=dp)        :: res
    !    ! local variables
    !    real(kind=dp)           :: r1, r2
    !    
    !    ! call random_number
    !    call random_number(r1)
    !    call random_number(r2)
    !    ! store result
    !    res = cmplx(a + (b-a)*r1,a + (b-a)*r2,kind=dp)
    !end subroutine CmplxRandomNumber
end program TMP_Driver