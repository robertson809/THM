program TMP_Driver
    use tmp_eigsolver3
    implicit none
    ! local variables
    type(trid_mp)                   :: mat_poly
    integer                         :: clock, j, k
    ! lapack variables
    integer                         :: iseed(4)
    ! eigsolver variables
    complex(kind=dp), allocatable   :: eigval(:), eigvec(:,:)
    real(kind=dp), allocatable      :: berr(:), cond(:)
    integer, parameter              :: itmax = 50
    integer, allocatable            :: conv(:)
    ! intrinsic subroutines
    intrinsic                       :: system_clock
    ! external subroutines
    external                        :: zlarnv
    
    ! initial iseed
    call system_clock(count = clock)
    iseed = clock + 37*(/ (j-1, j=1,4) /)
    iseed = (/ (mod(iseed(j),4095), j=1,4)/)
    if(mod(iseed(4),2)==0) iseed(4)=iseed(4)+1
    ! allocate memory for mat_poly
    mat_poly%size = 100
    mat_poly%degree = 2
    allocate(mat_poly%coeff(mat_poly%degree+1))
    do k=1,mat_poly%degree+1
        allocate(mat_poly%coeff(k)%du(mat_poly%size-1))
        allocate(mat_poly%coeff(k)%d(mat_poly%size))
        allocate(mat_poly%coeff(k)%dl(mat_poly%size-1))
    end do
    ! create random tridiagonal matrix polynomial
    do k=1,mat_poly%degree+1
        call zlarnv(5,iseed,mat_poly%size-1,mat_poly%coeff(k)%du)
        call zlarnv(5,iseed,mat_poly%size,mat_poly%coeff(k)%d)
        call zlarnv(5,iseed,mat_poly%size-1,mat_poly%coeff(k)%dl)
    end do
    ! call eigenvalue solver
    allocate(eigval(mat_poly%size*mat_poly%degree))
    allocate(eigvec(mat_poly%size,mat_poly%size*mat_poly%degree))
    allocate(berr(mat_poly%size*mat_poly%degree))
    allocate(cond(mat_poly%size*mat_poly%degree))
    allocate(conv(mat_poly%size*mat_poly%degree))
    call EigSolver(mat_poly,eigval,eigvec,berr,cond,conv,itmax)
    deallocate(eigval,eigvec,berr,cond,conv)
    ! deallocate mat_poly memory
    do k=1,mat_poly%degree+1
        deallocate(mat_poly%coeff(k)%du)
        deallocate(mat_poly%coeff(k)%d)
        deallocate(mat_poly%coeff(k)%dl)
    end do
    deallocate(mat_poly%coeff)
end program TMP_Driver