module TMP_EigSolver3
    use fpml
    implicit none
    real(kind=dp), parameter            :: mu=2.0_dp**(-53)
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
        real(kind=dp)                   :: alpha(mat_poly%degree+1)
        complex(kind=dp)                :: lag_term1, lag_term2, z
        ! intrinsic functions
        intrinsic                       :: abs
        
        ! open file to store eigenvalue estimates
        open(unit=1,file="data_files/eigvals.csv")
        ! store norm of matrix coefficients
        do i=1,mat_poly%degree+1
            alpha(i) = FroNorm(mat_poly%size,mat_poly%coeff(i))*(3.8_dp*(i-1)+1.0_dp)
        end do
        ! initial estimates
        call InitEst(mat_poly,eigval)
        ! write eigenvalue estimates to file
        total_eig = mat_poly%size*mat_poly%degree
        do j=1,total_eig-1
            if(aimag(eigval(j))>0) then
                write(1,'(F10.3,A,F10.3,A)', advance='no') real(eigval(j)), ' + ', abs(aimag(eigval(j))), '*1j, '
            else
                write(1,'(F10.3,A,F10.3,A)', advance='no') real(eigval(j)), ' - ', abs(aimag(eigval(j))), '*1j, '
            end if
        end do
        j=total_eig
        if(aimag(eigval(j))>0) then
            write(1,'(F10.3,A,F10.3,A)') real(eigval(j)), ' + ', abs(aimag(eigval(j))), '*1j'
        else
            write(1,'(F10.3,A,F10.3,A)') real(eigval(j)), ' - ', abs(aimag(eigval(j))), '*1j'
        end if
        ! main loop
        conv = 0
        num_eig = 0
        do i=1,itmax
            do j=1,total_eig
                if(conv(j)==0) then
                    z = eigval(j)
                    r = abs(z)
                    if(r>1.0_dp) then
                        ! reversal Hyman to compute lag terms
                        call RevHyman(mat_poly,alpha,z,r,lag_term1,lag_term2,berr(j),conv(j),total_eig)
                    else
                        ! standard Hyman to compute lag terms
                        call Hyman(mat_poly,alpha,z,r,lag_term1,lag_term2,berr(j),conv(j))
                    end if
                    if(conv(j)==0) then
                        ! modify lag terms and update eigenvalue approximation
                        call ModifyLaguerre(total_eig,lag_term1,lag_term2,z,j,eigval)
                        eigval(j) = z - lag_term2
                    else
                        ! compute eigenvectors, backward error, and condition numbers
                        num_eig = num_eig + 1
                        if(num_eig==total_eig) then
                            write(*,*) 'itnum = ', i
                            go to 10
                        end if
                    end if
                end if
            end do
            ! write eigenvalue estimates to file
            do j=1,total_eig-1
                if(aimag(eigval(j))>0) then
                    write(1,'(F10.3,A,F10.3,A)', advance='no') real(eigval(j)), ' + ', abs(aimag(eigval(j))), '*1j, '
                else
                    write(1,'(F10.3,A,F10.3,A)', advance='no') real(eigval(j)), ' - ', abs(aimag(eigval(j))), '*1j, '
                end if
            end do
            j=total_eig
            if(aimag(eigval(j))>0) then
                write(1,'(F10.3,A,F10.3,A)') real(eigval(j)), ' + ', abs(aimag(eigval(j))), '*1j'
            else
                write(1,'(F10.3,A,F10.3,A)') real(eigval(j)), ' - ', abs(aimag(eigval(j))), '*1j'
            end if
        end do
        ! final steps
        10 continue
        write(*,*) 'berr = ', maxval(berr)
        ! write last eigenvalue estimates to file
        do j=1,total_eig-1
            if(aimag(eigval(j))>0) then
                write(1,'(F10.3,A,F10.3,A)', advance='no') real(eigval(j)), ' + ', abs(aimag(eigval(j))), '*1j, '
            else
                write(1,'(F10.3,A,F10.3,A)', advance='no') real(eigval(j)), ' - ', abs(aimag(eigval(j))), '*1j, '
            end if
        end do
        j=total_eig
        if(aimag(eigval(j))>0) then
            write(1,'(F10.3,A,F10.3,A)') real(eigval(j)), ' + ', abs(aimag(eigval(j))), '*1j, '
        else
            write(1,'(F10.3,A,F10.3,A)') real(eigval(j)), ' - ', abs(aimag(eigval(j))), '*1j, '
        end if
        ! close file
        close(1)
    end subroutine EigSolver
    !****************************************************************
    !				           Modify Laguerre                      *
    !****************************************************************
    !   Modifys the Laguerre correction terms computed by Hyman's 
    !   method by taking into account the deflation strategy. The
    !   final Laguerre correction term is returned in lag_term2.
    !****************************************************************
    subroutine ModifyLaguerre(total_eig,lag_term1,lag_term2,z,j,eigval)
        implicit none
        ! argument variables
        integer, intent(in)             :: total_eig, j
        complex(kind=dp), intent(in)    :: eigval(:), z
        complex(kind=dp), intent(inout) :: lag_term1, lag_term2
        ! local variables
        integer                         :: k
        complex(kind=dp)                :: temp
        ! intrinsic functions
        intrinsic                       :: abs, sqrt
        
        ! deflation strategy
        do k=1,j-1
            temp = 1.0_dp/(z - eigval(k))
            lag_term1 = lag_term1 - temp
            lag_term2 = lag_term2 - temp**2
        end do
        do k=j+1,total_eig
            temp = 1.0_dp/(z - eigval(k))
            lag_term1 = lag_term1 - temp
            lag_term2 = lag_term2 - temp**2
        end do
        ! Laguerre correction
        temp = sqrt((total_eig-1)*(total_eig*lag_term2-lag_term1**2))
        lag_term2 = lag_term1 + temp
        lag_term1 = lag_term1 - temp
        if(abs(lag_term1)>abs(lag_term2)) then
            lag_term2 = total_eig/lag_term1
        else
            lag_term2 = total_eig/lag_term2
        end if
    end subroutine ModifyLaguerre
    !****************************************************************
    !				           RevHyman                             *
    !****************************************************************
    !   Applies the reversal polynomial and Hyman's method to compute
    !   the Laguerre correction terms. In the process, the backward
    !   error is computed and stopping criterion is checked. If the
    !   stopping criterion is true, then conv is set to 1 and the
    !   eigenvalue approximation will not be updated.
    !****************************************************************
    subroutine RevHyman(mat_poly,alpha,z,r,lag_term1,lag_term2,berr,conv,total_eig)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        integer, intent(in)             :: total_eig
        integer, intent(out)            :: conv
        real(kind=dp), intent(in)       :: alpha(mat_poly%degree+1), r
        real(kind=dp), intent(out)      :: berr
        complex(kind=dp), intent(in)    :: z
        complex(kind=dp), intent(out)   :: lag_term1, lag_term2
        ! local variables
        integer                         :: k
        type(trid)                      :: mat, triHorner(3)
        real(kind=dp)                   :: rcond, rr
        complex(kind=dp)                :: y(mat_poly%size-1,3), q(3), sum, zz
        ! lapack variables
        complex(kind=dp)                :: du2(mat_poly%size-2), work(mat_poly%size*2)
        integer                         :: ipiv(mat_poly%size), info
        ! external subroutine
        external                        :: zgttrf, zgtcon
        ! external function
        real(kind=dp)                   :: dznrm2
        external                        :: dznrm2
        
        ! reversal polynomial evaluation (note triHorner(1) gets allocated memory seperately from mat_poly%coeff)
        zz = 1/z
        triHorner(1) = mat_poly%coeff(1)
        ! polynomial derivative evaluation
        allocate(triHorner(2)%du(mat_poly%size-1),triHorner(2)%d(mat_poly%size),triHorner(2)%dl(mat_poly%size-1))
        triHorner(2)%du = cmplx(0,0,kind=dp)
        triHorner(2)%d = cmplx(0,0,kind=dp)
        triHorner(2)%dl = cmplx(0,0,kind=dp)
        ! polynomial 2nd derivative evaluation
        allocate(triHorner(3)%du(mat_poly%size-1),triHorner(3)%d(mat_poly%size),triHorner(3)%dl(mat_poly%size-1))
        triHorner(3)%du = cmplx(0,0,kind=dp)
        triHorner(3)%d = cmplx(0,0,kind=dp)
        triHorner(3)%dl = cmplx(0,0,kind=dp)
        ! Horner's method
        do k=2,mat_poly%degree+1
            ! update polynomial 2nd derivative evaluation
            triHorner(3)%du = triHorner(3)%du*zz + triHorner(2)%du
            triHorner(3)%d = triHorner(3)%d*zz + triHorner(2)%d
            triHorner(3)%dl = triHorner(3)%dl*zz + triHorner(2)%dl    
            ! update polynomial derivative evaluation
            triHorner(2)%du = triHorner(2)%du*zz + triHorner(1)%du
            triHorner(2)%d = triHorner(2)%d*zz + triHorner(1)%d
            triHorner(2)%dl = triHorner(2)%dl*zz + triHorner(1)%dl
            ! update polynomial evaluation
            triHorner(1)%du = triHorner(1)%du*zz + mat_poly%coeff(k)%du
            triHorner(1)%d = triHorner(1)%d*zz + mat_poly%coeff(k)%d
            triHorner(1)%dl = triHorner(1)%dl*zz + mat_poly%coeff(k)%dl
        end do
        ! scale 2nd derivative 
        triHorner(3)%du = 2.0_dp*triHorner(3)%du
        triHorner(3)%d = 2.0_dp*triHorner(3)%d
        triHorner(3)%dl = 2.0_dp*triHorner(3)%dl
        ! backward error
        rr = 1/r
        berr = alpha(1)
        do k=2,mat_poly%degree+1
            berr = rr*berr + alpha(k)
        end do
        mat = triHorner(1)
        call zgttrf(mat_poly%size,mat%dl,mat%d,mat%du,du2,ipiv,info)
        call zgtcon('1',mat_poly%size,mat%dl,mat%d,mat%du,du2,ipiv,berr,rcond,work,info)
        berr = rcond
        !  deallocate mat memory (allocated above)
        deallocate(mat%du,mat%d,mat%dl)
        if(berr > eps) then
            ! check for subdiagonal zeros in triHorner(1)
            do k=1,mat_poly%size-1
                if(abs(triHorner(1)%dl(k)) .le. mu) then
                    write(*,*) 'Hyman: small subdiagonal'
                    triHorner(1)%dl(k) = mu
                end if
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
            y(:,2) = y(:,2) - HymanMatMult(mat_poly%size-1,triHorner(2)%dl(:),triHorner(2)%d(2:mat_poly%size-1), &
                                triHorner(2)%du(2:mat_poly%size-2),y(:,1))
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,2))
            ! solve linear system Rx''=y''-R''x-2R'x'
            y(:,3) = y(:,3) - HymanMatMult(mat_poly%size-1,triHorner(3)%dl(:),triHorner(3)%d(2:mat_poly%size-1), &
                                triHorner(3)%du(2:mat_poly%size-2),y(:,1)) - 2.0_dp*HymanMatMult(mat_poly%size-1, &
                                triHorner(2)%dl(:),triHorner(2)%d(2:mat_poly%size-1),triHorner(2)%du(2:mat_poly%size-2),y(:,2))
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,3))
            ! compute q=-h^Tx
            q(1) = -triHorner(1)%d(1)*y(1,1)-triHorner(1)%du(1)*y(2,1)
            ! compute q'=-h'^Tx-h^Tx'
            q(2) = -triHorner(2)%d(1)*y(1,1)-triHorner(2)%du(1)*y(2,1) &
                    -triHorner(1)%d(1)*y(1,2)-triHorner(1)%du(1)*y(2,2)
            ! compute q''=-h''^tx-2h'^Tx'-h^Tx''
			q(3) = -triHorner(3)%d(1)*y(1,1)-triHorner(3)%du(1)*y(2,1) &
                    -2.0_dp*(triHorner(2)%d(1)*y(1,2)+triHorner(2)%du(1)*y(2,2)) &
                    -triHorner(1)%d(1)*y(1,3)-triHorner(1)%du(1)*y(2,3)
            ! laguerre terms
            lag_term1 = q(2)/q(1)
            lag_term2 = lag_term1**2 - q(3)/q(1)
            ! update laguerre terms with log derivatives
            do k=1,mat_poly%size-1
                sum = triHorner(2)%dl(k)/triHorner(1)%dl(k)
                lag_term1 = lag_term1 + sum
                lag_term2 = lag_term2 - triHorner(3)%dl(k)/triHorner(1)%dl(k) + sum**2
            end do
            ! modify for reversal
            lag_term2 = zz**2*(total_eig-2.0_dp*zz*lag_term1+zz**2*lag_term2)
            lag_term1 = zz*(total_eig-zz*lag_term1)
        else
            ! stop updating eigenvalue
            conv = 1
        end if
        ! deallocate triHorner memory (allocated in Horner)
        do k=1,3
            deallocate(triHorner(k)%du,triHorner(k)%d,triHorner(k)%dl)
        end do
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
    subroutine Hyman(mat_poly,alpha,z,r,lag_term1,lag_term2,berr,conv)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        integer, intent(out)            :: conv
        real(kind=dp), intent(in)       :: alpha(mat_poly%degree+1), r
        real(kind=dp), intent(out)      :: berr
        complex(kind=dp), intent(in)    :: z
        complex(kind=dp), intent(out)   :: lag_term1, lag_term2
        ! local variables
        integer                         :: k
        type(trid)                      :: mat, triHorner(3)
        real(kind=dp)                   :: rcond
        complex(kind=dp)                :: y(mat_poly%size-1,3), q(3), sum
        ! lapack variables
        complex(kind=dp)                :: du2(mat_poly%size-2), work(mat_poly%size*2)
        integer                         :: ipiv(mat_poly%size), info
        ! external subroutine
        external                        :: zgttrf, zgtcon
        ! external function
        real(kind=dp)                   :: dznrm2
        external                        :: dznrm2
        
        ! polynomial evaluation (note triHorner(1) gets allocated memory seperately from mat_poly%coeff)
        triHorner(1) = mat_poly%coeff(mat_poly%degree+1)
        ! polynomial derivative evaluation
        allocate(triHorner(2)%du(mat_poly%size-1),triHorner(2)%d(mat_poly%size),triHorner(2)%dl(mat_poly%size-1))
        triHorner(2)%du = cmplx(0,0,kind=dp)
        triHorner(2)%d = cmplx(0,0,kind=dp)
        triHorner(2)%dl = cmplx(0,0,kind=dp)
        ! polynomial 2nd derivative evaluation
        allocate(triHorner(3)%du(mat_poly%size-1),triHorner(3)%d(mat_poly%size),triHorner(3)%dl(mat_poly%size-1))
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
        ! scale 2nd derivative 
        triHorner(3)%du = 2.0_dp*triHorner(3)%du
        triHorner(3)%d = 2.0_dp*triHorner(3)%d
        triHorner(3)%dl = 2.0_dp*triHorner(3)%dl
        ! backward error
        berr = alpha(mat_poly%degree+1)
        do k=mat_poly%degree,1,-1
            berr = r*berr + alpha(k)
        end do
        mat = triHorner(1)
        call zgttrf(mat_poly%size,mat%dl,mat%d,mat%du,du2,ipiv,info)
        call zgtcon('1',mat_poly%size,mat%dl,mat%d,mat%du,du2,ipiv,berr,rcond,work,info)
        berr = rcond
        !  deallocate mat memory (allocated above)
        deallocate(mat%du,mat%d,mat%dl)
        if(berr > eps) then
            ! check for subdiagonal zeros in triHorner(1)
            do k=1,mat_poly%size-1
                if(abs(triHorner(1)%dl(k)) .le. mu) then
                    write(*,*) 'Hyman: small subdiagonal'
                    triHorner(1)%dl(k) = mu
                end if
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
            y(:,2) = y(:,2) - HymanMatMult(mat_poly%size-1,triHorner(2)%dl(:),triHorner(2)%d(2:mat_poly%size-1), &
                                triHorner(2)%du(2:mat_poly%size-2),y(:,1))
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,2))
            ! solve linear system Rx''=y''-R''x-2R'x'
            y(:,3) = y(:,3) - HymanMatMult(mat_poly%size-1,triHorner(3)%dl(:),triHorner(3)%d(2:mat_poly%size-1), &
                                triHorner(3)%du(2:mat_poly%size-2),y(:,1)) - 2.0_dp*HymanMatMult(mat_poly%size-1, &
                                triHorner(2)%dl(:),triHorner(2)%d(2:mat_poly%size-1),triHorner(2)%du(2:mat_poly%size-2),y(:,2))
            call HymanLinearSolve(mat_poly%size-1,triHorner(1)%dl(:),triHorner(1)%d(2:mat_poly%size-1), &
                                triHorner(1)%du(2:mat_poly%size-2),y(:,3))
            ! compute q=-h^Tx
            q(1) = -triHorner(1)%d(1)*y(1,1)-triHorner(1)%du(1)*y(2,1)
            ! compute q'=-h'^Tx-h^Tx'
			q(2) = -triHorner(2)%d(1)*y(1,1)-triHorner(2)%du(1)*y(2,1) &
                    -triHorner(1)%d(1)*y(1,2)-triHorner(1)%du(1)*y(2,2)
            ! compute q''=-h''^tx-2h'^Tx'-h^Tx''
			q(3) = -triHorner(3)%d(1)*y(1,1)-triHorner(3)%du(1)*y(2,1) &
                    -2.0_dp*(triHorner(2)%d(1)*y(1,2)+triHorner(2)%du(1)*y(2,2)) &
                    -triHorner(1)%d(1)*y(1,3)-triHorner(1)%du(1)*y(2,3)
            ! laguerre terms
            lag_term1 = q(2)/q(1)
            lag_term2 = lag_term1**2 - q(3)/q(1)
            ! update laguerre terms with log derivatives
            do k=1,mat_poly%size-1
                sum = triHorner(2)%dl(k)/triHorner(1)%dl(k)
                lag_term1 = lag_term1 + sum
                lag_term2 = lag_term2 - triHorner(3)%dl(k)/triHorner(1)%dl(k) + sum**2
            end do
        else
            ! stop updating eigenvalue
            conv = 1
        end if
        ! deallocate triHorner memory (allocated in Horner)
        do k=1,3
            deallocate(triHorner(k)%du,triHorner(k)%d,triHorner(k)%dl)
        end do
    end subroutine Hyman
    !****************************************************************
    !				           HymanMatMult                         *   
    !****************************************************************
    !   Performs matrix multiplication with the matrix R in Hymans
    !   method and an arbitrary vector. We reference R by its three
    !   diagonals: d, du, du2.
    !****************************************************************
    function HymanMatMult(n,d,du,du2,x) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        complex(kind=dp), intent(in)    :: d(:), du(:), du2(:), x(:)
        ! local variables
        integer                         :: i
        complex(kind=dp)                :: res(n)
        
        ! rows 1,...,n-2
        do i=1,n-2
            res(i) = d(i)*x(i)+du(i)*x(i+1)+du2(i)*x(i+2)
        end do
        ! row n-1
        res(n-1) = d(n-1)*x(n-1)+du(n-1)*x(n)
        ! row n
        res(n) = d(n)*x(n)
        return
    end function HymanMatMult
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
    !   Computes the initial eigenvalue estimates for the tridiagonal 
    !   matrix polynomial stored in mat_poly.
    !****************************************************************
    subroutine InitEst(mat_poly,eigval)
        implicit none
        ! argument variables
        type(trid_mp), intent(in)       :: mat_poly
        complex(kind=dp), intent(out)   :: eigval(:)
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
        ! testing: random initial estimates
        ! call zlarnv(5,iseed,mat_poly%size*mat_poly%degree,eigval)
        ! return
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
            do k=1,mat_poly%degree
                ! multiply random vector by kth coefficient
                y = TriMult(mat_poly%size,mat_poly%coeff(k)%dl,mat_poly%coeff(k)%d,mat_poly%coeff(k)%du,x)
                ! solve linear system Am*x=y
                call zgttrs('N',mat_poly%size,1,mat%dl,mat%d,mat%du,du2,ipiv,y,mat_poly%size,info)
                ! compute dot product x^H*y
                poly(k) = zdotc(mat_poly%size,x,1,y,1)
            end do
            poly(mat_poly%degree+1) = cmplx(1,0,kind=dp)
            ! compute roots of polynomial
            call fpml_main(poly,mat_poly%degree,eigval((j-1)*mat_poly%degree+1:j*mat_poly%degree),berr,cond,conv,itmax)
        end do
        ! deallocate mat memory
         deallocate(mat%dl,mat%d,mat%du)
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
    !****************************************************************
    !				           FroNorm                              *
    !****************************************************************
    !   Computes the Frobenius Norm of a tridiagonal matrix.
    !****************************************************************
    function FroNorm(n,mat) result(res)
        implicit none
        ! argument variables
        integer, intent(in)             :: n
        type(trid), intent(in)          :: mat
        ! local variables
        integer                         :: j
        real(kind=dp)                   :: res, scale, sum
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
        call zlassq(2, tmp(1), 1, scale, sum)
        ! columns 2,...,n-1
        do j=2,n-1
            tmp(1) = mat%du(j-1)
            tmp(2) = mat%d(j)
            tmp(3) = mat%dl(j)
            call zlassq(3, tmp(1), 1, scale, sum)
        end do
        ! column n
        tmp(1) = mat%du(n-1)
        tmp(2) = mat%d(n)
        call zlassq(2, tmp(1), 1, scale, sum)
        ! store result
        res = scale*sqrt(sum)
        return
    end function FroNorm
end module TMP_EigSolver3