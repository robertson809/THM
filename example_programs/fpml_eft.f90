!********************************************************************************
!   FPML_EFT: Error-Free Transformations for FPML
!   Author: Thomas R. Cameron
!   Institution: Davidson College, Mathematics and Computer Science Department
!   Last Modified: 13 December 2018
!********************************************************************************
module fpml_eft
    implicit none
    integer, parameter                          :: dp = kind(1.d0), factor = 2**27 + 1
    !********************************************************
    !   REFT: Real Error Free Transformation type.
    !********************************************************
    !   Made up of two real(kind=dp) variables x and y
    !   called primary and secondary, respectively. 
    !   For TwoSum and TwoProd, x will denote floating-
    !   point result and y will denote error in result.
    !   For Split, x and y are the two parts of the
    !   floating-point number. 
    !********************************************************                 
    type REFT
        real(kind=dp)                           :: x        ! primary
        real(kind=dp)                           :: y        ! secondary
    end type REFT
    !********************************************************
    !   CEFTSum: Complex Error Free Transformation Sum.
    !********************************************************
    !   Made up of two complex(kind=dp) variables x and y
    !   called primary and secondary, respectively. 
    !   For TwoSumCplx, x denotes the floating-point 
    !   result and y denotes the error.
    !********************************************************
    type CEFTSum
        complex(kind=dp)                        :: x        ! primary
        complex(kind=dp)                        :: y        ! secondary
    end type CEFTSum
    !********************************************************
    !   CEFTProd: Complex Error Free Transformation Product.
    !********************************************************
    !   Made up of four complex(kind=dp) variables p, e, f,
    !   and g, which we call primary, secondary, tertiary,
    !   and quaternary, respectively. For TwoProductCplx,
    !   p denotes the floating-point result, and e, f, and g
    !   denote errors associated with the sum and product
    !   of real and imaginary parts.                        
    !********************************************************
    type CEFTProd
        complex(kind=dp)                        :: p        ! primary
        complex(kind=dp)                        :: e        ! secondary
        complex(kind=dp)                        :: f        ! tertiary
        complex(kind=dp)                        :: g        ! quaternary
    end type CEFTProd
    !********************************************************
    !   CEFTHorner: Complex Error Free Transformation Horner.
    !********************************************************
    !   Made up of a complex(kind=dp) variable h to store 
    !   result of standard Horner method and four allocatable
    !   complex(kind=dp) variables to store polynomial 
    !   coefficients p, q, r, and s, where pi, qi, ri is the
    !   error in the product term and si is the error in the
    !   sum term on the ith iteration of Horner's method.   
    !********************************************************
    type CEFTHorner
        complex(kind=dp)                        :: h        ! result
        complex(kind=dp), allocatable           :: p(:)     ! error in product
        complex(kind=dp), allocatable           :: q(:)     ! error in product
        complex(kind=dp), allocatable           :: r(:)     ! error in product
        complex(kind=dp), allocatable           :: s(:)     ! error in sum
    end type CEFTHorner
contains
    !********************************************************
    !                       TwoSum                          *
    !********************************************************
    ! Computes the sum of two floating-point numbers
    ! a and b. The result and error is returned in
    ! comp which has type REFT. 
    !********************************************************
    function TwoSum(a,b) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)       :: a, b
        ! local variables
        real(kind=dp)       :: z
        type(REFT)          :: comp

        ! compute floating-point sum
        comp%x = a + b
        z = comp%x - a
        ! compute error in floating-point sum
        comp%y = (a - (comp%x - z)) + (b - z)
        return
    end function TwoSum
    !****************************************************
    !                       AccSum                      *
    !****************************************************
    ! Computes the accurate sum of the floating-point 
    ! numbers a, b, c, and d. Result is returned in comp 
    ! which has type real(kind=dp).          
    !****************************************************
    function AccSum(a,b,c,d) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)           :: a, b, c, d
        ! local variables
        real(kind=dp)           :: comp, cor
        type(REFT)              :: sum
            
        ! compute eft sum of a and b
        sum = TwoSum(a,b)
        ! store result and error
        comp = sum%x
        cor = sum%y
        ! compute eft sum of comp and c
        sum = TwoSum(comp,c)
        ! update result and error
        comp = sum%x
        cor = cor + sum%y
        ! compute eft sum of comp and d
        sum = TwoSum(comp,d)
        ! update result and error
        comp = sum%x
        cor = cor + sum%y
        ! add error back into result
        comp = comp + cor
        return
    end function AccSum
    !********************************************************
    !                       TwoSumCplx                      *
    !********************************************************
    ! Computes the sum of two complex floating-point 
    ! numbers a and b. The result and error is 
    ! returned in comp which has type CEFTSum. 
    !********************************************************
    function TwoSumCplx(a,b) result(comp)
        implicit none
        ! argument variables
        complex(kind=dp)        :: a, b
        ! local variables
        type(REFT)              :: rcomp, icomp
        type(CEFTSum)           :: comp
            
        ! compute sum of real and imaginary parts
        rcomp = TwoSum(real(a),real(b))
        icomp = TwoSum(aimag(a),aimag(b))
        ! store primary and secondary complex numbers in comp
        comp%x = cmplx(rcomp%x,icomp%x,kind=dp)
        comp%y = cmplx(rcomp%y,icomp%y,kind=dp)
        return
    end function TwoSumCplx
    !****************************************************
    !                       AccSumCplx                  *
    !****************************************************
    ! Computes the accurate sum of the complex floating-
    ! point numbers a, b, c, and d. Result is returned
    ! in comp which has type complex(kind=dp).          
    !****************************************************
    function AccSumCplx(a,b,c,d) result(comp)
        implicit none
        ! argument variables
        complex(kind=dp)        :: a, b, c, d
        ! local variables
        complex(kind=dp)        :: comp, cor
        type(CEFTSum)           :: sum
            
        ! compute eft sum of a and b
        sum = TwoSumCplx(a,b)
        ! store result and error
        comp = sum%x
        cor = sum%y
        ! compute eft sum of comp and c
        sum = TwoSumCplx(comp,c)
        ! update result and error
        comp = sum%x
        cor = cor + sum%y
        ! compute eft sum of comp and d
        sum = TwoSumCplx(comp,d)
        ! update result and error
        comp = sum%x
        cor = cor + sum%y
        ! add error back into result
        comp = comp + cor
        return
    end function AccSumCplx
    !********************************************************
    !                       Split                           *
    !********************************************************
    ! Computes the splitting of a floating-point 
    ! number a. Both parts have at most 26 nonzero bits 
    ! and are returned in comp which has type REFT.
    !********************************************************
    function Split(a) result(comp)
        implicit none
        ! argument variables
        real(kind=dp)       :: a
        ! local variables   
        real(kind=dp)       :: c
        type(REFT)          :: comp
        
        ! factor a
        c = factor * a
        ! split into two 26-bit numbers
        comp%x = c - (c - a)
        comp%y = a - comp%x
        return
    end function Split
    !********************************************************
    !                       TwoProductCplx                  *
    !********************************************************
    ! Computes the product of two complex floating-point 
    ! numbers a and b. The result and error is returned 
    ! in comp which has type CEFTP. 
    !********************************************************
    function TwoProductCplx(a,b) result(comp)
        implicit none
        ! argument variables
        complex(kind=dp)    :: a, b
        ! local variables
        real(kind=dp)       :: z1, z2, z3, z4
        type(REFT)          :: sra, sia, srb, sib, sumr, sumi
        type(CEFTProd)      :: comp
            
        ! perform splitting of real(a), imag(a), real(b), imag(b)
        sra = Split(real(a)); sia = Split(aimag(a))
        srb = Split(real(b)); sib = Split(aimag(b))
        ! compute product of real and imaginary parts
        z1 = real(a) * real(b); z2 = aimag(a) * aimag(b)
        z3 = real(a) * aimag(b); z4 = aimag(a) * real(b)
        ! error in z1 and z3
        comp%e = cmplx(sra%y * srb%y - (((z1 - sra%x * srb%x) - sra%y * srb%x) - sra%x * srb%y),&
                    sra%y * sib%y - (((z3 - sra%x * sib%x) - sra%y * sib%x) - sra%x * sib%y),kind=dp)
        ! error in z2 and z4
        comp%f = cmplx(-(sia%y * sib%y - (((z2 - sia%x * sib%x) - sia%y * sib%x) - sia%x * sib%y)),&
                    sia%y * srb%y - (((z4 - sia%x * srb%x) - sia%y * srb%x) - sia%x * srb%y),kind=dp)
        ! compute sum z1-z2
        sumr = TwoSum(z1,-z2)
        ! compute sum z3+z4
        sumi = TwoSum(z3,z4)
        ! store floating-point result
        comp%p = cmplx(sumr%x,sumi%x,kind=dp)
        ! error in sums
        comp%g = cmplx(sumr%y,sumi%y,kind=dp)
        return
    end function TwoProductCplx
end module fpml_eft