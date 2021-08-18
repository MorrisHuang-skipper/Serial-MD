! --------------------------------------------------------------------
! MODULE  math

!    Subroutines:
!    (1) mag2              - calulate |A|^2
!
!    Functions:
!    (1) 
! TODO:
! --------------------------------------------------------------------
module math
    implicit none
    double precision, parameter :: sqrtpi = dsqrt(4.d0*atan(1.d0))
    double precision, parameter :: pi = 4.d0*atan(1.d0)
    double precision, parameter :: twopi = 2.d0 * pi
    double precision, parameter :: sqrt2 = dsqrt(2.d0)
    ! some constant
    double precision, parameter :: c0=2.99792458d8, me=9.10938d-31, mi=1836.d0*me
    double precision, parameter :: e=1.602d-19, k=8.9875d9
    double precision, parameter :: mu0=1.256637d-6, eps0=8.854187d-12
    double precision, parameter :: kb=1.380649e-23, NA=6.02214d23
    double precision, parameter :: Da=1.660539066d-27
    double precision, parameter :: mu=1d-6, nano=1.d-9, ang=1.d-10, pico=1.d-12
    double precision, parameter :: femto=1.d-15, atto=1.d-18
    double precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0
    double precision, parameter :: five=5.d0, six=6.d0, seven=7.d0, eight=8.d0, nine=9.d0
    double precision, parameter :: ten=10.d0, half=.5d0
    logical, parameter :: yes=.True., no=.False.
    logical, parameter :: on=.True., off=.False.
    contains

! --------------------------------------------------------------------
!    Purpose: |A|^2
! --------------------------------------------------------------------
    function mag2(A)
        implicit none
        double precision :: mag2, A(3) 
        mag2 = dot_product(A, A)
    end function


! --------------------------------------------------------------------
!    Purpose: |A|
! --------------------------------------------------------------------
    function mag(A)
        implicit none
        double precision :: mag, A(3) 
        mag = dsqrt(dot_product(A, A))
    end function


! --------------------------------------------------------------------
!    Purpose: unit vector
! --------------------------------------------------------------------
    function hat(A)
        implicit none
        double precision :: hat(3), A(3)
        hat = A / mag(A)
    end function


! --------------------------------------------------------------------
!    Purpose: randomly distribute in a sphere
!    
!    input:
!          (1) r: radius of sphere
!          (2) cx, cy, cz: center of sphere in x, y, z
!          (3) dim: dimension of sphere
!    output:
!          (1) (x, y, z)
! --------------------------------------------------------------------
    function RDsphere(r, cx, cy, cz, dim)
        implicit none
        double precision :: RDsphere(3), r, a, b, c
        double precision :: cx, cy, cz
        integer :: dim
        !> 1d sphere
        if (dim == 1) then
            call random_number(a)
            a = (a - .5d0) * 2.d0 * r
            RDsphere = (/ a+cx, cy, cz /)
        end if
        !> 2d sphere
        if (dim == 2) then
            do
                call random_number(a)
                call random_number(b)
                a = (a - .5d0) * 2.d0 * r
                b = (b - .5d0) * 2.d0 * r
                if ((a**2+b**2) .le. r**2) then
                    exit
                end if
            end do
            RDsphere = (/ a+cx, b+cy, cz /)
        end if
        !> 3d sphere
        if (dim == 3) then
            do
                call random_number(a)
                call random_number(b)
                call random_number(c)
                a = (a - .5d0) * 2.d0 * r
                b = (b - .5d0) * 2.d0 * r
                c = (c - .5d0) * 2.d0 * r
                if ((a**2+b**2+c**2) .le. r**2) then
                    exit
                end if
            end do
            RDsphere = (/ a+cx, b+cy, c+cz /)
        end if
    end function RDsphere


! --------------------------------------------------------------------
!    Purpose: randomly distribute in a rectangle
!    
!    input:
!          (1) Lx, Ly, Lz: size in x, y, z direction
!          (2) cx, cy, cz: center of rectangle
!          (3) dim: dimension
!    output:
!          (1) (x, y, z) of random coordinate
! --------------------------------------------------------------------
    function RDrect(Lx, Ly, Lz, cx, cy, cz, dim)
        implicit none
        integer :: dim
        double precision :: RDrect(3), Lx, Ly, Lz, x, y, z
        double precision :: cx, cy, cz
        call random_seed()
        if (dim == 1) then
            call random_number(x)
            x = (x - .5d0) * Lx
            RDrect = (/ x+cx, cy, cz /)
        end if
        if (dim == 2) then
            call random_number(x)
            call random_number(y)
            x = (x - .5d0) * Lx
            y = (y - .5d0) * Ly
            RDrect = (/ x+cx, y+cy, cz /)
        end if
        if (dim == 3) then
            call random_number(x)
            call random_number(y)
            call random_number(z)
            x = (x - .5d0) * Lx
            y = (y - .5d0) * Ly
            z = (z - .5d0) * Lz
            RDrect = (/ x+cx, y+cy, z+cz /)
        end if
    end function RDrect

! --------------------------------------------------------------------
!    Purpose: create gaussian random distribution
!    
!    input:
!          (1) mu: mean
!          (2) sig: sigma, standard deviation
!    output:
!          (1) gaurand: one gaussian distribtution number
! --------------------------------------------------------------------
    function gaurand(mu, sig)
        implicit none
        double precision :: a, b, gaurand, sig, mu
        call random_number(a)
        call random_number(b)
        gaurand = sig * dsqrt(-2.d0 * log(a)) * dcos(twopi * b) + mu
    end function

! --------------------------------------------------------------------
!    Purpose: create maxwell distribution for specific dimension
!    
!    input:
!          (1) mass: mass of species
!          (2) temperature: temperature of species
!          (3) dim: dimension
!    output:
!          (1) maxwell: (/ a, b, c /)
! --------------------------------------------------------------------
    function maxwell(temperature, mass, dim)
        implicit none
        integer :: dim, i
        double precision :: maxwell(3), mass, temperature
        maxwell = zero
        select case (dim)
        case(1)
            do i = 1, 1
                maxwell(i) = dsqrt(kb * temperature / mass) * gaurand(zero, one)
            end do
        case(2)
            do i = 1, 2
                maxwell(i) = dsqrt(kb * temperature / mass) * gaurand(zero, one)
            end do
        case(3)
            do i = 1, 3
                maxwell(i) = dsqrt(kb * temperature / mass) * gaurand(zero, one)
            end do
        end select
    end function maxwell


end module