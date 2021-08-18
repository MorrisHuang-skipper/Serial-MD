module species
    use math
    implicit none

    contains

! ! --------------------------------------------------------------------
! !    Purpose: species particle number for sperical droplet
! ! --------------------------------------------------------------------
!         function SpeNum(radius, sp)
!             implicit none
!             !> output number of particle of specific species within sphere
!             integer :: SpeNum
!             !> radius of target (SI unit)
!             double precision :: radius
!             !> species name
!             character(len=30) :: sp
!             !> atomic mass
!             double precision :: AM
!             !> density (g/cm3)
!             double precision :: density

!             if (sp .eq. 'Tin') then
!                 AM = 118.71d0
!                 density = 7.287d0
!                 SpeNum = int(four / three * pi * radius * radius * radius * &
!                              density * 1.d3 / (AM * DA))
!             end if

!             if (sp .eq. 'H') then
!                 AM = 1.008d0
!                 density = 0.899d0 / 1.d3
!                 SpeNum = int(four / three * pi * radius * radius * radius * &
!                              density * 1.d3 / (AM * DA))
!             end if

!         end function SpeNum


! ! --------------------------------------------------------------------
! !    Purpose: species particle atomic number
! ! --------------------------------------------------------------------
!         function SpeZ(sp)
!             implicit none
!             character(len=30) :: sp
!             integer :: SpeZ

!             if (sp .eq. 'Tin') then
!                 SpeZ = 50
!             else if (sp .eq. 'H') then
!                 SpeZ = 1
!             else
!                 print *, 'Invalid species'
!             end if
!         end function SpeZ
end module species