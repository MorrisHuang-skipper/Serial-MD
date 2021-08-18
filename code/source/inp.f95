! --------------------------------------------------------------------
!    Purpose: specified simulaton input parameter
! --------------------------------------------------------------------
module inp
    use math
    implicit none

    ! --------------------------------------------------------------------------------
    !> number of total particles (electron + hydrogen)
    integer :: dim
    !> target radius and density (Tin)
    double precision :: radius
    !> number of particles of different species and total particles
    integer :: NP1, NP2, NP
    !> simulation box length
    double precision :: Lx, Ly, Lz
    !> output data directory
    character(len=30) :: loc
    !> timestep
    double precision :: dt, maxt
    !> charge size
    double precision :: w0, wpic
    !> cutoff radius
    double precision :: rcut
    !> output data duration
    double precision :: dumptime
    !> initial target temperature
    double precision :: temperature
    !> laser parameter
    logical :: E_ext
    double precision :: freq, omega
    double precision :: lambda, waveN
    double precision :: Intensity, Emax
    double precision :: laser_start, laser_end
    !> periodic boudary
    logical :: periodic
    !> initial distribution (box, sphere)
    character(len=30) :: INIT_DIS
    !> initial velocity distribution (gaussian, uniform)
    character(len=30) :: INIT_v


    !> output data timestep
    integer :: dumpstep
    integer :: Numdump
    !> simulation time
    integer :: Nstep
    !> number of auxiliary cell in () direction
    integer :: Ncellx, Ncelly, Ncellz, totCell
    !> auxiliary cell length (>=rcut)
    double precision :: cLx, cLy, cLz
    !> simulation time
    double precision :: time=0.d0
    !> output data timestep
    integer :: idump=0
    !> particles position and velocity
    double precision, allocatable, dimension(:, :) :: r, v, vb
    !> heads of the cell list
    integer, allocatable, dimension(:) :: head, list
    ! --------------------------------------------------------------------------------

    ! --------------------------------------------------------------------------------
    namelist /input/ dim, radius, NP1, NP2, &
                     Lx, Ly, Lz, loc, dt, maxt, w0,  &
                     wpic, rcut, NumDump, temperature, &
                     periodic, INIT_DIS, INIT_v

    namelist /laser/ E_ext, freq, omega, lambda, waveN, &
                     Intensity, Emax, laser_start, laser_end

    namelist /derived/ NP, Nstep, dumpstep, dumptime, &
                       Ncellx, Ncelly, Ncellz, totCell, &
                       cLx, cLy, cLz
    ! --------------------------------------------------------------------------------
    ! !> number of total particles (electron + hydrogen)
    ! integer, parameter :: dimension=3
    ! !> target radius and density (Tin)
    ! double precision, parameter :: radius=5.d0*nano
    ! !> number of particles of different species
    ! integer, parameter :: NP1=1000, NP2=1000
    ! !> total number of paricles
    ! integer, parameter :: NP = NP1+NP2
    ! !> simulation box length
    ! double precision, parameter :: Lx=40*nano, Ly=40*nano, Lz=40*nano
    ! !> output data directory
    ! character(len=30), parameter :: loc = '../../data/ext/' 
    ! !> timestep
    ! double precision, parameter :: dt= 1.d0*atto, maxt=100.d0*femto 
    ! !> charge size
    ! double precision, parameter :: w0=1.d0*ang, wpic=5.d0*ang
    ! !> cutoff radius
    ! double precision, parameter :: rcut = 3.d0*wpic
    ! !> output data duration
    ! integer, parameter :: dmp=100
    ! double precision, parameter :: dmptime=maxt/dble(dmp)
    ! !> initial target temperature
    ! double precision, parameter :: temperature=zero*e
    ! !> laser parameter
    ! logical, parameter :: E_ext = off
    ! double precision, parameter :: freq=one/femto, omega=twopi*freq
    ! double precision, parameter :: lambda=c0/freq, waveN=twopi/lambda
    ! double precision, parameter :: Intensity=0.d16, Emax=dsqrt(two*Intensity/c0/eps0)
    ! double precision, parameter :: laser_start=50.d0*femto, laser_end=100.d0*femto
    ! !> periodic boudary
    ! logical, parameter :: periodic=off
    ! !> initial distribution
    ! ! character(len=30), parameter :: INIT_DIS = 'box'
    ! character(len=30), parameter :: INIT_DIS = 'sphere'
    ! !> initial velocity distribution
    ! ! character(len=30), parameter :: INIT_v = 'gaussian'
    ! character(len=30), parameter :: INIT_v = 'uniform'
    ! ! --------------------------------------------------------------------------------

    contains
    

        ! subroutine read_input(path)
        !     implicit none
        !     character(len=*), intent(in) :: path
        !     open(66, file=path)
        !     write(66, nml=input)
        ! end subroutine read_input

! --------------------------------------------------------------------
!    Purpose: species charge
! --------------------------------------------------------------------
        function charge(idx)
            implicit none
            double precision :: charge
            integer :: idx

            if (idx .le. NP1) then
                charge = -e
            else if ((idx .ge. NP1+1).and.(idx .le. NP1+NP2)) then
                charge = e
            else
                charge = zero
                print *, 'invalid index'
            end if
        end function charge


! --------------------------------------------------------------------
!    Purpose: species mass
! --------------------------------------------------------------------
        function mass(idx)
            implicit none
            double precision :: mass
            integer :: idx

            if (idx .le. NP1) then
                mass = me
            else if ((idx .ge. NP1+1).and.(idx .le. NP1+NP2)) then
                mass = mi
            else
                mass = zero
                print *, 'invalid index'
            end if
        end function mass

end module inp