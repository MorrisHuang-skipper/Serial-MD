! --------------------------------------------------------------------
! MODULE  Mic
!
!    Modules:
!    (1) var                  - contain all simulation variables
!
!    Subroutines:
!    (1) density              - calculate electron density
!
!    Functions:
!    (1) 
! TODO:
!    (1) openmp test
!    (2) applying external field
!    (3) periodic bd condition
!    (4) Boris
!    (5) EM PIC 3D
! --------------------------------------------------------------------
module mic
    use omp_lib
    use math
    use inp
    implicit none

    contains

! --------------------------------------------------------------------
!    Purpose: calculate relativistic speed from energy (ev)
! --------------------------------------------------------------------
    function calv(energy, m)
        implicit none
        double precision :: calv, energy, m
        ! calv = c0 * dsqrt(1.d0-1.d0/(1.d0+ev*e/mass/c0**2)**2)
        ! non-relativistic
        calv = dsqrt(two * energy / m)
    end function


! --------------------------------------------------------------------
!    Purpose: initialize position and velocity
! --------------------------------------------------------------------
    subroutine INIT
        implicit none
        integer :: i
        double precision :: vel_e, vel_i 
 
        call random_seed()

        vel_e = calv(temperature, me)
        vel_i = calv(temperature, mi)

        if (INIT_v .eq. 'uniform') then
            !> species one (eg. electron)
            do i = 1, NP1
                v(i, :) = RDsphere(r=vel_e, cx=zero, cy=zero, cz=zero, dim=dim)
            end do
            !> species two (eg. hydorgen)
            do i = NP1+1, NP1+NP2
                v(i, :) = RDsphere(r=vel_i, cx=zero, cy=zero, cz=zero, dim=dim)
            end do
        else if (INIT_v .eq. 'gaussian') then
            do i = 1, NP1
                v(i, :) = maxwell(temperature=temperature, mass=mass(i), dim=dim)
            end do
            do i = NP1+1, NP1+NP2
                v(i, :) = maxwell(temperature=temperature, mass=mass(i), dim=dim)
            end do
        end if

        !> spherical unifrom distribution for species 1 and species 2
        if (INIT_DIS .eq. 'sphere') then
            do i = 1, NP
                r(i, :) = RDsphere(r=radius, cx=Lx/two, cy=Ly/two, cz=Lz/two, dim=3)
            end do
        else if (INIT_DIS .eq. 'box') then
        !> random rectangle distribution
            select case (dim)
            case (1)
                do i = 1, NP1
                    r(i, :) = RDrect(Lx=Lx, Ly=zero, Lz=zero, cx=Lx/two, cy=Ly/two, &
                                    cz=Lz/two, dim=dim)
                end do
                do i = NP1+1, NP1+NP2
                    r(i, :) = RDrect(Lx=Lx, Ly=zero, Lz=zero, cx=Lx/two, cy=Ly/two, &
                                    cz=Lz/two, dim=dim)
                end do
            case (2)
                do i = 1, NP1
                    r(i, :) = RDrect(Lx=Lx, Ly=Ly, Lz=zero, cx=Lx/two, cy=Ly/two, &
                                    cz=Lz/two, dim=dim)
                end do
                do i = NP1+1, NP1+NP2
                    r(i, :) = RDrect(Lx=Lx, Ly=Ly, Lz=zero, cx=Lx/two, cy=Ly/two, &
                                    cz=Lz/two, dim=dim)
                end do
            case (3)
                do i = 1, NP1
                    r(i, :) = RDrect(Lx=Lx, Ly=Ly, Lz=Lz, cx=Lx/two, cy=Ly/two, &
                                    cz=Lz/two, dim=dim)
                end do
                do i = NP1+1, NP1+NP2
                    r(i, :) = RDrect(Lx=Lx, Ly=Ly, Lz=Lz, cx=Lx/two, cy=Ly/two, &
                                    cz=Lz/two, dim=dim)
                end do
            end select
        else
            print *, 'did not assign INIT distribution'
        end if

    end subroutine INIT


! --------------------------------------------------------------------
!    Purpose: leapfrog method first leap -1/2 step
! --------------------------------------------------------------------
    subroutine halfleap
        implicit none
        v = v - accel() * dt / 2.d0
    end subroutine halfleap


! --------------------------------------------------------------------
!    Purpose: particle mover
! --------------------------------------------------------------------
    subroutine move(method)
        implicit none
        character(len=2) :: method
        !> particle mover: leapfrog
        if (method .eq. 'LF') then
            !> v(n-1/2)
            vb = v
            !> v(n+1/2)
            v = v + accel() * dt
            r = r + v * dt
        end if
        !> apply periodic boundary conditin
        if (periodic .eqv. yes) then
            if (dim == 1) then
                call periodx
            end if

            if (dim == 2) then
                call periodx
                call periody
            end if

            if (dim == 3) then
                call periodx
                call periody
                call periodz
            end if
        end if
        time = time + dt
    end subroutine move


! --------------------------------------------------------------------
!    Purpose: apply periodic boundary condition for particles
! --------------------------------------------------------------------
    subroutine periodx
        implicit none
        integer :: i
        do i = 1, NP
            if (r(i, 1) .gt. Lx) then
                r(i, 1) = r(i, 1) - Lx
            else if (r(i, 1) .lt. zero) then
                r(i, 1) = r(i, 1) + Lx
            end if
        end do
    end subroutine periodx

    subroutine periody
        implicit none
        integer :: i
        do i = 1, NP
            if (r(i, 2) .gt. Ly) then
                r(i, 2) = r(i, 2) - Ly
            else if (r(i, 2) .lt. zero) then
                r(i, 2) = r(i, 2) + Ly
            end if
        end do
    end subroutine periody

    subroutine periodz
        implicit none
        integer :: i
        do i = 1, NP
            if (r(i, 3) .gt. Lz) then
                r(i, 3) = r(i, 3) - Lz
            else if (r(i, 3) .lt. zero) then
                r(i, 3) = r(i, 3) + Lz
            end if
        end do
    end subroutine periodz


! --------------------------------------------------------------------
!    Purpose: Create a cell link
! --------------------------------------------------------------------
    subroutine links
        implicit none
        integer :: i, xcell, ycell, zcell, icell
        head = 0
        list = 0
        do i = 1, NP
            xcell = ceiling(r(i, 1)/cLx)
            ycell = ceiling(r(i, 2)/cLy)
            zcell = ceiling(r(i, 3)/cLz)
            !> if no periodic, consider particle only in the simulation box
            if (periodic .eqv. no) then
                if ((xcell .lt. 1).or.(xcell .gt. Ncellx)) cycle
                if ((ycell .lt. 1).or.(ycell .gt. Ncelly)) cycle
                if ((zcell .lt. 1).or.(zcell .gt. Ncellz)) cycle
            end if
            icell = (xcell-1) * Ncelly * Ncellz + (ycell-1) * Ncellz + zcell
            list(i) = head(icell)
            head(icell) = i
        end do
    end subroutine links


! --------------------------------------------------------------------
!    purpose: calculate force from particles within rcut 
! --------------------------------------------------------------------
    function brute_force(i)
        implicit none
        integer i, j
        double precision :: brute_force(3), D

        brute_force = zero
        do j = 1, NP
            ! calculate force within the cutoff radius
            if (i .ne. j) then
                D = mag(r(i, :)-r(j, :))
                if (D .le. rcut) then
                    brute_force = brute_force + &
                    Fmic(r(i, :), r(j, :), charge(i), charge(j))
                end if
            end if
        end do
        if (E_ext .eqv. on) then
            brute_force = brute_force + charge(i) * external_field(i)
        end if
    end function
 

! --------------------------------------------------------------------
!    purpose: calculate force from particles within rcut 
! --------------------------------------------------------------------
    function force(i)
        implicit none
        integer :: i, j, xcell, ycell, zcell, icell, x, y, z
        double precision :: force(3), D

        xcell = 0
        ycell = 0
        zcell = 0
        force = zero
        !> find the cell particle i lives in
        x = ceiling(r(i, 1)/cLx)
        y = ceiling(r(i, 2)/cLy)
        z = ceiling(r(i, 3)/cLz)
        !> 27 nearest cells to checks
        do xcell = x-1, x+1, 1
            do ycell = y-1, y+1, 1
                do zcell = z-1, z+1, 1
                    if ((xcell .lt. 1).or.(xcell .gt. Ncellx)) cycle
                    if ((ycell .lt. 1).or.(ycell .gt. Ncelly)) cycle
                    if ((zcell .lt. 1).or.(zcell .gt. Ncellz)) cycle
                    icell = (xcell-1) * Ncelly * Ncellz + (ycell-1) * Ncellz + zcell
                    j = head(icell)
                    do 
                        if (j .eq. 0) exit
                        D = mag(r(i, :) - r(j, :))
                        if ((D .le. rcut).and.(i .ne. j)) then
                            force = force + Fmic(r(i, :), r(j, :), charge(i), charge(j))
                        end if
                        j = list(j)
                    end do
                end do
            end do
        end do
        if (E_ext .eqv. on) then
            !> add external force
            force = force + charge(i) * external_field(i)
        end if
    end function force


! --------------------------------------------------------------------
!    purpose: calculate force from nearest image
! --------------------------------------------------------------------
    function PDforce(i)
        implicit none
        integer :: i, j, xcell, ycell, zcell, icell, x, y, z
        integer :: xc, yc, zc
        double precision :: PDforce(3), D, dx, dy, dz
        logical :: xbd, ybd, zbd

        xcell = 0
        ycell = 0
        zcell = 0
        PDforce = zero
        !> find the cell particle i lives in
        x = ceiling(r(i, 1)/cLx)
        y = ceiling(r(i, 2)/cLy)
        z = ceiling(r(i, 3)/cLz)
        !> 27 nearest cells to checks
        do xcell = x-1, x+1, 1
            do ycell = y-1, y+1, 1
                do zcell = z-1, z+1, 1
                    !> x boundary
                    if (xcell .eq. 0) then
                        xc = Ncellx
                        xbd = yes
                    else if (xcell .eq. Ncellx+1) then
                        xc = 1
                        xbd = yes
                    else 
                        xc = xcell
                        xbd = no
                    end if
                    !> y boundary
                    if (ycell .eq. 0) then
                        yc = Ncelly
                        ybd = yes
                    else if (ycell .eq. Ncelly+1) then
                        yc = 1
                        ybd = yes
                    else
                        yc = ycell
                        ybd = no
                    end if
                    !> z boundary
                    if (zcell .eq. 0) then
                        zc = Ncellz
                        zbd = yes
                    else if (zcell .eq. Ncellz+1) then
                        zc = 1
                        zbd = yes
                    else
                        zc = zcell
                        zbd = no
                    end if
                    !> normal case
                    icell = (xc-1) * Ncelly * Ncellz + (yc-1) * Ncellz + zc
                    j = head(icell)
                    do 
                        if (j .eq. 0) exit
                        if (xbd .eqv. yes) then
                            dx = Lx - abs(r(j, 1) - r(i, 1)) 
                        else
                            dx = r(j, 1) - r(i, 1)
                        end if
                        if (ybd .eqv. yes) then
                            dy = Ly - abs(r(j, 2) - r(i, 2))
                        else
                            dy = r(j, 2) - r(i, 2)
                        end if
                        if (zbd .eqv. yes) then
                            dz = Lz - abs(r(j, 3) - r(i, 3))
                        else
                            dz = r(j, 3) - r(i, 3)
                        end if
                        D = dsqrt(dx * dx + dy * dy + dz * dz)
                        if ((D .le. rcut).and.(i .ne. j)) then
                            PDforce = PDforce + Fmic(r(i, :), r(j, :), charge(i), charge(j))
                        end if
                        j = list(j)
                    end do
                end do
            end do
        end do
        !> add external force
        if (E_ext .eqv. on) then
            PDforce = PDforce + charge(i) * external_field(i)
        end if
    end function PDforce


! --------------------------------------------------------------------
!    purpose: calculate external laser field of particle i
! --------------------------------------------------------------------
    function external_field(i)
        implicit none
        integer :: i
        double precision :: external_field(3), Ey
        !> eg. sinusoidal plane wave (k in y axis)
        if ((time .ge. laser_start).and.(time .le. laser_end)) then
            Ey = Emax*dsin(waveN * r(i, 1) - omega * time)
            external_field = (/ zero, Ey, zero /)
        else
            external_field = (/ zero, zero, zero /)
        end if
    end function external_field


! --------------------------------------------------------------------
!    Purpose: calculate mic force
! --------------------------------------------------------------------
    function Fmic(ri, rj, qi, qj)
        implicit none 
        double precision :: Fmic(3), x, ri(3), rj(3), qi, qj
        x = mag(ri - rj)
        Fmic = (-sqrt2 * dexp(-x * x / two / w0 / w0) / &
               (sqrtpi * w0 * x) + erf(x / sqrt2 / w0) / x / x) &
               * hat(ri - rj) * qi * qj * k
    end function Fmic


! --------------------------------------------------------------------
!    Purpose: calculate mic potential between particle i and j
! --------------------------------------------------------------------
    function Umic(x, i, j)
        implicit none
        integer :: i, j
        double precision :: Umic, x
        Umic = erf(x / sqrt2 / w0) * k * charge(i) * charge(j) / x
    end function


! --------------------------------------------------------------------
!    Purpose: calculate accelaration of particles
! --------------------------------------------------------------------
    function accel()
        implicit none
        double precision :: accel(NP, 3)
        integer :: i
        accel = zero
        !$omp parallel
        !$omp do
        do i = 1, NP
            !> using linked list to find cutoff force
            if (periodic .eqv. yes) then
                !> calcuulate force for periodic bd condition
                accel(i, :) = PDforce(i) / mass(i)
            else
                !> calculate force within the simulation box
                accel(i, :) = force(i) / mass(i)
            end if
            !> using brute comparison to find cutoff force
            ! accel(i, :) = brute_force(i) / mass(i)
        end do
        !$omp end do
        !$omp end parallel
    end function


! --------------------------------------------------------------------
!    Purpose: calculate kinetic energy (move v back dt/2 step)
! --------------------------------------------------------------------
    function KE(start, finish)
        implicit none
        double precision :: KE
        integer :: i, start, finish
        KE = zero
        do i = start, finish
            KE = KE + mass(i) * mag2(half * (vb(i, :) + v(i, :)))
        end do
        KE = half * KE
    end function


! --------------------------------------------------------------------
!    Purpose: calculate total potential energy
! --------------------------------------------------------------------
    function Potential()
        implicit none
        double precision :: Potential, D
        integer :: i, j, x, y, z, xcell, ycell, zcell, icell
        Potential = zero
        xcell = 0
        ycell = 0
        zcell = 0
        icell = 0
        do i = 1, NP
            !> find the cell particle i in
            x = ceiling(r(i, 1)/cLx)
            y = ceiling(r(i, 2)/cLy)
            z = ceiling(r(i, 3)/cLz)
            !> 27 nearest cells to checks
            do xcell = x-1, x+1, 1
                do ycell = y-1, y+1, 1
                    do zcell = z-1, z+1, 1
                        if ((xcell .lt. 1).or.(xcell .gt. Ncellx)) cycle
                        if ((ycell .lt. 1).or.(ycell .gt. Ncelly)) cycle
                        if ((zcell .lt. 1).or.(zcell .gt. Ncellz)) cycle
                        icell = (xcell-1) * Ncelly * Ncellz + (ycell-1) * Ncellz + zcell
                        j = head(icell)
                        do 
                            if (j .eq. 0) exit
                            D = mag(r(i, :)-r(j, :))
                            if ((D .le. rcut).and.(i .ne. j)) then
                                Potential = Potential + Umic(D, i, j)
                            end if
                            j = list(j)
                        end do
                    end do
                end do
            end do
        end do
        Potential = half * Potential
    end function


! --------------------------------------------------------------------
!    Purpose: calculate total potential energy
! --------------------------------------------------------------------
    function PDpotential()
        implicit none
        double precision :: PDpotential, D
        integer :: i, j, x, y, z, xcell, ycell, zcell, icell
        integer :: xc, yc, zc
        double precision :: dx, dy, dz
        logical :: xbd, ybd, zbd
        ! type(container), pointer :: tmp
        PDpotential = zero
        xcell = 0
        ycell = 0
        zcell = 0
        icell = 0
        do i = 1, NP
            !> find the cell particle i in
            x = ceiling(r(i, 1)/cLx)
            y = ceiling(r(i, 2)/cLy)
            z = ceiling(r(i, 3)/cLz)
            !> 27 nearest cells to checks
            do xcell = x-1, x+1, 1
                do ycell = y-1, y+1, 1
                    do zcell = z-1, z+1, 1
                        !> x boundary
                        if (xcell .eq. 0) then
                            xc = Ncellx
                            xbd = yes
                        else if (xcell .eq. Ncellx+1) then
                            xc = 1
                            xbd = yes
                        else 
                            xc = xcell
                            xbd = no
                        end if
                        !> y boundary
                        if (ycell .eq. 0) then
                            yc = Ncelly
                            ybd = yes
                        else if (ycell .eq. Ncelly+1) then
                            yc = 1
                            ybd = yes
                        else
                            yc = ycell
                            ybd = no
                        end if
                        !> z boundary
                        if (zcell .eq. 0) then
                            zc = Ncellz
                            zbd = yes
                        else if (zcell .eq. Ncellz+1) then
                            zc = 1
                            zbd = yes
                        else
                            zc = zcell
                            zbd = no
                        end if
                        !> normal case
                        icell = (xc-1) * Ncelly * Ncellz + (yc-1) * Ncellz + zc
                        j = head(icell)
                        do 
                            if (j .eq. 0) exit
                            if (xbd .eqv. yes) then
                                dx = Lx - abs(r(j, 1) - r(i, 1)) 
                            else
                                dx = r(j, 1) - r(i, 1)
                            end if
                            if (ybd .eqv. yes) then
                                dy = Ly - abs(r(j, 2) - r(i, 2))
                            else
                                dy = r(j, 2) - r(i, 2)
                            end if
                            if (zbd .eqv. yes) then
                                dz = Lz - abs(r(j, 3) - r(i, 3))
                            else
                                dz = r(j, 3) - r(i, 3)
                            end if
                            D = dsqrt(dx * dx + dy * dy + dz * dz)
                            if ((D .le. rcut).and.(i .ne. j)) then
                                PDpotential = PDpotential + Umic(D, i, j)
                            end if
                            j = list(j)
                        end do
                    end do
                end do
            end do
        end do
        PDpotential = half * PDpotential
    end function PDpotential


! --------------------------------------------------------------------
!    Purpose: calculate total momentum
! --------------------------------------------------------------------
    function Momentum()
        implicit none
        integer :: i
        double precision :: Momentum(3)
        Momentum = zero
        do i = 1, NP
            Momentum = Momentum + mass(i) * half * (v(i, :) + vb(i, :))
        end do
    end function


! --------------------------------------------------------------------
!    Purpose: calculate temperature
! --------------------------------------------------------------------
    function Temp(start, finish)
        implicit none
        double precision :: Temp
        integer :: start, finish, i
        Temp = zero
        do i = start, finish
            Temp = Temp + mass(i) * mag2(half * (v(i, :) + vb(i, :)))
        end do
        Temp = half * Temp
    end function Temp


! --------------------------------------------------------------------
!    Purpose: calculate particle with a sphere
! --------------------------------------------------------------------
    function within(start, finish, center)
        implicit none
        integer :: start, finish, i
        double precision :: within, D, center(3)
        within = zero
        do i = start, finish
            D = mag(r(i, :)-center)
            if (D .le. radius) then
                within = within + one
            end if
        end do
        ! within = within / (dble(finish) - dble(start) + one)
    end function within


! --------------------------------------------------------------------
!    Purpose: avg. number of collision per particles
! --------------------------------------------------------------------
    function numcoll()
        implicit none
        integer :: i, j
        double precision :: D, numcoll
        numcoll = 0
        do i = 1, NP-1
            do j = i+1, NP
                D = mag(r(i, :)-r(j, :))
                if (D .le. rcut) then
                    numcoll = numcoll + 1
                end if
            end do
        end do
        numcoll = numcoll
    end function


! --------------------------------------------------------------------
!    Purpose: output data (r, v, energy, momentum, temperature, #coll)
! --------------------------------------------------------------------
    subroutine output(step)
        implicit none
        integer :: step
        double precision :: center(3)
        ! double precision :: accum=0.d0, P(3)
        open(10, file=trim(loc)//'rx.dat', status='unknown', &
             access='append')
        open(11, file=trim(loc)//'ry.dat', status='unknown', &
             access='append')
        open(12, file=trim(loc)//'rz.dat', status='unknown', &
             access='append')
        open(13, file=trim(loc)//'vx.dat', status='unknown', &
             access='append')
        open(14, file=trim(loc)//'vy.dat', status='unknown', &
             access='append')
        open(15, file=trim(loc)//'vz.dat', status='unknown', &
             access='append')
        open(16, file=trim(loc)//'info.dat', status='unknown', &
             access='append')
        open(17, file=trim(loc)//'P.dat', status='unknown', &
             access='append')
        if (mod(step, dumpstep) .eq. 0) then
            ! output position
            write(10, *) r(:, 1)
            write(11, *) r(:, 2)
            write(12, *) r(:, 3)
            ! output velocity
            write(13, *) v(:, 1)
            write(14, *) v(:, 2)
            write(15, *) v(:, 3)

            ! time, ion kinetic, electron kinetic,
            ! potential, temperature, #coll
            ! accum = accum + numcoll()
            ! P = Momentum()
            center = (/ Lx/two, Ly/two, Lz/two /)
            if (periodic .eqv. yes) then
                write(16, *) time, KE(1, NP1), KE(NP1+1, NP1+NP2), KE(1, NP),&
                            PDpotential(), Temp(1, NP1)/NP1, Temp(NP1+1, NP)/NP2, &
                            within(1, NP1, center)
            else
                write(16, *) time, KE(1, NP1), KE(NP1+1, NP1+NP2), KE(1, NP),&
                            Potential(), Temp(1, NP1)/NP1, Temp(NP1+1, NP)/NP2, &
                            within(1, NP, center)
            end if
            idump = idump + 1
        end if
    end subroutine output


! --------------------------------------------------------------------
!    Purpose: clear output file directory
! --------------------------------------------------------------------
    subroutine remove_file
        implicit none
        integer :: i
        call execute_command_line('python3 logo.py')
        write(6, *) 'Removing file in the output directory in...'
        do i = 2, 1, -1
            write(6, '(I3.1, A4)') i, 'sec'
            call sleep(1)
        end do
        call execute_command_line('rm '//trim(loc)//'*')
        write(6, '(A5)') 'start'
    end subroutine remove_file


! --------------------------------------------------------------------
!    Purpose: print run status
! --------------------------------------------------------------------
    subroutine pstatus(i)
        implicit none
        integer :: i
        if (mod(i, dumpstep) .eq. 0) then
            write(*, '(i5, 1a3, i5)') idump, '/', NumDump
        end if
    end subroutine pstatus


end module mic