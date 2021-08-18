! --------------------------------------------------------------------
!    Purpose: main program
! --------------------------------------------------------------------
program main
    use mic
    use math
    implicit none
    integer :: i=0
    double precision :: start, finish, t1, t2, t3, t0, t00

    !> read input file to initialize
    open(66, file='../InputFile/input.nml')
    read(66, nml=input)
    read(66, nml=laser)
    read(66, nml=derived)
    allocate(r(NP, 3))
    allocate(v(NP, 3))
    allocate(vb(NP, 3))
    allocate(head(totCell))
    allocate(list(NP))

    call remove_file
    call init
    call links
    call halfleap
    t1 = 0
    t2 = 0
    t3 = 0
    call cpu_time(t0)

    do i = 1, Nstep

        call cpu_time(start)
        call move(method='LF')
        call cpu_time(finish)
        t2 = t2 + finish-start

        call cpu_time(start)
        call links
        call cpu_time(finish)
        t1 = t1 + finish-start

        call cpu_time(start)
        call output(i)
        call cpu_time(finish)
        t3 = t3 + finish-start

        call pstatus(i)
    end do

    call cpu_time(t00)
    write(6, *) 'runtime', t00-t0
    write(6, *) 'links', t1, 'move', t2, 'output', t3

end program main
