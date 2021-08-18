program gau3d
    use math
    implicit none
    integer :: i
    double precision :: v(10000, 3)

    do i = 1, 10000
        v(i, :) = maxwell(temperature=10.d0*e, mass=me, dim=3)
    end do

    open(66, file='out3.dat', status='unknown')
    do i = 1, 10000
        write(66, *) mag(v(i, :))
    end do
end program gau3d