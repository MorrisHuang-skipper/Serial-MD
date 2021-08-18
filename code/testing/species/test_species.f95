program test
    use species_mod
    implicit none
    type(species) :: H

    call H % initialize(1.d0, 1.d-31, 1, 2)

    print*, H % charge
    print*, H % mass
end program test