module species_mod
    implicit none

    !> base species attribute
    type species
        double precision :: charge
        double precision :: mass
        integer :: id_start
        integer :: id_end

    contains
        procedure :: initialize
    end type species

    contains

    !> initialize species objects
    subroutine initialize(sp, charge, mass, id_start, id_end)
        class(species) :: sp
        double precision :: charge
        double precision :: mass
        integer :: id_start
        integer :: id_end

        sp % charge = charge
        sp % mass = mass
        sp % id_start = id_start
        sp % id_end = id_end
        
    end subroutine initialize

end module species_mod