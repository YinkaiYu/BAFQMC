module DQMC_Model_mod
    use NonInteract
    use Fields_mod
    use OperatorHubbard_mod
    implicit none
    
    public
    type(triangularLattice), allocatable :: Latt
    type(OperatorKinetic), allocatable  :: Op_T
    type(OperatorHubbard)               :: Op_U1,   Op_U2
    type(AuxConf), allocatable          :: Conf
    
contains
    subroutine Model_init(iseed)
        integer, intent(out) :: iseed
! read in parameters
        call read_input()
        call Params_set()
        call write_info()
! initiate lattice lists
        allocate(Latt)
        call Lattice_make(Latt)
! initiate auxiliary field configuration
        allocate(Conf)
        call Conf%make()
        call conf_in(Conf, iseed, Latt)
! set non-interacting exponential operator
        allocate(Op_T)
        call Op_T%make()
        call Op_T%set(Latt)
        write(50,*) 'highest energy level of single particle Ham     :', Op_T%energy_max
        write(50,*) 'lowest  energy level of single particle Ham     :', Op_T%energy_min
        write(50,*) 'bandwidth            of single particle Ham     :', Op_T%bandwidth
! set H-S exponential
        call Op_U1%set(RU1)
        call Op_U2%set(RU2)
        return
    end subroutine Model_init
    
    subroutine Model_clear(iseed)
        integer, intent(in) :: iseed
        call conf_out(Conf, iseed)
        deallocate(Op_T)
        deallocate(Conf)
        deallocate(Latt)
        return
    end subroutine Model_clear
end module DQMC_Model_mod
