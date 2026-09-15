module LocalSweep_mod
    use Dynamics_mod
    use LocalU_mod
    use ObserEqual_mod
    implicit none
    
    public
    private :: Dyn
    
    type :: LocalSweep
    contains
        procedure                   :: init         => Local_sweep_init
        procedure                   :: clear        => Local_sweep_clear
        procedure, private, nopass  :: reset        => Local_sweep_reset
        procedure                   :: therm        => Local_sweep_therm
        procedure                   :: pre          => Local_sweep_pre
        procedure, private, nopass  :: sweep_L      => Local_sweep_L
        procedure, private, nopass  :: sweep_R      => Local_sweep_R
        procedure                   :: sweep        => Local_sweep
        procedure, nopass           :: ctrl_print_l => Local_control_print
        procedure, nopass           :: ctrl_print_t => Therm_control_print
    end type LocalSweep
    
    type(ObserTau), allocatable :: Obs_tau
    type(ObserEqual), allocatable :: Obs_equal
    type(Dynamics) :: Dyn
    
contains
    subroutine Local_sweep_init(this)
        class(LocalSweep), intent(inout) :: this
        allocate(Obs_equal)
        call Obs_equal%make()
        if (is_tau) then
            allocate(Obs_tau)
            call Obs_tau%make()
            call Dyn%init()
        endif
        call LocalU_init(Op_U1)
        call LocalU_init(Op_U2)
        return
    end subroutine Local_sweep_init
    
    subroutine Local_sweep_clear(this)
        class(LocalSweep), intent(inout) :: this
        deallocate(Obs_equal)
        call LocalU_clear(Op_U1)
        call LocalU_clear(Op_U2)
        if (is_tau) then
            deallocate(Obs_tau)
            call Dyn%clear()
        endif
        return
    end subroutine Local_sweep_clear
    
    subroutine Local_sweep_reset(toggle)
        logical, intent(in) :: toggle
        call LocalU_reset(Op_U1)
        call LocalU_reset(Op_U2)
        call Obs_equal%reset()
        if (toggle) call Obs_tau%reset()
        return
    end subroutine Local_sweep_reset
    
    subroutine Local_sweep_therm(this, iseed)
        class(LocalSweep), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer :: ii, nt
        call this%reset(.false.)
        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call LocalU_prop_therm(Op_U1, iseed, 1, nt)
            if (abs(RU2) > Zero) call LocalU_prop_therm(Op_U2, iseed, 2, nt)
        enddo
        if (abs(RU1) > Zero) call Op_U1%Acc_U_therm%ratio()
        if (abs(RU2) > Zero) call Op_U2%Acc_U_therm%ratio()
        return
    end subroutine Local_sweep_therm
    
    subroutine Local_sweep_pre(this, Prop, WrList)
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer :: nt
        call this%reset(.false.)
        call Wrap_pre(Prop, WrList, 0)
        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call propU_pre(Op_U1, Prop, 1, nt)
            if (abs(RU2) > Zero) call propU_pre(Op_U2, Prop, 2, nt)
            call propT_pre(Prop)
            if (mod(nt, Nwrap) == 0) call Wrap_pre(Prop, WrList, nt)
        enddo
        return
    end subroutine Local_sweep_pre
    
    subroutine Local_sweep_L(Prop, WrList, iseed, Nobs)
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: iseed, Nobs
        integer :: nt
        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0) call Wrap_L(Prop, WrList, nt, "S")
            call Obs_equal%calc(Prop, nt)
            Nobs = Nobs + 1
            call propT_L(Prop)
            if (abs(RU2) > Zero) call LocalU_prop_L(Op_U2, Prop, iseed, 2, nt)
            if (abs(RU1) > Zero) call LocalU_prop_L(Op_U1, Prop, iseed, 1, nt)
        enddo
        call Wrap_L(Prop, WrList, 0, "S")
        return
    end subroutine Local_sweep_L
    
    subroutine Local_sweep_R(Prop, WrList, iseed, toggle, Nobs, Nobst)
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        logical, intent(in) :: toggle
        integer, intent(inout) :: iseed, Nobs, Nobst
        integer :: nt
        if (toggle) then ! calculating time-sliced Green's function before sweep_right
            call Dyn%reset(Prop)
            call Dyn%sweep_R(Obs_tau, WrList)
            Nobst = Nobst + 1
        endif
        call Wrap_R(Prop, WrList, 0, "S")
        do nt = 1, Ltrot
            if (abs(RU1) > Zero) call LocalU_prop_R(Op_U1, Prop, iseed, 1, nt)
            if (abs(RU2) > Zero) call LocalU_prop_R(Op_U2, Prop, iseed, 2, nt)
            call propT_R(Prop)
            if (mod(nt, Nwrap) == 0) call Wrap_R(Prop, WrList, nt, "S")
            call Obs_equal%calc(Prop, nt)
            Nobs = Nobs + 1
        enddo
        return
    end subroutine Local_sweep_R
    
    subroutine Local_sweep(this, Prop, WrList, iseed, is_beta, toggle)
        class(LocalSweep), intent(inout) :: this
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: iseed
        logical, intent(inout) :: is_beta
        logical, intent(in) :: toggle
        integer :: Nobs, Nobst, nsw

        call this%reset(toggle)
        Nobs = 0; Nobst = 0
        do nsw = 1, Nsweep
            if(is_beta) then
                call this%sweep_L(Prop, WrList, iseed, Nobs)
                call this%sweep_R(Prop, WrList, iseed, toggle, Nobs, Nobst)
            else
                call this%sweep_R(Prop, WrList, iseed, toggle, Nobs, Nobst)
                call this%sweep_L(Prop, WrList, iseed, Nobs)
            endif
        enddo
        call Obs_equal%ave(Nobs)
        if (toggle) call Obs_tau%ave(Nobst)
        if (abs(RU1) > Zero) call Op_U1%Acc_U_local%ratio()
        if (abs(RU2) > Zero) call Op_U2%Acc_U_local%ratio()
        return
    end subroutine Local_sweep
    
    subroutine Local_control_print(toggle)
        include 'mpif.h'
        logical, intent(in) :: toggle
        real(kind=8) :: collect
        if (abs(RU1) > Zero) then
            collect = 0.d0
            call MPI_Reduce(Op_U1%Acc_U_local%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) then
                Op_U1%Acc_U_local%acc = collect / dble(ISIZE * Nbin)
                write(50,*) 'Accept_U1_shift                                :', Op_U1%Acc_U_local%acc
            endif
        endif
        if (abs(RU2) > Zero) then
            collect = 0.d0
            call MPI_Reduce(Op_U2%Acc_U_local%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) then
                Op_U2%Acc_U_local%acc = collect / dble(ISIZE * Nbin)
                write(50,*) 'Accept_U2_shift                                :', Op_U2%Acc_U_local%acc
            endif
        endif
        if (toggle) call Dyn%ctrl_print()
    end subroutine Local_control_print
    
    subroutine Therm_control_print()
        include 'mpif.h'
        real(kind=8) :: collect
        if (abs(RU1) > Zero) then
            collect = 0.d0
            call MPI_Reduce(Op_U1%Acc_U_therm%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) then
                Op_U1%Acc_U_therm%acc = collect / dble(ISIZE * Nwarm)
                write(50,*) "Thermalize U1 Accept Ratio                     :", Op_U1%Acc_U_therm%acc
            endif
        endif
        if (abs(RU2) > Zero) then
            collect = 0.d0
            call MPI_Reduce(Op_U2%Acc_U_therm%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) then
                Op_U2%Acc_U_therm%acc = collect / dble(ISIZE * Nwarm)
                write(50,*) "Thermalize U2 Accept Ratio                     :", Op_U2%Acc_U_therm%acc
            endif
        endif
        return
    end subroutine Therm_control_print
end module LocalSweep_mod