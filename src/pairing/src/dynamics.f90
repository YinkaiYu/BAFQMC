module Dynamics_mod
    use Stabilize_mod
    use Multiply_mod
    use ObserTau_mod
    implicit none
    
    public
    private :: Prop_d, PropGr
    
    type :: Dynamics
    contains
        procedure, nopass :: init => Dyn_init
        procedure, nopass :: clear => Dyn_clear
        procedure, nopass :: reset => Dyn_reset
        procedure, nopass :: sweep_R => Dyn_sweep_R
        procedure, nopass :: ctrl_print => Dyn_control_print
    end type Dynamics
    
    type(Propagator), allocatable :: Prop_d
    type(PropGreen), allocatable :: PropGr
    
contains
    subroutine Dyn_init()
        allocate(Prop_d)
        call Prop_d%make()
        allocate(PropGr)
        call PropGr%make()
        return
    end subroutine Dyn_init
    
    subroutine Dyn_clear()
        deallocate(Prop_d)
        deallocate(PropGr)
        return
    end subroutine Dyn_clear
    
    subroutine Dyn_reset(Prop)
        class(Propagator), intent(in) :: Prop
        call Prop_d%asgn(Prop)
        call PropGr%reset(Prop%Gr) ! the starting point of the time-sliced Green's function
        return
    end subroutine Dyn_reset
    
    subroutine Dyn_sweep_R(Obs, WrList)
        class(ObserTau), intent(inout) :: Obs
        class(WrapList), intent(in) :: WrList
        integer :: nt
        do nt = 1, Ltrot
            call Obs%calc(PropGr, nt)
            if (abs(U1) > Zero) call propgrU_R(Op_U1, Prop_d, PropGr, 1, nt)
            if (abs(U2) > Zero) call propgrU_R(Op_U2, Prop_d, PropGr, 1, nt)
            call propgrT_R(Prop_d, PropGr)
            if (mod(nt, Nwrap) == 0) call Wrap_tau(Prop_d, PropGr, WrList, nt)
        enddo
        return
    end subroutine Dyn_sweep_R
    
    subroutine Dyn_control_print()
        include 'mpif.h'
        integer :: no, N
        real(kind=8) :: collect

        do no = 1, 3
            collect = 0.d0
            call MPI_Reduce(PropGr%Xmaxm(no), collect, 1, MPI_Real8, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) PropGr%Xmaxm(no) = collect
            N = ISIZE * Nbin * Nst * Nsweep
            collect = 0.d0
            call MPI_Reduce(PropGr%Xmeanm(no), collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
            if (IRANK == 0) PropGr%Xmeanm(no) = collect / dble(N)
        enddo
        if (IRANK == 0) then
            write(50,*) 'Max diff GRtt                                  :', PropGr%Xmaxm(1)
            write(50,*) 'Mean diff GRtt                                 :', PropGr%Xmeanm(1)
            write(50,*) 'Max diff GRt0                                  :', PropGr%Xmaxm(2)
            write(50,*) 'Mean diff GRt0                                 :', PropGr%Xmeanm(2)
            write(50,*) 'Max diff GR0t                                  :', PropGr%Xmaxm(3)
            write(50,*) 'Mean diff GR0t                                 :', PropGr%Xmeanm(3)
        endif
        return
    end subroutine Dyn_control_print
end module Dynamics_mod