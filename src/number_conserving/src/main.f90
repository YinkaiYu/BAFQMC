program bosonDQMC
    use LocalSweep_mod
    ! use GlobalUpdate_mod
    use FourierTrans_mod
    use PoleDiagnostics_mod
    implicit none
    include 'mpif.h'
    
    integer :: status(MPI_STATUS_SIZE)
    integer:: iseed, nth, nbc, N
    logical :: is_beta, istau_tmp
    real(kind=8) :: collect, CPUT
    integer(kind=8) :: ICPU_1, ICPU_2, N_P_SEC
    
    ! type(GlobalUpdate) :: Sweep_global
    type(LocalSweep) :: Sweep_local
    type(FourierTrans) :: Fourier
    type(PoleDiagnostics) :: PoleDiag
    type(Propagator), allocatable :: Prop
    type(WrapList), allocatable :: WrList

    call MPI_INIT(IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ISIZE, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, IERR)
    
    call SYSTEM_CLOCK(COUNT_RATE = N_P_SEC)
    call SYSTEM_CLOCK(COUNT = ICPU_1)
! initiate
    call Model_init(iseed)
    allocate(Prop)
    call Prop%make()
    allocate(WrList)
    call WrList%make()
    call Stabilize_init()
    call Sweep_local%init()
    ! if (is_global) call Sweep_global%init()
! boson warm-up
    if (is_warm) then
        do nth = 1, Nwarm
            call Sweep_local%therm(iseed)
        enddo
        call Sweep_local%ctrl_print_t()
    else
        if (IRANK == 0) write(50,*) "Skipping Bosonic warm-up"
    endif
! Sweep
    call Sweep_local%pre(Prop, WrList)
    is_beta = .true.; istau_tmp = .false.
    do nbc = 1, Nbin
        if (nbc .gt. Nthermal) istau_tmp = is_tau
        ! if (is_global) call Sweep_global%sweep(Prop, WrList, iseed, is_beta)
        call Sweep_local%sweep(Prop, WrList, iseed, is_beta, istau_tmp)
        call Fourier%preq(Obs_equal)
        call PoleDiag%write(Prop)
        if (istau_tmp) call Fourier%prtau(Obs_tau)
    enddo
! control print
    ! if (is_global) call Sweep_global%ctrl_print()
    call Sweep_local%ctrl_print_l(istau_tmp)
    collect = 0.d0
    call MPI_Reduce(Prop%Xmaxm, collect, 1, MPI_Real8, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) Prop%Xmaxm = collect
    N = 2* ISIZE * Nbin * Nst * Nsweep
    collect = 0.d0
    call MPI_Reduce(Prop%Xmeanm, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) Prop%Xmeanm = collect / dble(N)
    call SYSTEM_CLOCK(COUNT = ICPU_2)
    CPUT = dble(ICPU_2 - ICPU_1) / dble(N_P_SEC)
    collect = 0.d0
    call MPI_Reduce(CPUT, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
    if (IRANK == 0) CPUT = collect / dble(ISIZE)
    if (IRANK == 0) then
        write(50,*) 'Max diff Matrix                                :', Prop%Xmaxm
        write(50,*) 'Mean diff Matrix                               :', Prop%Xmeanm
        write(50,*) 'Tot CPU time                                   :', CPUT
    endif
! deallocate
    ! if (is_global) call Sweep_global%clear()
    call Sweep_local%clear()
    call Stabilize_clear()
    deallocate(Prop)
    deallocate(WrList)
    call Model_clear(iseed) ! conf-out
    
    call MPI_FINALIZE(IERR)
end program bosonDQMC
