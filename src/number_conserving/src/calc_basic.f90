module CalcBasic ! Global parameters
    implicit none

! constants
    real(kind=8),           parameter           :: Zero = 1.0d-10
    real(kind=8),           parameter           :: PI = acos(-1.d0)
    real(kind=8),           parameter           :: upbound = 1.0d+200
! lattice parameters
    integer,                parameter           :: Norb  = 1  ! triangular lattice
    integer,                parameter           :: Nsub  = 1  ! triangular lattice
    integer,                parameter           :: Nbond = 3  ! triangular lattice
    integer,                parameter           :: Naux  = 2  ! flavor number of auxiliary field, respectively for U1 term and U2 term
    integer,                public              :: Nlx, Nly, NlxTherm, NlyTherm
    integer,                public              :: Lq, LqTherm
    integer,                public              :: Ndim, NdimTherm
    real(kind=8),           public              :: Dtau
    real(kind=8),           public              :: Beta
    integer,                public              :: Ltrot, LtrotTherm
! Hamiltonian parameters
    real(kind=8),           public,     save    :: RT
    real(kind=8),           public,     save    :: RU1, RU2
    real(kind=8),           public,     save    :: mu
! update parameters
    real(kind=8),           public              :: shiftLoc
    ! logical,                public              :: is_global ! global update switch: Wolff-shift joint update
    ! integer,                public              :: Nglobal
    ! real(kind=8),           public              :: shiftGlb(Naux)
! initial state parameters
    integer,                public              :: iniType ! type of initial phonon field configuration
    real(kind=8),           public              :: iniAmpl ! Gaussian amplitude of initial phonon fields
    real(kind=8),           public              :: iniBias(Naux) ! initial phonon field balance position
! process control parameters
    logical,                    public          :: is_tau ! whether to calculate time-sliced Green function
    integer,                    public          :: Nthermal ! calculate time-sliced Green function from the (Nthermal + 1)-th bin
    logical,                    public          :: is_warm ! bosonic warm-up switch
    integer,                    public          :: Nwarm
    real(kind=8),               public          :: shiftWarm(Naux)
    integer,                    public          :: Nst
    integer,                    public          :: Nwrap
    integer,                    public          :: Nbin
    integer,                    public          :: Nsweep
    integer,                    public          :: ISIZE, IRANK, IERR
    
contains
    subroutine read_input()
        include 'mpif.h'
        if (IRANK == 0) then
            open(unit=20, file='paramC_sets.txt', status='unknown')
            read(20,*) RU1, RU2, mu
            read(20,*) Nlx, Nly, Ltrot, Beta
            read(20,*) NlxTherm, NlyTherm, LtrotTherm
            read(20,*) Nwrap, Nbin, Nsweep, shiftLoc
            read(20,*) is_tau, Nthermal
            read(20,*) is_warm, Nwarm, shiftWarm(1), shiftWarm(2)
            ! read(20,*) is_global, Nglobal, shiftGlb(1), shiftGlb(2)
            read(20,*) iniType, iniAmpl, iniBias(1), iniBias(2)
            close(20)
        endif 
!   MPI process: parallelization
        call MPI_BCAST(Beta, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RU1, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(RU2, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(mu, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(shiftLoc, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        ! call MPI_BCAST(shiftGlb, Naux, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(shiftWarm, Naux, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniAmpl, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniBias, Naux, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(iniType, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nlx, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nly, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Ltrot, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlxTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlyTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(LtrotTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwrap, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nbin, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwarm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        ! call MPI_BCAST(Nglobal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nthermal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nsweep, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_tau, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_warm, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        ! call MPI_BCAST(is_global, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        return
    end subroutine read_input
    
    subroutine Params_set()
        RT = 1.d0 
        Lq = Nlx * Nly
        LqTherm = NlxTherm * NlyTherm
        Dtau = Beta / dble(Ltrot)
        Ndim = Lq * Norb
        NdimTherm = LqTherm * Norb
        if (mod(Ltrot, Nwrap) == 0) then
            Nst = Ltrot / Nwrap
        else
            write(6,*) "Ltrot is not a multiple of Nwrap"; stop
        endif
        return
    end subroutine Params_set
    
    integer function nranf(iseed, N)
        integer, intent(inout) :: iseed
        integer, intent(in) :: N
        real (kind=8), external :: ranf
        nranf  = nint(ranf(iseed)*dble(N) + 0.5)
        if (nranf .lt. 1 ) nranf = 1
        if (nranf .gt. N ) nranf = N 
        return
    end function nranf
    
    ! define the periodic boundary condition in x,y-directions.
    integer function npbc(nr, L)
        integer, intent(in) :: nr, L
        npbc = nr
        if (nr .gt. L) npbc = nr - L
        if (nr .lt. 1) npbc = nr + L
        return
    end function npbc
    
    pure function sqr_vec(vec)
        real(kind=8) :: sqr_vec
        real(kind=8), dimension(Naux), intent(in) :: vec
        integer :: ns
        sqr_vec = 0.d0
        do ns = 1, Naux
            sqr_vec = sqr_vec + vec(ns) * vec(ns)
        enddo
        return
    end function sqr_vec
    
    subroutine write_info()
        if (IRANK == 0) then
            open (unit=50, file='info.txt', status='unknown', action="write")
            write(50,*) '========================='
            write(50,*) 'DQMC for boson Hubbard on triangular lattice'
            write(50,*) 'Linear lengh Lx                                :', Nlx
            write(50,*) 'Linear lengh Ly                                :', Nly
            write(50,*) 'Hopping t                                      :', RT
            write(50,*) 'Hubbard U1                                     :', RU1
            write(50,*) 'Hubbard U2                                     :', RU2
            write(50,*) 'chemical potential                             :', mu
            write(50,*) 'Local update auxiliary field magnitude Shift   :', shiftLoc
            ! if (is_global) then
            ! write(50,*) '# Global                                       :', Nglobal ! global Metropolis algorithm
            ! write(50,*) 'Global shift auxiliary field magnitude for U1  :', shiftGlb(1)
            ! write(50,*) 'Global shift auxiliary field magnitude for U2  :', shiftGlb(2)
            ! endif
            if (is_warm) then
            write(50,*) '# Warm                                         :', Nwarm
            write(50,*) 'Thermalize auxiliary field magnitude for U1    :', shiftWarm(1)
            write(50,*) 'Thermalize auxiliary field magnitude for U2    :', shiftWarm(2)
            endif
            write(50,*) 'Choosing initial distribution type             :', iniType
            write(50,*) 'Magnitude of Gaussian distribution             :', iniAmpl
            write(50,*) 'Beta                                           :', Beta
            write(50,*) 'Trotter number                                 :', Ltrot
            write(50,*) '=>Dtau                                         :', Dtau
            write(50,*) 'N_Ortho                                        :', Nwrap
            write(50,*) '# Bins                                         :', Nbin
            write(50,*) '# Nsweep in one bin                            :', Nsweep
            write(50,*) '# Cores                                        :', ISIZE
            if (mod(Ltrot, Nwrap) .NE. 0) then
                write(50,*) 'Ltrot is not a multiple of Nwrap'; stop
            endif
        endif
        return
    end subroutine write_info
    
end module CalcBasic
