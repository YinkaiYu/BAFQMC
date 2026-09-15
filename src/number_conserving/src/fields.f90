module Fields_mod
    use MyLattice
    implicit none
    
    public :: AuxConf, conf_in, conf_out
    private
    
    type :: AuxConf
        real(kind=8), dimension(:,:,:), public, allocatable :: phi_list
    contains
        procedure   :: make =>  AuxConf_make
        procedure   :: log_weight => AuxConf_log_weight
        final       ::          AuxConf_clear
    end type AuxConf

contains
    subroutine AuxConf_make(this)
        class(AuxConf), intent(inout) :: this
        allocate(this%phi_list(Naux, Ndim, Ltrot))
        return
    end subroutine AuxConf_make
    
    subroutine AuxConf_clear(this)
        type(AuxConf), intent(inout) :: this
        deallocate(this%phi_list)
        return
    end subroutine AuxConf_clear

    real(kind=8) function AuxConf_log_weight(this) result(log_weight)
        class(AuxConf), intent(in) :: this
        log_weight = -0.5d0 * sum(this%phi_list * this%phi_list)
        return
    end function AuxConf_log_weight
    
    subroutine conf_in(Conf, iseed, Latt)
        class(AuxConf), intent(inout) :: Conf
        integer, intent(out) :: iseed
        class(triangularLattice), intent(in) :: Latt
! Local: 
        real(kind=8), dimension(Naux, NdimTherm, LtrotTherm) :: phi_list_therm
        
        call conf_therm_in(phi_list_therm, iseed)
        call conf_transfer(Conf%phi_list, phi_list_therm, Latt)
        return
    end subroutine conf_in
    
    subroutine conf_therm_in(phi_list_therm, iseed)
!#define DEC
        include 'mpif.h'
! Arguments: 
        real(kind=8), dimension(Naux, NdimTherm, LtrotTherm), intent(inout) :: phi_list_therm
        integer, intent(out) :: iseed
! Local: 
        integer :: status(MPI_STATUS_SIZE)
        integer :: iseed0, itmp
        integer :: ii, ns, nt, N
        real(kind=8), external :: ranf
        real(kind=8), dimension(Naux, NdimTherm, LtrotTherm) :: phi_list_itmp
        
        if (IRANK == 0 ) then
            open(unit=30, file='confin.txt', status='unknown')
            open(unit=10, file='seeds.txt', status='unknown')
        endif
        if ( IRANK == 0 ) then
            write(6,*) 'Number of process', ISIZE
            read(30,*) iseed ! read seeds from "confin"
            if (iseed == 0) then
                read(10,*) iseed0 ! if confin=0, read the first row of "seeds" as the seed of the master process.
                do N = 1, ISIZE - 1
!   Setup node I and send data.
                    read(10,*) itmp ! read the following rows of "seeds" as the seeds for other nodes. Each node has only one seed to generate the space-time auxiliary fields
                    call conf_set(phi_list_itmp, itmp)
                    call MPI_SEND(itmp, 1, MPI_Integer, N, N, MPI_COMM_WORLD, IERR)
					call MPI_SEND(phi_list_itmp, NdimTherm*Naux*LtrotTherm, MPI_Real8, N, N+1024, MPI_COMM_WORLD, IERR)
                    print '("after sending message from master ", i3.1, " to Rank ", i3.1)', IRANK, N
                enddo
! Set node zero.
                iseed = iseed0
                call conf_set(phi_list_therm, iseed)
                print '("initiating master rank ", i3.1)', IRANK
            else
! read all confins from NODE 0.
! Setup Node 0
                do nt = 1, LtrotTherm
                    do ii = 1, NdimTherm
                        do ns = 1, Naux
                            read(30,*) phi_list_therm(ns, ii, nt)
                        enddo
                    enddo
                enddo
                do N = 1, ISIZE - 1
                    read(30,*) itmp
                    do nt = 1, LtrotTherm
                        do ii = 1, NdimTherm
                            do ns = 1, Naux
                                read(30,*) phi_list_itmp(ns, ii, nt)
                            enddo
                        enddo
                    enddo
                    call MPI_SEND(itmp, 1, MPI_Integer, N, N, MPI_COMM_WORLD, IERR)
					call MPI_SEND(phi_list_itmp, Naux*NdimTherm*LtrotTherm, MPI_Real8, N, N+1024, MPI_COMM_WORLD, IERR)
                enddo
            endif
        else
            call MPI_RECV(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, STATUS, IERR)
			call MPI_RECV(phi_list_therm, Naux*NdimTherm*LtrotTherm, MPI_Real8, 0, IRANK + 1024, MPI_COMM_WORLD, STATUS, IERR)
            print '("Rank ", i3.1, " after receiving message. ")', IRANK
        endif
        if (IRANK == 0 ) then
            close(30); close(10)
        endif
        return
    end subroutine conf_therm_in

    subroutine conf_set(phi_list, itmp)
        real(kind=8), dimension(Naux, NdimTherm, LtrotTherm), intent(out) :: phi_list
        integer, intent(inout) :: itmp
! Local: 
        integer :: ii, ns, nt
        real(kind=8) :: X
        real(kind=8), external :: ranf
        
        if (iniType == 1) then
            do nt = 1, LtrotTherm
                do ii = 1, NdimTherm
                    do ns = 1, Naux
                        X = ranf(itmp)
                        phi_list(ns, ii, nt) = iniBias(ns) + iniAmpl * (X - 0.5)
                    enddo
                enddo
            enddo
        elseif (iniType == 2) then
            do nt = 1, LtrotTherm
                do ii = 1, NdimTherm
                    do ns = 1, Naux
                        X = rng_Gaussian(itmp)
                        phi_list(ns, ii, nt) = iniBias(ns) + iniAmpl * X
                    enddo
                enddo
            enddo
        elseif (iniType == 3) then
            do nt = 1, LtrotTherm
                do ii = 1, NdimTherm
                    do ns = 1, Naux
                        X = rng_Gaussian(itmp)
                        phi_list(ns, ii, nt) = iniBias(ns) + iniAmpl * abs(X)
                    enddo
                enddo
            enddo
        else
            write(6,*) "incorrect initype input", initype, " on rank ", IRANK
        endif
        return
    end subroutine conf_set
    
    subroutine conf_transfer(phi_list, phi_list_therm, Latt)
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(inout) :: phi_list
        real(kind=8), dimension(Naux, NdimTherm, LtrotTherm), intent(in) :: phi_list_therm
        class(triangularLattice), intent(in) :: Latt
        integer :: nt, ntt, nx, ny, ii, iit, nc, no, n
        integer :: cell_list_therm(LqTherm, 1:2), inv_cell_list_therm(NlxTherm, NlyTherm)
        integer :: dim_list_therm(NdimTherm, 1:2), inv_dim_list_therm(LqTherm, Norb)
! set therm list
        nc = 0
        do no = 1, Norb
            do n = 1, LqTherm
                nc = nc + 1
                dim_list_therm(nc, 1) = n
                dim_list_therm(nc, 2) = no
                inv_dim_list_therm(n, no) = nc
            enddo
        enddo
        nc = 0
        do ny = 1, NlyTherm
            do nx = 1, NlxTherm
                nc = nc + 1
                cell_list_therm(nc, 1) = nx
                cell_list_therm(nc, 2) = ny
                inv_cell_list_therm(nx, ny) = nc
            enddo
        enddo
! copy Conf
        do no = 1, Norb
            do nt = 1, LtrotTherm
                do ny = 1, NlyTherm
                    do nx = 1, NlxTherm
                        iit = inv_dim_list_therm(inv_cell_list_therm(nx, ny), no)
                        ii  = Latt%inv_dim_list(  Latt%inv_cell_list(nx, ny), no)
                        phi_list(:, ii, nt) = phi_list_therm(:, iit, nt)
                    enddo
                enddo
            enddo
        enddo
        do no = 1, Norb
            do nt = 1, LtrotTherm
                do ny = 1, NlyTherm
                    do nx = NlxTherm+1, Nlx
                        iit = inv_dim_list_therm(inv_cell_list_therm(nx - NlxTherm, ny), no)
                        ii  = Latt%inv_dim_list(  Latt%inv_cell_list(nx, ny), no)
                        phi_list(:, ii, nt) = phi_list_therm(:, iit, nt)
                    enddo
                enddo
            enddo
        enddo
        do no = 1, Norb
            do nt = 1, LtrotTherm
                do ny = NlyTherm+1, Nly
                    do nx = 1, NlxTherm
                        iit = inv_dim_list_therm(inv_cell_list_therm(nx, ny - NlyTherm), no)
                        ii  = Latt%inv_dim_list(  Latt%inv_cell_list(nx, ny), no)
                        phi_list(:, ii, nt) = phi_list_therm(:, iit, nt)
                    enddo
                enddo
            enddo
        enddo
        do no = 1, Norb
            do nt = 1, LtrotTherm
                do ny = NlyTherm+1, Nly
                    do nx = NlxTherm+1, Nlx
                        iit = inv_dim_list_therm(inv_cell_list_therm(nx - NlxTherm, ny - NlyTherm), no)
                        ii  = Latt%inv_dim_list(  Latt%inv_cell_list(nx, ny), no)
                        phi_list(:, ii, nt) = phi_list_therm(:, iit, nt)
                    enddo
                enddo
            enddo
        enddo
        do ii = 1, Ndim
            do nt = LtrotTherm+1, Ltrot
                ntt = nt - LtrotTherm
                phi_list(:, ii, nt) = phi_list(:, ii, ntt)
            enddo
        enddo
        return
    end subroutine conf_transfer

    subroutine conf_out(Conf, iseed)
!#define DEC
        include 'mpif.h'
! Arguments: 
        class(AuxConf), intent(in) :: Conf
        integer, intent(in) :: iseed
! Local:
        integer :: ii, ns, nt, N
        integer :: status(MPI_STATUS_SIZE)

        if (IRANK == 0) open(unit=35, file='confout.txt', status='unknown')
        if (IRANK .NE. 0) then
            call MPI_SEND(iseed, 1, MPI_Integer, 0, IRANK, MPI_COMM_WORLD, IERR)
			call MPI_SEND(Conf%phi_list, Ndim*Naux*Ltrot, MPI_Real8, 0, IRANK+1024, MPI_COMM_WORLD, IERR)
            print '("Rank", i3.1, " after sending message. ")', IRANK
        endif
        if (IRANK == 0) then
            write(35,*) iseed
            do nt = 1, Ltrot
                do ii = 1, Ndim
                    do ns = 1, Naux
                        write(35,*) Conf%phi_list(ns, ii, nt)
                    enddo
                enddo
            enddo
            print '("Rank ", i3.1, " after writing down own auxfields. ")', IRANK
            do N = 1, ISIZE - 1
                call MPI_RECV(iseed, 1, MPI_Integer, N, N, MPI_COMM_WORLD, STATUS, IERR)
				call MPI_RECV(Conf%phi_list, Ndim*Naux*Ltrot, MPI_Real8, N, N+1024, MPI_COMM_WORLD, STATUS, IERR)
                write(35,*) iseed
                do nt = 1, Ltrot
                    do ii = 1, Ndim
                        do ns = 1, Naux
                            write(35,*) Conf%phi_list(ns, ii, nt)
                        enddo
                    enddo
                enddo
                print '("Master rank ", i3.1, " after writing down auxfields on RANK ",i3.1 ,".")', IRANK, N
            enddo
        endif
        if (IRANK == 0) close(35)
        return
    end subroutine conf_out
    
    real(kind=8) function rng_Gaussian(iseed) result(X)
        integer, intent(inout) :: iseed
        real(kind=8), external :: ranf
        real(kind=8) :: X1, X2
        X1 = ranf(iseed)
        X2 = ranf(iseed)
        X = sqrt(-2.0*log(X1))*cos(2.0*Pi*X2) ! Gaussian distribution
        return
    end function rng_Gaussian
end module Fields_mod
