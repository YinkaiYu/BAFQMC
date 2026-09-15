module PoleDiagnostics_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none

    private
    public :: PoleDiagnostics

    type :: PoleDiagnostics
    contains
        procedure :: write => PoleDiagnostics_write
    end type PoleDiagnostics

contains
    subroutine compute_eigenvalues(Gr, evals)
        complex(kind=8), dimension(:,:), intent(in) :: Gr
        complex(kind=8), dimension(:), intent(out) :: evals
        complex(kind=8), dimension(:,:), allocatable :: matrix
        complex(kind=8), dimension(:), allocatable :: work
        complex(kind=8) :: work_query(1)
        complex(kind=8) :: vl(1, 1), vr(1, 1)
        real(kind=8), dimension(:), allocatable :: rwork
        integer :: info, lwork, n
        external :: zgeev

        n = size(Gr, 1)
        allocate(matrix(n, n))
        allocate(rwork(max(1, 2*n)))
        matrix = Gr
        lwork = -1
        call zgeev('N', 'N', n, matrix, n, evals, vl, 1, vr, 1, &
                   work_query, lwork, rwork, info)
        if (info .ne. 0) then
            write(6,*) 'ERROR: ZGEEV workspace query failed in pole diagnostics', &
                       ' on rank ', IRANK, ' info=', info
            stop 1
        endif

        lwork = max(1, ceiling(real(work_query(1), kind=8)))
        allocate(work(lwork))
        matrix = Gr
        call zgeev('N', 'N', n, matrix, n, evals, vl, 1, vr, 1, &
                   work, lwork, rwork, info)
        if (info .ne. 0) then
            write(6,*) 'ERROR: ZGEEV failed in pole diagnostics on rank ', &
                       IRANK, ' info=', info
            stop 1
        endif
        deallocate(matrix, work, rwork)
        return
    end subroutine compute_eigenvalues

    subroutine compute_singular_values(Gr, svals)
        complex(kind=8), dimension(:,:), intent(in) :: Gr
        real(kind=8), dimension(:), intent(out) :: svals
        complex(kind=8), dimension(:,:), allocatable :: matrix
        complex(kind=8), dimension(:), allocatable :: work
        complex(kind=8) :: work_query(1)
        complex(kind=8) :: u(1, 1), vt(1, 1)
        real(kind=8), dimension(:), allocatable :: rwork
        integer :: info, lwork, n
        external :: zgesvd

        n = size(Gr, 1)
        allocate(matrix(n, n))
        allocate(rwork(max(1, 5*n)))
        matrix = Gr
        lwork = -1
        call zgesvd('N', 'N', n, n, matrix, n, svals, u, 1, vt, 1, &
                    work_query, lwork, rwork, info)
        if (info .ne. 0) then
            write(6,*) 'ERROR: ZGESVD workspace query failed in pole diagnostics', &
                       ' on rank ', IRANK, ' info=', info
            stop 1
        endif

        lwork = max(1, ceiling(real(work_query(1), kind=8)))
        allocate(work(lwork))
        matrix = Gr
        call zgesvd('N', 'N', n, n, matrix, n, svals, u, 1, vt, 1, &
                    work, lwork, rwork, info)
        if (info .ne. 0) then
            write(6,*) 'ERROR: ZGESVD failed in pole diagnostics on rank ', &
                       IRANK, ' info=', info
            stop 1
        endif
        deallocate(matrix, work, rwork)
        return
    end subroutine compute_singular_values

    subroutine calc_local_diagnostics(Prop, zvals, scalar_local)
        class(Propagator), intent(in) :: Prop
        complex(kind=8), dimension(Ndim), intent(out) :: zvals
        real(kind=8), dimension(5), intent(out) :: scalar_local
        complex(kind=8), dimension(Ndim) :: evals
        real(kind=8), dimension(Ndim) :: svals
        real(kind=8) :: eval_abs, small
        real(kind=8) :: pole_distance, pole_x, rho_g, smax_g
        real(kind=8) :: log_weight
        integer :: ia

        small = tiny(1.d0)
        call compute_eigenvalues(Prop%Gr, evals)
        call compute_singular_values(Prop%Gr, svals)

        do ia = 1, Ndim
            eval_abs = abs(evals(ia))
            if (eval_abs .le. small) then
                zvals(ia) = dcmplx(1.d0 / small, 0.d0)
            else
                zvals(ia) = 1.d0 / evals(ia)
            endif
        enddo

        pole_distance = max(minval(abs(zvals)), small)
        pole_x = -log10(pole_distance)
        rho_g = maxval(abs(evals))
        smax_g = maxval(svals)
        log_weight = Conf%log_weight() + 2.d0 * sum(log(max(svals, small)))

        scalar_local = (/ pole_distance, pole_x, rho_g, smax_g, log_weight /)
        return
    end subroutine calc_local_diagnostics

    subroutine append_pole_z(z_collect)
        complex(kind=8), dimension(Ndim, ISIZE), intent(in) :: z_collect
        integer :: ia, irank_out

        open(unit=90, file='pole_z', status='unknown', action='write', &
             position='append')
        do irank_out = 1, ISIZE
            do ia = 1, Ndim
                write(90, '(1X,ES24.16E3,1X,ES24.16E3)', advance='no') &
                    real(z_collect(ia, irank_out), kind=8), &
                    aimag(z_collect(ia, irank_out))
            enddo
            write(90,*)
        enddo
        close(90)
        return
    end subroutine append_pole_z

    subroutine append_scalar_column(filename, scalar_collect, column)
        character(len=*), intent(in) :: filename
        real(kind=8), dimension(5, ISIZE), intent(in) :: scalar_collect
        integer, intent(in) :: column
        integer :: irank_out

        open(unit=91, file=filename, status='unknown', action='write', &
             position='append')
        do irank_out = 1, ISIZE
            write(91, '(ES24.16E3)') scalar_collect(column, irank_out)
        enddo
        close(91)
        return
    end subroutine append_scalar_column

    subroutine PoleDiagnostics_write(this, Prop)
        include 'mpif.h'
        class(PoleDiagnostics), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        complex(kind=8), dimension(Ndim) :: zvals
        complex(kind=8), dimension(:,:), allocatable :: z_collect
        real(kind=8), dimension(5) :: scalar_local
        real(kind=8), dimension(:,:), allocatable :: scalar_collect

        call calc_local_diagnostics(Prop, zvals, scalar_local)

        allocate(z_collect(Ndim, ISIZE))
        allocate(scalar_collect(5, ISIZE))
        call MPI_GATHER(zvals, Ndim, MPI_complex16, z_collect, Ndim, &
                        MPI_complex16, 0, MPI_COMM_WORLD, IERR)
        call MPI_GATHER(scalar_local, 5, MPI_real8, scalar_collect, 5, &
                        MPI_real8, 0, MPI_COMM_WORLD, IERR)

        if (IRANK == 0) then
            call append_pole_z(z_collect)
            call append_scalar_column('pole_distance', scalar_collect, 1)
            call append_scalar_column('pole_x', scalar_collect, 2)
            call append_scalar_column('green_spectral_radius', scalar_collect, 3)
            call append_scalar_column('green_smax', scalar_collect, 4)
            call append_scalar_column('log_weight', scalar_collect, 5)
        endif

        deallocate(z_collect, scalar_collect)
        return
    end subroutine PoleDiagnostics_write
end module PoleDiagnostics_mod
