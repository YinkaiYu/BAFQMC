! Small BLAS/LAPACK implementation of the MyMats interface used by BAFQMC.
! Written for this distribution. No Lib_90_new, EISPACK, LINPACK or NAG source
! is included. See README.md for numerical compatibility and validation.
module MyMats
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: mmult, diag, inv, udv
    integer, parameter :: dp = real64
contains
    subroutine require_success(info, operation)
        integer, intent(in) :: info
        character(len=*), intent(in) :: operation
        if (info /= 0) then
            write(*, *) 'LAPACK error in ', operation, ': info = ', info
            error stop 1
        end if
    end subroutine

    subroutine mmult(c, a, b)
        complex(dp), intent(in) :: a(:, :), b(:, :)
        complex(dp), intent(out) :: c(:, :)
        integer :: m, n, k
        m = size(a, 1); k = size(a, 2); n = size(b, 2)
        if (size(b, 1) /= k .or. size(c, 1) /= m .or. size(c, 2) /= n) &
            error stop 'mmult: inconsistent dimensions'
        call zgemm('N', 'N', m, n, k, (1.0_dp, 0.0_dp), a, m, b, k, &
                   (0.0_dp, 0.0_dp), c, m)
    end subroutine

    subroutine diag(a, u, w)
        complex(dp), intent(in) :: a(:, :)
        complex(dp), intent(out) :: u(:, :)
        real(dp), intent(out) :: w(:)
        complex(dp), allocatable :: work(:)
        real(dp), allocatable :: rwork(:)
        complex(dp) :: query(1)
        integer :: n, lwork, info
        n = size(a, 1)
        if (size(a, 2) /= n .or. any(shape(u) /= shape(a)) .or. size(w) /= n) &
            error stop 'diag: inconsistent dimensions'
        allocate(rwork(max(1, 3*n-2)))
        u = a
        call zheev('V', 'L', n, u, n, w, query, -1, rwork, info)
        call require_success(info, 'zheev workspace query')
        lwork = max(1, int(real(query(1), dp)))
        allocate(work(lwork))
        call zheev('V', 'L', n, u, n, w, work, lwork, rwork, info)
        call require_success(info, 'zheev')
    end subroutine

    subroutine inv(a, ainv, det)
        complex(dp), intent(in) :: a(:, :)
        complex(dp), intent(out) :: ainv(:, :), det
        complex(dp), allocatable :: work(:)
        integer, allocatable :: ipiv(:)
        complex(dp) :: query(1)
        integer :: n, i, lwork, info
        n = size(a, 1)
        if (size(a, 2) /= n .or. any(shape(ainv) /= shape(a))) &
            error stop 'inv: inconsistent dimensions'
        allocate(ipiv(n))
        ainv = a
        call zgetrf(n, n, ainv, n, ipiv, info)
        call require_success(info, 'zgetrf')
        det = (1.0_dp, 0.0_dp)
        do i = 1, n
            det = det * ainv(i, i)
            if (ipiv(i) /= i) det = -det
        end do
        call zgetri(n, ainv, n, ipiv, query, -1, info)
        call require_success(info, 'zgetri workspace query')
        lwork = max(1, int(real(query(1), dp)))
        allocate(work(lwork))
        call zgetri(n, ainv, n, ipiv, work, lwork, info)
        call require_success(info, 'zgetri')
    end subroutine

    subroutine udv(a, u, d, v, ncon)
        ! A = U diag(D) V, U has orthonormal columns, D is positive,
        ! and V is upper triangular. Used by the alternative stabilization.
        complex(dp), intent(in) :: a(:, :)
        complex(dp), intent(out) :: u(:, :), d(:), v(:, :)
        integer, intent(in) :: ncon
        complex(dp), allocatable :: tau(:), work(:), scaled_v(:, :)
        complex(dp) :: query(1)
        real(dp) :: scale
        integer :: m, n, i, lwork, info
        m = size(a, 1); n = size(a, 2)
        if (m < n .or. any(shape(u) /= shape(a)) .or. size(d) /= n .or. &
            any(shape(v) /= [n, n])) error stop 'udv: inconsistent dimensions'
        allocate(tau(n))
        u = a
        call zgeqrf(m, n, u, m, tau, query, -1, info)
        call require_success(info, 'zgeqrf workspace query')
        lwork = max(1, int(real(query(1), dp)))
        allocate(work(lwork))
        call zgeqrf(m, n, u, m, tau, work, lwork, info)
        call require_success(info, 'zgeqrf')
        v = (0.0_dp, 0.0_dp)
        do i = 1, n
            scale = abs(u(i, i))
            if (scale == 0.0_dp) error stop 'udv: zero diagonal'
            d(i) = cmplx(scale, 0.0_dp, dp)
            v(i, i:n) = u(i, i:n) / scale
        end do
        call zungqr(m, n, n, u, m, tau, query, -1, info)
        call require_success(info, 'zungqr workspace query')
        lwork = max(1, int(real(query(1), dp)))
        deallocate(work)
        allocate(work(lwork))
        call zungqr(m, n, n, u, m, tau, work, lwork, info)
        call require_success(info, 'zungqr')
        if (ncon == 1) then
            allocate(scaled_v(n, n))
            do i = 1, n
                scaled_v(i, :) = d(i) * v(i, :)
            end do
            write(*, *) 'UDV reconstruction residual: ', maxval(abs(matmul(u, scaled_v)-a))
        end if
    end subroutine
end module MyMats
