program test_numerics
    use MyMats
    use, intrinsic :: iso_fortran_env, only: real64, int64
    implicit none
    integer, parameter :: dp = real64, n = 3
    complex(dp) :: a(n,n), u(n,n), v(n,n), d(n), inverse(n,n), product(n,n), identity(n,n), det
    real(dp) :: w(n), residual, value
    real(dp), external :: ranf
    integer :: i, seed
    integer(int64) :: expected
    a = reshape([cmplx(4.0_dp,0.0_dp,dp), cmplx(1.0_dp,-0.7_dp,dp), cmplx(0.2_dp,0.1_dp,dp), &
                 cmplx(1.0_dp,0.7_dp,dp), cmplx(3.0_dp,0.0_dp,dp), cmplx(-0.3_dp,-0.2_dp,dp), &
                 cmplx(0.2_dp,-0.1_dp,dp), cmplx(-0.3_dp,0.2_dp,dp), cmplx(2.0_dp,0.0_dp,dp)], [n,n])
    identity = (0.0_dp,0.0_dp)
    do i = 1,n
        identity(i,i) = (1.0_dp,0.0_dp)
    end do
    call diag(a,u,w)
    call mmult(product,conjg(transpose(u)),u)
    if (maxval(abs(product-identity)) > 1.0e-12_dp) error stop 'eigenvector orthogonality'
    v = u
    do i = 1,n
        v(:,i) = v(:,i)*w(i)
    end do
    call mmult(product,v,conjg(transpose(u)))
    if (maxval(abs(product-a)) > 1.0e-12_dp) error stop 'Hermitian eigendecomposition'
    ! Nonsymmetric complex inversion is needed by the pairing kinetic operator.
    a(1,2) = a(1,2) + cmplx(0.3_dp,0.2_dp,dp)
    call inv(a,inverse,det)
    call mmult(product,a,inverse)
    if (maxval(abs(product-identity)) > 1.0e-12_dp) error stop 'matrix inverse'
    if (abs(det-(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
                 -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
                 +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))) > 1.0e-12_dp) error stop 'determinant'
    call udv(a,u,d,v,0)
    do i = 1,n
        v(i,:) = d(i)*v(i,:)
    end do
    call mmult(product,u,v)
    residual = maxval(abs(product-a))
    if (residual > 1.0e-12_dp) error stop 'UDV reconstruction'
    call mmult(product,conjg(transpose(u)),u)
    if (maxval(abs(product-identity)) > 1.0e-12_dp) error stop 'QR orthogonality'
    ! Literal checkpoints independently evaluated for the historical recurrence.
    seed = 1
    value = ranf(seed)
    if (seed /= 48828125) error stop 'RNG first checkpoint'
    value = ranf(seed)
    if (seed /= 52882121) error stop 'RNG second checkpoint'
    expected = int(seed,int64)
    do i = 1,10000
        expected = modulo(48828125_int64*expected,2147483648_int64)
        value = ranf(seed)
        if (int(seed,int64) /= expected) error stop 'RNG recurrence'
        if (value < 0.0_dp .or. value >= 1.0_dp) error stop 'RNG range'
    end do
    print *, 'PASS: eigenvectors, complex inverse/determinant, UDV, and 10002 RNG draws'
end program
