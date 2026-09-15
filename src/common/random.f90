! Original BAFQMC uniform RNG recurrence, expressed without integer overflow.
! This produces the same seed sequence as the historical 32-bit ranf routine:
! seed <- (48828125 * seed) mod 2^31; u = seed / 2^31.
function ranf(iq) result(value)
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    integer, intent(inout) :: iq
    real(real64) :: value
    integer(int64), parameter :: modulus = 2147483648_int64
    iq = int(modulo(48828125_int64 * int(iq, int64), modulus))
    value = real(iq, real64) / real(modulus, real64)
end function ranf
