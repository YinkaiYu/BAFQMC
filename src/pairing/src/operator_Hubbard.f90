module OperatorHubbard_mod
    use MyLattice
    implicit none
    public

    ! The public Hamiltonian convention is
    !   U1 (n_b-n_c)^2 + U2 (n_b+n_c)^2,
    ! with U1 >= 0 and U2 <= 0.  Keep the channel identifiers explicit:
    ! the sign of a coupling no longer determines which operator it is.
    integer, parameter :: CHANNEL_RELATIVE = 1
    integer, parameter :: CHANNEL_TOTAL    = 2
    
    type :: AccCounter
        real(kind=8), private :: NC_eff_up, ACC_eff_up
        real(kind=8), public :: acc
    contains
        procedure :: init  => Acc_init
        procedure :: reset => Acc_reset
        procedure :: count => Acc_count
        procedure :: ratio => Acc_calc_ratio
    end type AccCounter

    type :: PhaseCounter
        real(kind=8), private :: NC_prop_bin, NC_acc_bin
        real(kind=8), private :: sum_abs_prop_bin, sum_abs_acc_bin
        real(kind=8), private :: branch_pi_count_bin
        real(kind=8), public :: count_proposed, count_accepted
        real(kind=8), public :: sum_abs_arg_proposed, sum_abs_arg_accepted
        real(kind=8), public :: branch_pi_count
        real(kind=8), public :: mean_abs_arg_proposed
        real(kind=8), public :: mean_abs_arg_accepted
        real(kind=8), public :: branch_pi_fraction
    contains
        procedure :: init  => Phase_init
        procedure :: reset => Phase_reset
        procedure :: count => Phase_count
        procedure :: ratio => Phase_calc_ratio
    end type PhaseCounter
    
    type :: OperatorHubbard
        integer, private :: IUflag
        real(kind=8), private :: Uvalue
        complex(kind=8), private :: alpha ! = sqrt(-2UΔτ)
        complex(kind=8), private :: gaussian
        complex(kind=8), private :: expalpha
        complex(kind=8), private :: expalpha_minus

        complex(kind=8), public :: Delta(Nsec, Nsec)
        complex(kind=8), public :: ratio_gaussian
        complex(kind=8), public :: ratio_constant
        type(AccCounter), public :: Acc_U_local, Acc_U_therm
        type(PhaseCounter), public :: Phase_U_local
    contains
        procedure           :: set          => opU_set
        procedure, private  :: get_exp      => opU_get_exp
        procedure           :: get_delta    => opU_get_delta
        procedure           :: mmult_R      => opU_mmult_R
        procedure           :: mmult_L      => opU_mmult_L
    end type OperatorHubbard
    
contains
    subroutine opU_set(this, Uvalue, channel)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: Uvalue
        integer, intent(in) :: channel
        this%IUflag = 0
        this%Uvalue = 0.0d0
        this%alpha = dcmplx( 0.d0, 0.d0 )
        if (channel /= CHANNEL_RELATIVE .and. channel /= CHANNEL_TOTAL) then
            error stop "invalid pairing HS channel"
        endif
        if (channel == CHANNEL_RELATIVE .and. Uvalue < -Zero) then
            error stop "paper U1 (relative-density coupling) must be non-negative"
        endif
        if (channel == CHANNEL_TOTAL .and. Uvalue > Zero) then
            error stop "paper U2 (total-density coupling) must be non-positive"
        endif
        if (abs(Uvalue) > Zero) then
            this%IUflag = channel
            this%Uvalue = Uvalue
            if (channel == CHANNEL_TOTAL) then
                this%alpha = dcmplx( sqrt(-2.d0 * Uvalue * Dtau), 0.d0 )
            else
                this%alpha = dcmplx( 0.d0, sqrt(2.d0 * Uvalue * Dtau) )
            endif
        endif
        return
    end subroutine opU_set
    
    subroutine opU_get_exp(this, phi, nflag)
        class(OperatorHubbard), intent(inout) :: this
        integer, intent(in) :: nflag ! +1 or -1; propagating direction
        real(kind=8), intent(in) :: phi ! space time local auxiliary field value, for phi_1 or phi_2
        this%gaussian = dcmplx( exp(-0.5d0 * phi * phi), 0.d0 )
        this%expalpha = exp( this%alpha * phi * nflag )
        this%expalpha_minus = exp( - this%alpha * phi * nflag )
        return
    end subroutine opU_get_exp
    
    subroutine opU_get_delta(this, phi_old, phi_new)
        class(OperatorHubbard), intent(inout) :: this
        real(kind=8), intent(in) :: phi_old, phi_new
        real(kind=8) :: delta_phi
        complex(kind=8) :: gaussian_old, gaussian_new
        complex(kind=8) :: expalpha_old, expalpha_new
        complex(kind=8) :: expdelta_p, expdelta_m
        call this%get_exp(phi_old, 1)
        gaussian_old = this%gaussian
        expalpha_old = this%expalpha
        call this%get_exp(phi_new, 1)
        gaussian_new = this%gaussian
        expalpha_new = this%expalpha
        this%ratio_gaussian = gaussian_new / gaussian_old
        this%ratio_constant = dcmplx(1.d0, 0.d0)
        if (this%IUflag == CHANNEL_TOTAL) then
            this%ratio_constant = expalpha_old / expalpha_new ! Nambu normal-ordering constant
        endif
        ! calculate the Delta matrix for ratio_det
        this%Delta = dcmplx(0.d0,0.d0)
        delta_phi = phi_new - phi_old
        call this%get_exp(delta_phi, 1)
        expdelta_p = this%expalpha        ! exp(   alpha * (phi'-phi) )
        expdelta_m = this%expalpha_minus  ! exp( - alpha * (phi'-phi) )
        if (this%IUflag == CHANNEL_TOTAL) then ! U2 total-density term
            this%Delta(1,1) = expdelta_p - dcmplx(1.d0, 0.d0)
            this%Delta(2,2) = expdelta_p - dcmplx(1.d0, 0.d0)
            this%Delta(3,3) = expdelta_m - dcmplx(1.d0, 0.d0)
            this%Delta(4,4) = expdelta_m - dcmplx(1.d0, 0.d0)
        endif
        if (this%IUflag == CHANNEL_RELATIVE) then ! U1 relative-density term
            this%Delta(1,1) = expdelta_p - dcmplx(1.d0, 0.d0)
            this%Delta(2,2) = expdelta_m - dcmplx(1.d0, 0.d0)
            this%Delta(3,3) = expdelta_m - dcmplx(1.d0, 0.d0)
            this%Delta(4,4) = expdelta_p - dcmplx(1.d0, 0.d0)
        endif
        return
    end subroutine opU_get_delta
    
    subroutine opU_mmult_R(this, Mat, Latt, phi, i_site, nflag)
! Arguments: 
        class(OperatorHubbard), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(pairingTriangularLattice), intent(in) :: Latt
        real(kind=8), intent(in) :: phi
        integer, intent(in) :: i_site, nflag
! Local: 
        integer :: i1, i2, i3, i4

        i1 = Latt%inv_dim_list(i_site, 1)
        i2 = Latt%inv_dim_list(i_site, 2)
        i3 = Latt%inv_dim_list(i_site, 3)
        i4 = Latt%inv_dim_list(i_site, 4)

        call this%get_exp(phi, nflag)
        if (this%IUflag == CHANNEL_TOTAL) then ! U2 total-density term
            Mat(i1,1:Ndim) = this%expalpha * Mat(i1,1:Ndim)
            Mat(i2,1:Ndim) = this%expalpha * Mat(i2,1:Ndim)
            Mat(i3,1:Ndim) = this%expalpha_minus * Mat(i3,1:Ndim)
            Mat(i4,1:Ndim) = this%expalpha_minus * Mat(i4,1:Ndim)
        endif
        if (this%IUflag == CHANNEL_RELATIVE) then ! U1 relative-density term
            Mat(i1,1:Ndim) = this%expalpha * Mat(i1,1:Ndim)
            Mat(i2,1:Ndim) = this%expalpha_minus * Mat(i2,1:Ndim)
            Mat(i3,1:Ndim) = this%expalpha_minus * Mat(i3,1:Ndim)
            Mat(i4,1:Ndim) = this%expalpha * Mat(i4,1:Ndim)
        endif
        return
    end subroutine opU_mmult_R
    
    subroutine opU_mmult_L(this, Mat, Latt, phi, j_site, nflag)
! Arguments: 
        class(OperatorHubbard), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(pairingTriangularLattice), intent(in) :: Latt
        real(kind=8), intent(in) :: phi
        integer, intent(in) :: j_site, nflag
! Local: 
        integer :: j1, j2, j3, j4

        j1 = Latt%inv_dim_list(j_site, 1)
        j2 = Latt%inv_dim_list(j_site, 2)
        j3 = Latt%inv_dim_list(j_site, 3)
        j4 = Latt%inv_dim_list(j_site, 4)

        call this%get_exp(phi, nflag)
        if (this%IUflag == CHANNEL_TOTAL) then ! U2 total-density term
            Mat(1:Ndim,j1) = Mat(1:Ndim,j1) * this%expalpha
            Mat(1:Ndim,j2) = Mat(1:Ndim,j2) * this%expalpha
            Mat(1:Ndim,j3) = Mat(1:Ndim,j3) * this%expalpha_minus
            Mat(1:Ndim,j4) = Mat(1:Ndim,j4) * this%expalpha_minus
        endif
        if (this%IUflag == CHANNEL_RELATIVE) then ! U1 relative-density term
            Mat(1:Ndim,j1) = Mat(1:Ndim,j1) * this%expalpha
            Mat(1:Ndim,j2) = Mat(1:Ndim,j2) * this%expalpha_minus
            Mat(1:Ndim,j3) = Mat(1:Ndim,j3) * this%expalpha_minus
            Mat(1:Ndim,j4) = Mat(1:Ndim,j4) * this%expalpha
        endif
        return
    end subroutine opU_mmult_L
    
    subroutine Acc_init(this)
        class(AccCounter), intent(inout) :: this
        this%acc = 0.d0
        return
    end subroutine Acc_init
    
    subroutine Acc_reset(this)
        class(AccCounter), intent(inout) :: this
        this%NC_eff_up = 0.d0
        this%ACC_eff_up = 0.d0
        return
    end subroutine Acc_reset
    
    subroutine Acc_count(this, toggle)
        class(AccCounter), intent(inout) :: this
        logical, intent(in) :: toggle
        this%NC_eff_up = this%NC_eff_up + 1
        if (toggle) this%ACC_eff_up = this%ACC_eff_up + 1
        return
    end subroutine Acc_count
    
    subroutine Acc_calc_ratio(this)
        class(AccCounter), intent(inout) :: this
        this%acc = this%acc + this%ACC_eff_up / this%NC_eff_up
        return
    end subroutine Acc_calc_ratio

    subroutine Phase_init(this)
        class(PhaseCounter), intent(inout) :: this
        this%count_proposed = 0.d0
        this%count_accepted = 0.d0
        this%sum_abs_arg_proposed = 0.d0
        this%sum_abs_arg_accepted = 0.d0
        this%branch_pi_count = 0.d0
        this%mean_abs_arg_proposed = 0.d0
        this%mean_abs_arg_accepted = 0.d0
        this%branch_pi_fraction = 0.d0
        call this%reset()
        return
    end subroutine Phase_init

    subroutine Phase_reset(this)
        class(PhaseCounter), intent(inout) :: this
        this%NC_prop_bin = 0.d0
        this%NC_acc_bin = 0.d0
        this%sum_abs_prop_bin = 0.d0
        this%sum_abs_acc_bin = 0.d0
        this%branch_pi_count_bin = 0.d0
        return
    end subroutine Phase_reset

    subroutine Phase_count(this, phase, accepted)
        class(PhaseCounter), intent(inout) :: this
        real(kind=8), intent(in) :: phase
        logical, intent(in) :: accepted
        real(kind=8) :: abs_phase

        abs_phase = abs(phase)
        this%NC_prop_bin = this%NC_prop_bin + 1.d0
        this%sum_abs_prop_bin = this%sum_abs_prop_bin + abs_phase
        if (abs_phase > 0.9d0 * PI) this%branch_pi_count_bin = this%branch_pi_count_bin + 1.d0
        if (accepted) then
            this%NC_acc_bin = this%NC_acc_bin + 1.d0
            this%sum_abs_acc_bin = this%sum_abs_acc_bin + abs_phase
        endif
        return
    end subroutine Phase_count

    subroutine Phase_calc_ratio(this)
        class(PhaseCounter), intent(inout) :: this
        this%count_proposed = this%count_proposed + this%NC_prop_bin
        this%count_accepted = this%count_accepted + this%NC_acc_bin
        this%sum_abs_arg_proposed = this%sum_abs_arg_proposed + this%sum_abs_prop_bin
        this%sum_abs_arg_accepted = this%sum_abs_arg_accepted + this%sum_abs_acc_bin
        this%branch_pi_count = this%branch_pi_count + this%branch_pi_count_bin
        if (this%count_proposed > 0.d0) then
            this%mean_abs_arg_proposed = this%sum_abs_arg_proposed / this%count_proposed
            this%branch_pi_fraction = this%branch_pi_count / this%count_proposed
        endif
        if (this%count_accepted > 0.d0) then
            this%mean_abs_arg_accepted = this%sum_abs_arg_accepted / this%count_accepted
        endif
        return
    end subroutine Phase_calc_ratio
end module OperatorHubbard_mod
