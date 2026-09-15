module ObserTau_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserTau
    contains
        procedure :: make   => Obs_tau_make
        procedure :: reset  => Obs_tau_reset
        procedure :: ave    => Obs_tau_ave
        procedure :: calc   => Obs_tau_calc
        final :: Obs_tau_clear
    end type ObserTau
    
contains
    subroutine Obs_tau_make(this)
        class(ObserTau), intent(inout) :: this
        return
    end subroutine Obs_tau_make
    
    subroutine Obs_tau_reset(this)
        class(ObserTau), intent(inout) :: this
        return
    end subroutine Obs_tau_reset
    
    subroutine Obs_tau_ave(this, Nobst)
        class(ObserTau), intent(inout) :: this
        integer, intent(in) :: Nobst
        complex(kind=8) :: znorm
        znorm = dcmplx(1.d0 / dble(Nobst), 0.d0)
        return
    end subroutine Obs_tau_ave
    
    subroutine Obs_tau_clear(this)
        type(ObserTau), intent(inout) :: this
        return
    end subroutine Obs_tau_clear
    
    subroutine Obs_tau_calc(this, PropGr, ntau) ! here ntau ranges from 1 to Ltrot, with ntau=1 denotes equal time correlation
        class(ObserTau), intent(inout) :: this
        class(PropGreen), intent(in) :: PropGr
        integer, intent(in) :: ntau
! Local: 
        ! complex(kind=8), dimension(Ndim, Ndim) :: Gr00up, Gr00upc, Grttup, Grttupc, Gr0tup, Grt0up
        
        ! Gr00up = PropGr%Gr00 ! Gr00(i, j) = <c_i(0) * c^+_j(0)>
        ! Gr00upc = ZKRON - transpose(Gr00up) ! Gr00c(i, j) = <c^+_i(0) * c_j(0)>
        ! Grttup = PropGr%Grtt ! Grtt(i, j) = <c_i(\tau) * c^+_j(\tau)>
        ! Grttupc = ZKRON - transpose(Grttup) ! Grttc(i, j) = <c^+_i(\tau) * c_j(\tau)>
        ! Grt0up = PropGr%Grt0 ! Grt0(i, j) = <c_i(\tau) * c^+_j(0)>
        ! Gr0tup = PropGr%Gr0t ! Gr0t(i, j) = - <c^+_j(\tau) * c_i(0)>

        return
    end subroutine Obs_tau_calc
end module ObserTau_mod