module Multiply_mod
    use DQMC_Model_mod
    use ProcessMatrix
    implicit none
    
contains
    subroutine propU_pre(Op_U, Prop, nf, nt)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(in) :: nt, nf
        integer :: ii
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
        enddo
        return
    end subroutine propU_pre
    
    subroutine propU_L(Op_U, Prop, nf, nt)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(in) :: nf, nt
        integer :: ii
        do ii = Ndim, 1, -1
            call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
            call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list(nf, ii, nt), ii, -1)
        enddo
        do ii = Ndim, 1, -1
            call Op_U%mmult_L(Prop%UUL, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
        enddo
        return
    end subroutine propU_L
    
    subroutine propU_R(Op_U, Prop, nf, nt)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(in) :: nf, nt
        integer :: ii
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
            call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list(nf, ii, nt), ii, -1)
        enddo
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
        enddo
        return
    end subroutine propU_R
    
    subroutine propgrU_R(Op_U, Prop, Propgr, nf, nt)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: Propgr
        integer, intent(in) :: nf, nt
        integer :: ii
        do ii = 1, Ndim
            call Op_U%mmult_R(Propgr%Grt0, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
            call Op_U%mmult_L(Propgr%Gr0t, Latt, Conf%phi_list(nf, ii, nt), ii, -1)
            call Op_U%mmult_R(Propgr%Grtt, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
            call Op_U%mmult_L(Propgr%Grtt, Latt, Conf%phi_list(nf, ii, nt), ii, -1)
        enddo
        do ii = 1, Ndim
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, nt), ii, 1)
        enddo
        return
    end subroutine propgrU_R
    
    subroutine propT_pre(Prop)
        class(Propagator), intent(inout) :: Prop
        call Op_T%mmult_R(Prop%UUR, 1)
        return
    end subroutine propT_pre
    
    subroutine propT_L(Prop)
        class(Propagator), intent(inout) :: Prop
        call Op_T%mmult_L(Prop%Gr, 1)
        call Op_T%mmult_R(Prop%Gr, -1)
        call Op_T%mmult_L(Prop%UUL, 1)
        return
    end subroutine propT_L
    
    subroutine propT_R(Prop)
        class(Propagator), intent(inout) :: Prop
        call Op_T%mmult_R(Prop%Gr, 1)
        call Op_T%mmult_L(Prop%Gr, -1)
        call Op_T%mmult_R(Prop%UUR, 1)
        return
    end subroutine propT_R
    
    subroutine propgrT_R(Prop, Propgr)
        class(Propagator), intent(inout) :: Prop
        class(PropGreen), intent(inout) :: Propgr
        call Op_T%mmult_R(Propgr%Grt0, 1)
        call Op_T%mmult_L(Propgr%Gr0t, -1)
        call Op_T%mmult_R(Propgr%Grtt, 1)
        call Op_T%mmult_L(Propgr%Grtt, -1)
        call Op_T%mmult_R(Prop%UUR, 1)
        return
    end subroutine propgrT_R
end module Multiply_mod