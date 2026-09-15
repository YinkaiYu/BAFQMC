module LocalU_mod
    use Multiply_mod
    implicit none
    
    public
    private :: LocalU_metro, LocalU_metro_therm, phi_new
    
    
    real(kind=8) :: phi_new
    
contains
    subroutine LocalU_init(Op_U)
        type(OperatorHubbard), intent(inout) :: Op_U
        call Op_U%Acc_U_local%init()
        call Op_U%Acc_U_therm%init()
        call Op_U%Phase_U_local%init()
        return
    end subroutine LocalU_init
    
    subroutine LocalU_clear(Op_U)
        type(OperatorHubbard), intent(inout) :: Op_U
        return
    end subroutine LocalU_clear
    
    subroutine LocalU_reset(Op_U)
        type(OperatorHubbard), intent(inout) :: Op_U
        call Op_U%Acc_U_local%reset()
        call Op_U%Acc_U_therm%reset()
        call Op_U%Phase_U_local%reset()
        return
    end subroutine LocalU_reset
    
    subroutine LocalU_metro(Op_U, Gr, iseed, nf, i_site, ntau)
        use MyMats
! Arguments:
        type(OperatorHubbard), intent(inout) :: Op_U
	    complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        integer, intent(inout) :: iseed
        integer, intent(in) :: i_site, ntau, nf
!   Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: phi_old, phi_new
        real(kind=8) :: xflip, Xdif, random
        real(kind=8) :: ratio_abs, ratio_phase
        complex(kind=8) :: ratio_det, ratio_exp, ratio_total, det_Pblock
        integer :: mm, nn, ii, jj, isec(Nsec), ns, ns1, ns2, ns3
        complex(kind=8), dimension(Ndim, Ndim) :: Gr_new
        complex(kind=8), dimension(Nsec, Nsec) :: Pblock, Pblock_inv
        complex(kind=8) :: tempLeft(Ndim,Nsec), tempLeftMid(Ndim,Nsec), tempRight(Nsec,Ndim), tempMid(Nsec,Nsec)

! Local update on space-time (i_site, ntau) for auxiliary field flavor (nf)
        phi_old = Conf%phi_list(nf, i_site, ntau)
        xflip = ranf(iseed)
        Xdif = dble((xflip - 0.5) * abs(shiftLoc))
        phi_new = phi_old + Xdif
! Calculate Metropolis ratio   
        call Op_U%get_delta(phi_old, phi_new)
        forall(ns=1:Nsec) isec(ns) = Latt%inv_dim_list(i_site, ns)
        do ns1 = 1, Nsec
            do ns2 = 1, Nsec
                Pblock(ns1,ns2) = ZKRON(ns1,ns2) + Op_U%Delta(ns1,ns1) * ( ZKRON(ns1,ns2) - Gr(isec(ns1),isec(ns2)) )
            enddo
        enddo
        call inv_and_det_cmplx(Pblock, Pblock_inv, det_Pblock)
        ! write(6,*) 'arg(Det):', atan2(aimag(det_Pblock), real(det_Pblock)) * 180.0d0 / PI
        ratio_det = det_Pblock**(-0.5d0)
        ratio_exp = Op_U%ratio_gaussian * Op_U%ratio_constant
        ratio_total = ratio_exp * ratio_det
        ratio_phase = atan2(aimag(ratio_total), real(ratio_total))
        ratio_abs = abs(ratio_total)
! Upgrade Green's function and phi
        random = ranf(iseed)
        if (ratio_abs .gt. random) then
            call Op_U%Acc_U_local%count(.true.)
            call Op_U%Phase_U_local%count(ratio_phase, .true.)

            ! ! Gr = Gr - Gr(:,isec) * Pblock_inv * Op_U%Delta * (ZKRON(isec,:) - Gr(isec,:))  ! old version
            ! do ii = 1, Ndim
            !     do ns1 = 1, Nsec
            !         do ns2 = 1, Nsec
            !             do jj = 1, Ndim
            !                 Gr(ii,jj) = Gr(ii,jj) - Gr(ii,isec(ns1)) * Pblock_inv(ns1,ns2) * Op_U%Delta(ns2,ns2) * (ZKRON(isec(ns2),jj) - Gr(isec(ns2),jj))
            !             enddo
            !         enddo
            !     enddo
            ! enddo

            ! ! Gr = Gr - Gr(:,isec) * Pblock_inv * Op_U%Delta * (ZKRON(isec,:) - Gr(isec,:))
            do ns = 1, Nsec
                tempLeft(:,ns) = Gr(:,isec(ns))
            enddo
            ! tempMid = Pblock_inv * Op_U%Delta
            tempMid = Pblock_inv
            do ns = 1, Nsec
                call zscal(Nsec, Op_U%Delta(ns,ns), tempMid(1,ns), 1)
            enddo
            ! tempRight = ZKRON(isec,:) - Gr(isec,:)
            do ns = 1, Nsec
                tempRight(ns,1:Ndim) = ZKRON(isec(ns),1:Ndim) - Gr(isec(ns),1:Ndim)
            enddo
            ! Gr = Gr - tempLeft * tempMid * tempRight
            Gr_new = Gr
            call zgemm('N', 'N', Ndim, Nsec, Nsec, (1.d0,0.d0), tempLeft, Ndim, tempMid, Nsec, (0.d0,0.d0), tempLeftMid, Ndim)
            call zgemm('N', 'N', Ndim, Ndim, Nsec, (-1.d0,0.d0), tempLeftMid, Ndim, tempRight, Nsec, (1.d0,0.d0), Gr_new, Ndim)
            Gr = Gr_new

            Conf%phi_list(nf, i_site, ntau) = phi_new
        else
            call Op_U%Acc_U_local%count(.false.)
            call Op_U%Phase_U_local%count(ratio_phase, .false.)
        endif
        return
    end subroutine LocalU_metro
    
    subroutine LocalU_prop_L(Op_U, Prop, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = Nsite, 1, -1
            call LocalU_metro(Op_U, Prop%Gr, iseed, nf, ii, ntau)
            call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
        enddo
        ! wrap the left
        do ii = Nsite, 1, -1
            call Op_U%mmult_L(Prop%UUL, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
        enddo
        return
    end subroutine LocalU_prop_L
    
    subroutine LocalU_prop_R(Op_U, Prop, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = 1, Nsite
            call Op_U%mmult_R(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
            call Op_U%mmult_L(Prop%Gr, Latt, Conf%phi_list(nf, ii, ntau), ii, -1)
            call LocalU_metro(Op_U, Prop%Gr, iseed, nf, ii, ntau)
        enddo
        ! wrap the right
        do ii = 1, Nsite
            call Op_U%mmult_R(Prop%UUR, Latt, Conf%phi_list(nf, ii, ntau), ii, 1)
        enddo
        return
    end subroutine LocalU_prop_R
    
    subroutine LocalU_metro_therm(Op_U, nf, ii, ntau, iseed)
! Arguments: 
        type(OperatorHubbard), intent(inout) :: Op_U
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau, nf
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: phi_old, phi_new
        real(kind=8) :: xflip, Xdif, random
        real(kind=8) :: ratio_abs
        complex(kind=8) :: ratio_exp

! Local update on space-time (ii, ntau) for auxiliary field flavor (nf)
        phi_old = Conf%phi_list(nf, ii, ntau)
        xflip = ranf(iseed)
        Xdif = dble((xflip - 0.5) * abs(shiftLoc))
        phi_new = phi_old + Xdif
! Calculate auxiliary Gaussian ratio   
        call Op_U%get_delta(phi_old, phi_new)
        ratio_exp = Op_U%ratio_gaussian
        ratio_abs = abs(ratio_exp)
! Upgrade phi
        random = ranf(iseed)
        if (ratio_abs .gt. random) then
            call Op_U%Acc_U_therm%count(.true.)
            Conf%phi_list(nf, ii, ntau) = phi_new
        else
            call Op_U%Acc_U_therm%count(.false.)
        endif
        return
    end subroutine LocalU_metro_therm
    
    subroutine LocalU_prop_therm(Op_U, iseed, nf, ntau)
        type(OperatorHubbard), intent(inout) :: Op_U
        integer, intent(inout) :: iseed
        integer, intent(in) :: ntau, nf
        integer :: ii
        do ii = 1, Nsite
            call LocalU_metro_therm(Op_U, nf, ii, ntau, iseed)
        enddo
        return
    end subroutine LocalU_prop_therm

    subroutine inv_and_det_cmplx(A, invA, det)
        implicit none
        complex(kind=8), intent(in) :: A(:,:)
        complex(kind=8), intent(out) :: invA(size(A,1), size(A,2))
        complex(kind=8), intent(out) :: det
        integer :: n, info, i
        integer, allocatable :: ipiv(:)
        complex(kind=8), allocatable :: work(:)
        integer :: lwork

        n = size(A, 1)
        allocate(ipiv(n))
        invA = A

        call zgetrf(n, n, invA, n, ipiv, info)

        det = (1.0d0, 0.0d0)
        do i = 1, n
            det = det * invA(i, i)
            if (ipiv(i) /= i) det = -det
        end do

        lwork = 4 * n
        allocate(work(lwork))
        call zgetri(n, invA, n, ipiv, work, lwork, info)
        deallocate(work, ipiv)
    end subroutine inv_and_det_cmplx

end module LocalU_mod
