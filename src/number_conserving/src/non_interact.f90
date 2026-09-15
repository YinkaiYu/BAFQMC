module NonInteract
    use MyLattice
    implicit none

    public
    private :: def_hamT
    
    type :: OperatorKinetic
        complex(kind=8), dimension(:,:), allocatable :: expT_P, expT_M
        real(kind=8) :: energy_min, energy_max, bandwidth
    contains
        procedure :: make       => opT_make
        procedure :: set        => opT_set
        procedure :: mmult_R    => opT_mmult_R
        procedure :: mmult_L    => opT_mmult_L
        final :: opT_clear
    end type OperatorKinetic
    
contains
    subroutine opT_make(this)
        class(OperatorKinetic), intent(inout) :: this
        allocate(this%expT_P(Ndim, Ndim), this%expT_M(Ndim, Ndim))
        this%expT_P = dcmplx(0.d0, 0.d0); this%expT_M = dcmplx(0.d0, 0.d0)
        return
    end subroutine opT_make
    
    subroutine opT_clear(this)
        type(OperatorKinetic), intent(inout) :: this
        deallocate(this%expT_P, this%expT_M)
        return
    end subroutine opT_clear
    
    subroutine def_hamT(HamT, Latt)
! Arguments: 
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: HamT
        class(triangularLattice), intent(in) :: Latt
! Local: 
        complex(kind=8) :: Z
        integer :: ii, jj, nb
! integer :: i, j
        
        HamT = dcmplx(0.d0, 0.d0)
        Z = dcmplx( RT, 0.d0) 
!  nearest bond hopping
        do ii = 1, Ndim
            do nb = 1, Nbond
                jj = Latt%L_bonds(ii, nb)
                HamT(ii,jj) = HamT(ii,jj) + Z
                HamT(jj,ii) = HamT(jj,ii) + dconjg(Z)
            enddo
        enddo

        do ii = 1, Ndim
            HamT(ii,ii) = HamT(ii,ii) - dcmplx(mu, 0.d0)
        enddo

! write(6,*) 'HamT'
! write(6, "(A)", advance="no") '        '  
! do j = 1, Ndim
!     write(6, "(I17)", advance="no") j  
! end do
! write(6, *)
! do i = 1, Ndim
!     write(6, "(I16)", advance="no") i  
!     do j = 1, Ndim
!         write(6, "(' (',F6.2,',',F6.2,') ')", advance="no") real(HamT(i,j), kind=8), aimag(HamT(i,j))
!     end do
!     write(6, *)
! end do

        return
    end subroutine def_hamT

    subroutine opT_set(this, Latt)
        use MyMats
! Arguments: 
        class(OperatorKinetic), intent(inout) :: this
        class(triangularLattice), intent(in) :: Latt
! Local: 
!        real(kind=8) :: degen, en_free
        integer :: i, nl, nr
        complex(kind=8), dimension(Ndim, Ndim) :: HamT, Hlp1, Hlp1dag, temp1, temp2
        real(kind=8), dimension(Ndim) :: WC
        complex(kind=8), dimension(Ndim) :: dmat1, dmat2
        
        call def_hamT(HamT, Latt)
        call diag(HamT, Hlp1, WC)

        this%energy_min = minval(WC)
        this%energy_max = maxval(WC)
        this%bandwidth = this%energy_max - this%energy_min

        dmat1 = dcmplx(0.d0, 0.d0)
        dmat2 = dcmplx(0.d0, 0.d0)
        forall( i = 1:Ndim )
            dmat1(i) = dcmplx(exp(- Dtau * WC(i)), 0.d0)
            dmat2(i) = dcmplx(exp( Dtau * WC(i)), 0.d0)
        endforall
        Hlp1dag = dconjg(transpose(Hlp1))
        do nr = 1, Ndim
            do nl = 1, Ndim
                temp1(nl, nr) = Hlp1(nl, nr) * dmat1(nr)
                temp2(nl, nr) = Hlp1(nl, nr) * dmat2(nr)
            enddo
        enddo
        call mmult(this%expT_P, temp1, Hlp1dag) ! output
        call mmult(this%expT_M, temp2, Hlp1dag) ! output
        return
    end subroutine opT_set
    
    subroutine opT_mmult_R(this, Mat, nflag)
        use MyMats
!	In Mat Out exp(-Dtau*T) * Mat for nflag = 1
!	In Mat Out exp( Dtau*T) * Mat for nflag = -1  
        class(OperatorKinetic), intent(in) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, intent(in) :: nflag
        complex(kind=8), dimension(Ndim, Ndim) :: temp
        if (nflag == 1) then
            call mmult(temp, this%expT_P, Mat)
        elseif (nflag == -1) then
            call mmult(temp, this%expT_M, Mat)
        else
            write(6,*) "incorrect nflag in opT_mmult_R"; stop
        endif
        Mat = temp
        return
    end subroutine opT_mmult_R
    
    subroutine opT_mmult_L(this, Mat, nflag)
        use MyMats
!	In Mat Out Mat * exp(-Dtau*T) for nflag = 1
!	In Mat Out Mat * exp( Dtau*T) for nflag = -1
        class(OperatorKinetic), intent(in) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        integer, intent(in) :: nflag
        complex(kind=8), dimension(Ndim, Ndim) :: temp
        if (nflag == 1) then
            call mmult(temp, Mat, this%expT_P)
        elseif (nflag == -1) then
            call mmult(temp, Mat, this%expT_M)
        else
            write(6,*) "incorrect nflag in opT_mmult_L"; stop
        endif
        Mat = temp
        return
    end subroutine opT_mmult_L
end module NonInteract
