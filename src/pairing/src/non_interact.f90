module NonInteract
    use MyLattice
    implicit none

    public
    private :: def_hamT, exp_general_matrix
    
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
        class(pairingTriangularLattice), intent(in) :: Latt
! Local: 
        complex(kind=8) :: Z, Zc
        integer :: i_site, j_site, nb
        integer :: i1, i2, i3, i4
        integer :: j1, j2, j3, j4
        
        HamT = dcmplx(0.d0, 0.d0)

! nearest bond hopping
        Z  = dcmplx( RT, 0.d0) 
        Zc = dconjg(Z)
        do i_site = 1, Nsite
            do nb = 1, Nbond
                j_site = Latt%L_bonds(i_site, nb)
                i1 = Latt%inv_dim_list(i_site, 1) ! b
                i2 = Latt%inv_dim_list(i_site, 2) ! c
                i3 = Latt%inv_dim_list(i_site, 3) ! b^dagger
                i4 = Latt%inv_dim_list(i_site, 4) ! c^dagger
                j1 = Latt%inv_dim_list(j_site, 1) ! b
                j2 = Latt%inv_dim_list(j_site, 2) ! c
                j3 = Latt%inv_dim_list(j_site, 3) ! b^dagger
                j4 = Latt%inv_dim_list(j_site, 4) ! c^dagger
                HamT(i1,j1) = HamT(i1,j1) + Z
                HamT(i2,j2) = HamT(i2,j2) + Z
                HamT(j3,i3) = HamT(j3,i3) - Z
                HamT(j4,i4) = HamT(j4,i4) - Z
                HamT(j1,i1) = HamT(j1,i1) + dconjg(Z)
                HamT(j2,i2) = HamT(j2,i2) + dconjg(Z)
                HamT(i3,j3) = HamT(i3,j3) - dconjg(Z)
                HamT(i4,j4) = HamT(i4,j4) - dconjg(Z)
            enddo
        enddo

! pairing
        Z  = dcmplx( RDelta, 0.d0) 
        do i_site = 1, Nsite
            i1 = Latt%inv_dim_list(i_site, 1) ! b
            i2 = Latt%inv_dim_list(i_site, 2) ! c
            i3 = Latt%inv_dim_list(i_site, 3) ! b^dagger
            i4 = Latt%inv_dim_list(i_site, 4) ! c^dagger
            HamT(i1,i4) = HamT(i1,i4) + Z
            HamT(i2,i3) = HamT(i2,i3) + Z
            HamT(i3,i2) = HamT(i3,i2) - Z
            HamT(i4,i1) = HamT(i4,i1) - Z
        enddo

! chemical potential
        Z  = dcmplx( - mu, 0.d0) 
        do i_site = 1, Nsite
            i1 = Latt%inv_dim_list(i_site, 1) ! b
            i2 = Latt%inv_dim_list(i_site, 2) ! c
            i3 = Latt%inv_dim_list(i_site, 3) ! b^dagger
            i4 = Latt%inv_dim_list(i_site, 4) ! c^dagger
            HamT(i1,i1) = HamT(i1,i1) + Z
            HamT(i2,i2) = HamT(i2,i2) + Z
            HamT(i3,i3) = HamT(i3,i3) - Z
            HamT(i4,i4) = HamT(i4,i4) - Z
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
        class(pairingTriangularLattice), intent(in) :: Latt
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: HamT
        complex(kind=8), dimension(Ndim) :: WC
        
        call def_hamT(HamT, Latt)
        call exp_general_matrix(HamT, WC, this%expT_P, this%expT_M)

        this%energy_min = minval(real(WC))
        this%energy_max = maxval(real(WC))
        this%bandwidth = this%energy_max - this%energy_min
        return
    end subroutine opT_set

    subroutine exp_general_matrix(HamT, WC, expT_P, expT_M)
        use MyMats
! The bosonic Nambu commutator matrix is non-Hermitian at finite pairing.
! Use the general eigensolver; the Hermitian diagonalizer gives the wrong spectrum.
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(in) :: HamT
        complex(kind=8), dimension(Ndim), intent(out) :: WC
        complex(kind=8), dimension(Ndim, Ndim), intent(out) :: expT_P, expT_M
! Local:
        integer :: info, lwork, nl, nr
        real(kind=8), dimension(2*Ndim) :: rwork
        complex(kind=8), dimension(1,1) :: vl
        complex(kind=8), dimension(Ndim, Ndim) :: work_mat, eigvec, inv_eigvec, temp1, temp2
        complex(kind=8), dimension(1) :: work_query
        complex(kind=8), dimension(:), allocatable :: work
        complex(kind=8) :: det

        work_mat = HamT
        lwork = -1
        call zgeev('N', 'V', Ndim, work_mat, Ndim, WC, vl, 1, eigvec, Ndim, &
                   work_query, lwork, rwork, info)
        if (info .ne. 0) then
            write(6,*) "zgeev workspace query failed in exp_general_matrix, info =", info
            stop
        endif

        lwork = max(1, int(real(work_query(1))))
        allocate(work(lwork))
        work_mat = HamT
        call zgeev('N', 'V', Ndim, work_mat, Ndim, WC, vl, 1, eigvec, Ndim, &
                   work, lwork, rwork, info)
        deallocate(work)
        if (info .ne. 0) then
            write(6,*) "zgeev failed in exp_general_matrix, info =", info
            stop
        endif

        call inv(eigvec, inv_eigvec, det)
        do nr = 1, Ndim
            do nl = 1, Ndim
                temp1(nl, nr) = eigvec(nl, nr) * exp(-Dtau * WC(nr))
                temp2(nl, nr) = eigvec(nl, nr) * exp( Dtau * WC(nr))
            enddo
        enddo
        call mmult(expT_P, temp1, inv_eigvec)
        call mmult(expT_M, temp2, inv_eigvec)
        return
    end subroutine exp_general_matrix
    
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
