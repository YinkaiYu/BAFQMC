module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserEqual
        complex(kind=8), dimension(:,:,:), allocatable  :: den_corr_up, den_corr_do
        complex(kind=8), dimension(:), allocatable      :: den_corr_updo
        real(kind=8), dimension(:), allocatable         :: density_site_total
        real(kind=8)                                    :: density_up,  density_do
        real(kind=8)                                    :: density_total
        real(kind=8)                                    :: kinetic, doubleOcc, squareOcc
        real(kind=8)                                    :: num_up, num_do, numsquare_up, numsquare_do
        real(kind=8)                                    :: onsite_n2_up, onsite_n2_do
        real(kind=8)                                    :: interaction_energy_density, energy_density
        complex(kind=8)                                 :: sf_K, psf_Gamma, dw_K
    contains
        procedure :: make   => Obs_equal_make
        procedure :: reset  => Obs_equal_reset
        procedure :: ave    => Obs_equal_ave
        procedure :: calc   => Obs_equal_calc
        final :: Obs_equal_clear
    end type ObserEqual
    
contains
    subroutine Obs_equal_make(this)
        class(ObserEqual), intent(inout) :: this
        allocate( this%den_corr_up(Lq, Norb, Norb), this%den_corr_do(Lq, Norb, Norb), this%den_corr_updo(Lq) )
        allocate( this%density_site_total(Lq) )
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        deallocate( this%den_corr_up, this%den_corr_do, this%den_corr_updo )
        deallocate( this%density_site_total )
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%den_corr_up   = dcmplx(0.d0,0.d0)
        this%den_corr_do   = dcmplx(0.d0,0.d0)
        this%den_corr_updo = dcmplx(0.d0,0.d0)
        this%density_site_total = 0.d0
        this%density_up  = 0.d0
        this%density_do  = 0.d0
        this%density_total = 0.d0
        this%kinetic     = 0.d0
        this%doubleOcc   = 0.d0
        this%squareOcc   = 0.d0
        this%num_up      = 0.d0
        this%num_do      = 0.d0
        this%numsquare_up = 0.d0
        this%numsquare_do = 0.d0
        this%onsite_n2_up = 0.d0
        this%onsite_n2_do = 0.d0
        this%interaction_energy_density = 0.d0
        this%energy_density = 0.d0
        this%sf_K = dcmplx(0.d0,0.d0)
        this%psf_Gamma = dcmplx(0.d0,0.d0)
        this%dw_K = dcmplx(0.d0,0.d0)
        return
    end subroutine Obs_equal_reset
    
    subroutine Obs_equal_ave(this, Nobs)
        class(ObserEqual), intent(inout) :: this
        integer, intent(in) :: Nobs
        real(kind=8) :: znorm
        znorm = 1.d0 / dble(Nobs)
        this%den_corr_up = this%den_corr_up * znorm
        this%den_corr_do = this%den_corr_do * znorm
        this%density_site_total = this%density_site_total * znorm
        this%density_up  = this%density_up  * znorm
        this%density_do  = this%density_do  * znorm
        this%density_total = this%density_total * znorm
        this%kinetic     = this%kinetic     * znorm
        this%doubleOcc   = this%doubleOcc   * znorm
        this%squareOcc   = this%squareOcc   * znorm
        this%den_corr_updo = this%den_corr_updo * znorm
        this%num_up      = this%num_up * znorm
        this%num_do      = this%num_do * znorm
        this%numsquare_up = this%numsquare_up * znorm
        this%numsquare_do = this%numsquare_do * znorm
        this%onsite_n2_up = this%onsite_n2_up * znorm
        this%onsite_n2_do = this%onsite_n2_do * znorm
        this%interaction_energy_density = this%interaction_energy_density * znorm
        this%energy_density = this%energy_density * znorm
        this%sf_K = this%sf_K * dcmplx(znorm,0.d0)
        this%psf_Gamma = this%psf_Gamma * dcmplx(znorm,0.d0)
        this%dw_K = this%dw_K * dcmplx(znorm,0.d0)
        return
    end subroutine Obs_equal_ave
    
    subroutine Obs_equal_calc(this, Prop, ntau)
!   Arguments: 
        class(ObserEqual), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Ndim, Ndim) :: Grupc, Grup
        complex(kind=8), dimension(Ndim, Ndim) :: Grdoc, Grdo
        complex(kind=8) :: phaseK, sf_K_sample, psf_Gamma_sample, dw_K_sample
        real(kind=8) :: onsite_n2_up_sample, onsite_n2_do_sample
        real(kind=8) :: doubleOcc_sample, kinetic_sample, interaction_energy_sample, density_local
        integer :: i, j, no1, no2, ii, jj, imj, nb, no, indexK
        logical :: hasK
        
        Grup    = Prop%Gr                           !   Gr(i, j)    = <b_i b^+_j >
        Grupc   = transpose(Grup) - ZKRON           !   Grc(i, j)   = <b^+_i b_j > = <b_j b^+_i > - δ(i,j)
        
        Grdo    = dconjg(Prop%Gr)                   !   Gr(i, j)    = <c_i c^+_j >
        Grdoc   = transpose(Grdo) - ZKRON           !   Grc(i, j)   = <c^+_i c_j > = <c_j c^+_i > - δ(i,j)

        onsite_n2_up_sample = 0.d0
        onsite_n2_do_sample = 0.d0
        doubleOcc_sample = 0.d0
        kinetic_sample = 0.d0
        interaction_energy_sample = 0.d0
        sf_K_sample = dcmplx(0.d0,0.d0)
        psf_Gamma_sample = dcmplx(0.d0,0.d0)
        dw_K_sample = dcmplx(0.d0,0.d0)
        hasK = (mod(Nlx,3) == 0) .and. (mod(Nly,3) == 0)
        indexK = 0
        if (hasK) indexK = Latt%inv_cell_list(2*Nlx/3+1, Nly/3+1)
        
        do ii = 1, Ndim
            i = Latt%dim_list(ii, 1)
            density_local = real( Grupc(ii,ii) + Grdoc(ii,ii) )
            this%density_up = this%density_up + real( Grupc(ii,ii) ) / dble(Lq)
            this%density_do = this%density_do + real( Grdoc(ii,ii) ) / dble(Lq)
            this%doubleOcc  = this%doubleOcc  + real( Grupc(ii,ii) * Grdoc(ii,ii) ) / dble(Lq)
            this%squareOcc  = this%squareOcc  + real( Grupc(ii,ii) * Grupc(ii,ii) + Grdoc(ii,ii) * Grdoc(ii,ii) ) / dble(Lq)
            this%num_up = this%num_up + real( Grupc(ii,ii) ) 
            this%num_do = this%num_do + real( Grdoc(ii,ii) ) 
            this%density_total = this%density_total + density_local / dble(Lq)
            this%density_site_total(i) = this%density_site_total(i) + density_local
            onsite_n2_up_sample = onsite_n2_up_sample + real( 2.d0 * Grupc(ii,ii) * Grupc(ii,ii) + Grupc(ii,ii) )
            onsite_n2_do_sample = onsite_n2_do_sample + real( 2.d0 * Grdoc(ii,ii) * Grdoc(ii,ii) + Grdoc(ii,ii) )
            doubleOcc_sample = doubleOcc_sample + real( Grupc(ii,ii) * Grdoc(ii,ii) ) / dble(Lq)
        enddo

        do ii = 1, Ndim
            do jj = 1, Ndim
                this%numsquare_up = this%numsquare_up + real( Grupc(ii,ii) * Grupc(jj,jj) + Grupc(ii,jj) * Grup(ii,jj) )
                this%numsquare_do = this%numsquare_do + real( Grdoc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,jj) * Grdo(ii,jj) )
            enddo
        enddo

        do i = 1, Lq
            do j = 1, Lq
                imj = Latt%imj(i, j)
                ii = Latt%inv_dim_list(i, 1)
                jj = Latt%inv_dim_list(j, 1)
                psf_Gamma_sample = psf_Gamma_sample + Grupc(ii,jj) * Grdoc(ii,jj) / dcmplx(dble(Lq)*dble(Lq), 0.d0)
                if (hasK) then
                    phaseK = exp( dcmplx(0.d0, Latt%k_dot_r(indexK, imj)) )
                    sf_K_sample = sf_K_sample + phaseK * (Grupc(ii,jj) + Grdoc(ii,jj)) / dcmplx(dble(Lq)*dble(Lq), 0.d0)
                    dw_K_sample = dw_K_sample + phaseK * ( &
                        Grupc(ii,ii) * Grupc(jj,jj) + Grupc(ii,jj) * Grup(ii,jj) + &
                        Grdoc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,jj) * Grdo(ii,jj) + &
                        Grupc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,ii) * Grupc(jj,jj) &
                        ) / dcmplx(dble(Lq)*dble(Lq), 0.d0)
                endif
                do no1 = 1, Norb
                    do no2 = 1, Norb
                        ii = Latt%inv_dim_list(i, no1)
                        jj = Latt%inv_dim_list(j, no2)
                        this%den_corr_up(imj, no1, no2) = this%den_corr_up(imj, no1, no2) + ( Grupc(ii,ii) * Grupc(jj,jj) + Grupc(ii,jj) * Grup(ii,jj) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_do(imj, no1, no2) = this%den_corr_do(imj, no1, no2) + ( Grdoc(ii,ii) * Grdoc(jj,jj) + Grdoc(ii,jj) * Grdo(ii,jj) ) / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_updo(imj) = this%den_corr_updo(imj) + ( Grupc(ii,ii) * Grdoc(jj,jj) ) / dcmplx(dble(Lq), 0.d0)
                    enddo
                enddo
            enddo
        enddo

        do ii = 1, Ndim
            do nb = 1, Nbond
                jj = Latt%L_bonds(ii, nb)
                kinetic_sample = kinetic_sample + RT * real( Grupc(ii,jj) + Grupc(jj,ii) + Grdoc(ii,jj) + Grdoc(jj,ii) ) / dble(Lq)
            enddo
        enddo

        ! U1 multiplies (n_b-n_c)^2 and U2 multiplies (n_b+n_c)^2.
        ! Expanding the two squares gives the normal-ordered estimator below.
        interaction_energy_sample = (U1 + U2) * (onsite_n2_up_sample + onsite_n2_do_sample) / dble(Lq) &
            + 2.d0 * (U2 - U1) * doubleOcc_sample
        this%kinetic = this%kinetic + kinetic_sample
        this%onsite_n2_up = this%onsite_n2_up + onsite_n2_up_sample
        this%onsite_n2_do = this%onsite_n2_do + onsite_n2_do_sample
        this%interaction_energy_density = this%interaction_energy_density + interaction_energy_sample
        this%energy_density = this%energy_density + kinetic_sample + interaction_energy_sample
        this%sf_K = this%sf_K + sf_K_sample
        this%psf_Gamma = this%psf_Gamma + psf_Gamma_sample
        this%dw_K = this%dw_K + dw_K_sample

        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod
