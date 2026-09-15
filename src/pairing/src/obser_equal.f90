module ObserEqual_mod
    use ProcessMatrix
    use DQMC_Model_mod
    implicit none
    
    type, public :: ObserEqual
        complex(kind=8), dimension(:,:,:), allocatable  :: den_corr_up, den_corr_do, single_corr
        complex(kind=8), dimension(:), allocatable      :: den_corr_updo
        real(kind=8), dimension(:), allocatable         :: density_site_total
        real(kind=8)                                    :: density_up, density_do, density, density_total
        real(kind=8)                                    :: kinetic, doubleOcc, squareOcc, local_numsquare, pair_equal
        real(kind=8)                                    :: num_up, num_do, numsquare_up, numsquare_do
        real(kind=8)                                    :: onsite_n2_up, onsite_n2_do
        real(kind=8)                                    :: interaction_energy_density, pairing_energy_density
        real(kind=8)                                    :: chemical_energy_density, grand_energy_density, energy_density
        complex(kind=8)                                 :: sf_K, dw_K, psf_Gamma
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
        allocate( this%den_corr_up(Lq, Norb, Norb), &
                  this%den_corr_do(Lq, Norb, Norb), &
                  this%single_corr(Lq, Norb, Norb), &
                  this%den_corr_updo(Lq) )
        allocate( this%density_site_total(Lq) )
        return
    end subroutine Obs_equal_make
    
    subroutine Obs_equal_clear(this)
        type(ObserEqual), intent(inout) :: this
        deallocate( this%den_corr_up, this%den_corr_do, this%single_corr, this%den_corr_updo )
        deallocate( this%density_site_total )
        return
    end subroutine Obs_equal_clear
    
    subroutine Obs_equal_reset(this)
        class(ObserEqual), intent(inout) :: this
        this%den_corr_up   = dcmplx(0.d0,0.d0)
        this%den_corr_do   = dcmplx(0.d0,0.d0)
        this%single_corr   = dcmplx(0.d0,0.d0)
        this%den_corr_updo = dcmplx(0.d0,0.d0)
        this%density_site_total = 0.d0
        this%density_up       = 0.d0
        this%density_do       = 0.d0
        this%density          = 0.d0
        this%density_total    = 0.d0
        this%kinetic          = 0.d0
        this%doubleOcc        = 0.d0
        this%squareOcc        = 0.d0
        this%local_numsquare  = 0.d0
        this%pair_equal       = 0.d0
        this%num_up           = 0.d0
        this%num_do           = 0.d0
        this%numsquare_up     = 0.d0
        this%numsquare_do     = 0.d0
        this%onsite_n2_up     = 0.d0
        this%onsite_n2_do     = 0.d0
        this%interaction_energy_density = 0.d0
        this%pairing_energy_density     = 0.d0
        this%chemical_energy_density    = 0.d0
        this%grand_energy_density       = 0.d0
        this%energy_density             = 0.d0
        this%sf_K        = dcmplx(0.d0,0.d0)
        this%dw_K        = dcmplx(0.d0,0.d0)
        this%psf_Gamma   = dcmplx(0.d0,0.d0)
        return
    end subroutine Obs_equal_reset
    
    subroutine Obs_equal_ave(this, Nobs)
        class(ObserEqual), intent(inout) :: this
        integer, intent(in) :: Nobs
        real(kind=8) :: znorm
        znorm = 1.d0 / dble(Nobs)
        this%den_corr_up = this%den_corr_up * znorm
        this%den_corr_do = this%den_corr_do * znorm
        this%single_corr = this%single_corr * znorm
        this%density_site_total = this%density_site_total * znorm
        this%density_up       = this%density_up       * znorm
        this%density_do       = this%density_do       * znorm
        this%density          = this%density          * znorm
        this%density_total    = this%density_total    * znorm
        this%kinetic          = this%kinetic          * znorm
        this%doubleOcc        = this%doubleOcc        * znorm
        this%squareOcc        = this%squareOcc        * znorm
        this%local_numsquare  = this%local_numsquare  * znorm
        this%pair_equal       = this%pair_equal       * znorm
        this%num_up           = this%num_up           * znorm
        this%num_do           = this%num_do           * znorm
        this%numsquare_up     = this%numsquare_up     * znorm
        this%numsquare_do     = this%numsquare_do     * znorm
        this%den_corr_updo = this%den_corr_updo * znorm
        this%onsite_n2_up     = this%onsite_n2_up     * znorm
        this%onsite_n2_do     = this%onsite_n2_do     * znorm
        this%interaction_energy_density = this%interaction_energy_density * znorm
        this%pairing_energy_density     = this%pairing_energy_density     * znorm
        this%chemical_energy_density    = this%chemical_energy_density    * znorm
        this%grand_energy_density       = this%grand_energy_density       * znorm
        this%energy_density             = this%energy_density             * znorm
        this%sf_K        = this%sf_K        * znorm
        this%dw_K        = this%dw_K        * znorm
        this%psf_Gamma   = this%psf_Gamma   * znorm
        return
    end subroutine Obs_equal_ave
    
    subroutine Obs_equal_calc(this, Prop, ntau)
!   Arguments: 
        class(ObserEqual), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        integer, intent(in) :: ntau
! Local: 
        complex(kind=8), dimension(Nsite, Nsite) :: G_b_bdag, G_b_cdag, G_b_b, G_b_c
        complex(kind=8), dimension(Nsite, Nsite) :: G_c_bdag, G_c_cdag, G_c_b, G_c_c
        complex(kind=8), dimension(Nsite, Nsite) :: G_bdag_bdag, G_bdag_cdag, G_bdag_b, G_bdag_c
        complex(kind=8), dimension(Nsite, Nsite) :: G_cdag_bdag, G_cdag_cdag, G_cdag_b, G_cdag_c
        complex(kind=8) :: Dbb, Dcc, Dbc, Dcb
        complex(kind=8) :: phase, sample_sf_K, sample_dw_K, sample_psf_Gamma
        complex(kind=8) :: nb_i, nc_i
        real(kind=8) :: sample_density_total, sample_doubleOcc, sample_pair_equal
        real(kind=8) :: sample_kinetic, sample_onsite_n2_up, sample_onsite_n2_do
        real(kind=8) :: sample_interaction_energy_density, sample_pairing_energy_density
        real(kind=8) :: sample_chemical_energy_density, sample_energy_density
        real(kind=8) :: kinetic_term, norm_k
        integer :: i, j, no1, no2, ii, jj, imj, nb, no, site, indexK
        logical :: has_K

        ! Gr:
        ! <bm_bP>, <bm_cP>, -<bm_bm>, -<bm_cm>
        ! <cm_bP>, <cm_cP>, -<cm_bm>, -<cm_cm>
        ! <bP_bP>, <bP_cP>, -<bP_bm>, -<bP_cm>
        ! <cP_bP>, <cP_cP>, -<cP_bm>, -<cP_cm>
        G_b_bdag       =   Prop%Gr(1:Nsite, 1:Nsite)
        G_b_cdag       =   Prop%Gr(1:Nsite, Nsite+1:2*Nsite)
        G_b_b          = - Prop%Gr(1:Nsite, 2*Nsite+1:3*Nsite)
        G_b_c          = - Prop%Gr(1:Nsite, 3*Nsite+1:4*Nsite)
        G_c_bdag       =   Prop%Gr(Nsite+1:2*Nsite, 1:Nsite)
        G_c_cdag       =   Prop%Gr(Nsite+1:2*Nsite, Nsite+1:2*Nsite)
        G_c_b          = - Prop%Gr(Nsite+1:2*Nsite, 2*Nsite+1:3*Nsite)
        G_c_c          = - Prop%Gr(Nsite+1:2*Nsite, 3*Nsite+1:4*Nsite)
        G_bdag_bdag    =   Prop%Gr(2*Nsite+1:3*Nsite, 1:Nsite)
        G_bdag_cdag    =   Prop%Gr(2*Nsite+1:3*Nsite, Nsite+1:2*Nsite)
        G_bdag_b       = - Prop%Gr(2*Nsite+1:3*Nsite, 2*Nsite+1:3*Nsite)
        G_bdag_c       = - Prop%Gr(2*Nsite+1:3*Nsite, 3*Nsite+1:4*Nsite)
        G_cdag_bdag    =   Prop%Gr(3*Nsite+1:4*Nsite, 1:Nsite)
        G_cdag_cdag    =   Prop%Gr(3*Nsite+1:4*Nsite, Nsite+1:2*Nsite)
        G_cdag_b       = - Prop%Gr(3*Nsite+1:4*Nsite, 2*Nsite+1:3*Nsite)
        G_cdag_c       = - Prop%Gr(3*Nsite+1:4*Nsite, 3*Nsite+1:4*Nsite)

        sample_density_total = 0.d0
        sample_doubleOcc = 0.d0
        sample_pair_equal = 0.d0
        sample_kinetic = 0.d0
        sample_onsite_n2_up = 0.d0
        sample_onsite_n2_do = 0.d0
        sample_sf_K = dcmplx(0.d0,0.d0)
        sample_dw_K = dcmplx(0.d0,0.d0)
        sample_psf_Gamma = dcmplx(0.d0,0.d0)
        norm_k = 1.d0 / (dble(Lq) * dble(Lq))
        has_K = (mod(Nlx, 3) .eq. 0) .and. (mod(Nly, 3) .eq. 0)
        if (has_K) indexK = Latt%inv_cell_list(2*Nlx/3 + 1, Nly/3 + 1)
        
        do ii = 1, Nsite
            site = Latt%site_list(ii, 1)
            nb_i = G_bdag_b(ii,ii)
            nc_i = G_cdag_c(ii,ii)
            Dbb = nb_i * nb_i + G_bdag_b(ii,ii) * G_b_bdag(ii,ii) + G_bdag_bdag(ii,ii) * G_b_b(ii,ii)
            Dcc = nc_i * nc_i + G_cdag_c(ii,ii) * G_c_cdag(ii,ii) + G_cdag_cdag(ii,ii) * G_c_c(ii,ii)
            Dbc = nb_i * nc_i + G_bdag_cdag(ii,ii) * G_b_c(ii,ii) + G_bdag_c(ii,ii) * G_b_cdag(ii,ii)

            this%density_up = this%density_up + real(nb_i) / dble(Nsite)
            this%density_do = this%density_do + real(nc_i) / dble(Nsite)
            this%density = this%density + real(nb_i + nc_i) / dble(Nsite)
            this%num_up = this%num_up + real(nb_i)
            this%num_do = this%num_do + real(nc_i)
            this%doubleOcc = this%doubleOcc + real(Dbc) / dble(Nsite)
            this%squareOcc = this%squareOcc + real(nb_i * nb_i + nc_i * nc_i) / dble(Nsite)
            this%local_numsquare = this%local_numsquare + real(Dbb + Dcc) / dble(Nsite)
            this%pair_equal = this%pair_equal + real(G_b_c(ii,ii) + G_bdag_cdag(ii,ii)) / dble(Nsite)
            this%density_site_total(site) = this%density_site_total(site) + real(nb_i + nc_i)

            sample_density_total = sample_density_total + real(nb_i + nc_i) / dble(Nsite)
            sample_doubleOcc = sample_doubleOcc + real(Dbc) / dble(Nsite)
            sample_pair_equal = sample_pair_equal + real(G_b_c(ii,ii) + G_bdag_cdag(ii,ii)) / dble(Nsite)
            sample_onsite_n2_up = sample_onsite_n2_up + real(Dbb)
            sample_onsite_n2_do = sample_onsite_n2_do + real(Dcc)
        enddo

        do i = 1, Lq
            do j = 1, Lq
                imj = Latt%imj(i, j)
                do no1 = 1, Norb
                    do no2 = 1, Norb
                        ii = Latt%inv_site_list(i, no1)
                        jj = Latt%inv_site_list(j, no2)
                        Dbb = G_bdag_b(ii,ii) * G_bdag_b(jj,jj) &
                            + G_bdag_b(ii,jj) * G_b_bdag(ii,jj) &
                            + G_bdag_bdag(ii,jj) * G_b_b(ii,jj)
                        Dcc = G_cdag_c(ii,ii) * G_cdag_c(jj,jj) &
                            + G_cdag_c(ii,jj) * G_c_cdag(ii,jj) &
                            + G_cdag_cdag(ii,jj) * G_c_c(ii,jj)
                        Dbc = G_bdag_b(ii,ii) * G_cdag_c(jj,jj) &
                            + G_bdag_cdag(ii,jj) * G_b_c(ii,jj) &
                            + G_bdag_c(ii,jj) * G_b_cdag(ii,jj)
                        Dcb = G_cdag_c(ii,ii) * G_bdag_b(jj,jj) &
                            + G_cdag_bdag(ii,jj) * G_c_b(ii,jj) &
                            + G_cdag_b(ii,jj) * G_c_bdag(ii,jj)

                        this%numsquare_up = this%numsquare_up + real(Dbb)
                        this%numsquare_do = this%numsquare_do + real(Dcc)
                        this%den_corr_up(imj, no1, no2) = this%den_corr_up(imj, no1, no2) + Dbb / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_do(imj, no1, no2) = this%den_corr_do(imj, no1, no2) + Dcc / dcmplx(dble(Lq), 0.d0)
                        this%den_corr_updo(imj) = this%den_corr_updo(imj) + Dbc / dcmplx(dble(Lq), 0.d0)
                        this%single_corr(imj, no1, no2) = this%single_corr(imj, no1, no2) + G_bdag_b(ii,jj) / dcmplx(dble(Lq), 0.d0)

                        if (has_K) then
                            phase = exp(dcmplx(0.d0, Latt%k_dot_r(indexK, imj)))
                            sample_sf_K = sample_sf_K + phase * (G_bdag_b(ii,jj) + G_cdag_c(ii,jj)) * norm_k
                            sample_dw_K = sample_dw_K + phase * (Dbb + Dcc + Dbc + Dcb) * norm_k
                        endif
                        sample_psf_Gamma = sample_psf_Gamma + ( &
                            G_bdag_cdag(ii,ii) * G_c_b(jj,jj) &
                            + G_bdag_c(ii,jj) * G_cdag_b(ii,jj) &
                            + G_bdag_b(ii,jj) * G_cdag_c(ii,jj) ) * norm_k
                    enddo
                enddo
            enddo
        enddo

        do ii = 1, Nsite
            do nb = 1, Nbond
                jj = Latt%L_bonds(ii, nb)
                kinetic_term = RT * real( &
                    G_bdag_b(ii,jj) + G_bdag_b(jj,ii) + &
                    G_cdag_c(ii,jj) + G_cdag_c(jj,ii) ) / dble(Nsite)
                this%kinetic = this%kinetic + kinetic_term
                sample_kinetic = sample_kinetic + kinetic_term
            enddo
        enddo

        sample_interaction_energy_density = &
            (RU1 + RU2) * (sample_onsite_n2_up + sample_onsite_n2_do) / dble(Nsite) &
            + 2.d0 * (RU1 - RU2) * sample_doubleOcc
        sample_pairing_energy_density = RDelta * sample_pair_equal
        sample_chemical_energy_density = -mu * sample_density_total
        sample_energy_density = sample_kinetic + sample_interaction_energy_density &
            + sample_pairing_energy_density

        this%density_total = this%density_total + sample_density_total
        this%onsite_n2_up = this%onsite_n2_up + sample_onsite_n2_up
        this%onsite_n2_do = this%onsite_n2_do + sample_onsite_n2_do
        this%interaction_energy_density = this%interaction_energy_density &
            + sample_interaction_energy_density
        this%pairing_energy_density = this%pairing_energy_density &
            + sample_pairing_energy_density
        this%chemical_energy_density = this%chemical_energy_density &
            + sample_chemical_energy_density
        this%energy_density = this%energy_density + sample_energy_density
        this%grand_energy_density = this%grand_energy_density &
            + sample_energy_density + sample_chemical_energy_density
        this%sf_K = this%sf_K + sample_sf_K
        this%dw_K = this%dw_K + sample_dw_K
        this%psf_Gamma = this%psf_Gamma + sample_psf_Gamma

        return
    end subroutine Obs_equal_calc
end module ObserEqual_mod
