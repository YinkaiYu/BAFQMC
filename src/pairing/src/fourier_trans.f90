module FourierTrans_mod
    use DQMC_Model_mod
    use ObserEqual_mod
    use ObserTau_mod
    implicit none
    
    type, public :: FourierTrans
    contains
        procedure, private, nopass :: m_write_real_1
        procedure, private, nopass :: m_write_real_2
        generic :: write_real => m_write_real_1
        generic :: write_real => m_write_real_2
        
        procedure, private, nopass :: write_cmplx => m_write_cmplx_3
        
        procedure, private, nopass :: m_write_reciprocal_1
        procedure, private, nopass :: m_write_reciprocal_2
        procedure, private, nopass :: m_write_reciprocal_3
        generic :: write_reciprocal => m_write_reciprocal_1
        generic :: write_reciprocal => m_write_reciprocal_2
        generic :: write_reciprocal => m_write_reciprocal_3
        
        procedure, private, nopass :: m_write_k_1
        procedure, private, nopass :: m_write_k_2
        procedure, private, nopass :: m_write_k_3
        generic :: write_k => m_write_k_1
        generic :: write_k => m_write_k_2
        generic :: write_k => m_write_k_3
        
        procedure, private, nopass :: integrate_susc => m_integrate_susc_2_mom
        procedure, private, nopass :: integrate_susc_freq => m_integrate_susc_2_freq
        
        procedure, private, nopass :: m_write_k_tau_2
        procedure, private, nopass :: m_write_k_tau_4
        generic :: write_k_tau => m_write_k_tau_2
        generic :: write_k_tau => m_write_k_tau_4
        
        procedure, private, nopass :: write_r_tau => m_write_r_tau_2
        
        procedure, private, nopass :: write_w => m_write_w
        
        procedure, private :: write_obs_equal => m_write_obs_equal
        procedure, private :: write_obs_tau => m_write_obs_tau
        
        procedure, public :: preq => m_process_obs_equal
        procedure, public :: prtau => m_process_obs_tau
    end type FourierTrans

contains
    subroutine m_write_real_1(gr, filek) ! overloading routine for other correlations
! Arguments:
        real(kind=8), dimension(Lq), intent(in) :: gr
        character(len=*), intent(in) :: filek
! Local: 
        integer :: nr
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            write(20,*) gr(nr)
        enddo
        close(20)
        return
    end subroutine m_write_real_1
    
    subroutine m_write_real_2(gr, filek) ! overloading routine
! Arguments:
        real(kind=8), dimension(:,:), intent(in) :: gr
        character(len=*), intent(in) :: filek
! Local: 
        integer :: nr, nf
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nr = 1, Lq
            write(20,*) Latt%aimj_v(nr, 1), Latt%aimj_v(nr, 2)
            if (size(gr, 1) == Lq) then ! (Lq, Nbond)
                do nf = 1, Nbond
                    write(20,*) gr(nr, nf)
                enddo
            elseif (size(gr, 2) == Lq) then ! (Naux, Lq)
                do nf = 1, Naux
                    write(20,*) gr(nf, nr)
                enddo
            else
                write(6,*) "ERROR: incorrect input size in write_real_2"; stop
            endif
        enddo
        close(20)
        return
    end subroutine m_write_real_2
    
    subroutine m_write_cmplx_3(gr, filek)
        complex(kind=8), dimension(Lq, Norb, Norb), intent(in) :: gr
        character(len=*), intent(in) :: filek
        integer :: nr, no1, no2
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do no2 = 1, Norb
            do no1 = 1, Norb
                do nr = 1, Lq
                    write(20, '(1X,E16.8)', advance='no') real(gr(nr, no1, no2))
                enddo
            enddo
        enddo
        do no2 = 1, Norb
            do no1 = 1, Norb
                do nr = 1, Lq
                    write(20, '(1X,E16.8)', advance='no') aimag(gr(nr, no1, no2))
                enddo
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_cmplx_3
   
    subroutine m_write_reciprocal_1(gk, filek)
        complex(kind=8), dimension(Lq), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            write(20,*) gk(nk)
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_1
    
    subroutine m_write_reciprocal_2(gk, filek)
        complex(kind=8), dimension(Lq, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk, no
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            do no = 1, Nbond
                write(20,*) gk(nk, no)
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_2
    
    subroutine m_write_reciprocal_3(gk, filek)
        complex(kind=8), dimension(Lq, Nbond, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer :: nk, no1, no2
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nk = 1, Lq
            write(20,*) Latt%xk_v(nk, 1), Latt%xk_v(nk, 2)
            do no2 = 1, Nbond
                do no1 = 1, Nbond
                    write(20,*) gk(nk, no1, no2)
                enddo
            enddo
        enddo
        close(20)
        return
    end subroutine m_write_reciprocal_3

    subroutine m_write_k_1(gk, filek, momindex)
        complex(kind=8), dimension(Lq), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30, *) real(gk(momindex)), imag(gk(momindex))
        close(30)
        return
    end subroutine m_write_k_1
    
    subroutine m_write_k_2(gk, filek, momindex, nf)
        complex(kind=8), dimension(Lq, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, nf
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30,*) real(gk(momindex, nf)), imag(gk(momindex, nf))
        close(30)
        return
    end subroutine m_write_k_2
    
    subroutine m_write_k_3(gk, filek, momindex, no1, no2)
        complex(kind=8), dimension(Lq, Nbond, Nbond), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, no1, no2
        open(unit=30, file=filek, status='unknown', action="write", position="append")
        write(30,*) real(gk(momindex, no1, no2)), imag(gk(momindex, no1, no2))
        close(30)
        return
    end subroutine m_write_k_3

    subroutine m_integrate_susc_2_mom(gr, gk)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        complex(kind=8), dimension(Lq), intent(out) :: gk
        integer :: nt, nk
        gk = dcmplx(0.d0, 0.d0)
        do nt = 1, Ltrot
            do nk = 1, Lq
                gk(nk) = gk(nk) + gr(nk, nt)
            enddo
        enddo
        gk = gk * dcmplx(Dtau, 0.d0)
        return
    end subroutine m_integrate_susc_2_mom
    
    subroutine m_integrate_susc_2_freq(gr, gk)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        complex(kind=8), dimension(Ltrot), intent(out) :: gk
        integer :: nt, nw
        gk = dcmplx(0.d0, 0.d0)
        do nw = 1, Ltrot
            do nt = 1, Ltrot
                gk(nw) = gk(nw) + exp( dcmplx(0.d0, 2.d0*Pi*dble((nt-1)*(nw-1))/dble(Ltrot)) ) * gr(1, nt)
            enddo
        enddo
        gk = gk * dcmplx(Dtau, 0.d0)
        return
    end subroutine m_integrate_susc_2_freq
    
    subroutine m_write_r_tau_2(gr, filek, momindex)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gr
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        integer :: nt
        complex(kind=8) :: tmp
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%aimj_v(momindex, 1), Latt%aimj_v(momindex, 2)
        do nt = 1, Ltrot
            tmp = gr(momindex, nt) / dble(Lq)
            write(20,*) real(tmp), imag(tmp)
        enddo
        close(20)
        return
    end subroutine m_write_r_tau_2
    
    subroutine m_write_w(gk, filek)
        complex(kind=8), dimension(Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        real(kind=8) :: omega
        integer :: nw
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        do nw = 1, Ltrot
            omega = 2.d0 * Pi * dble(nw-1) / beta ! Matsubara frequency
            write(20,*) omega
            write(20,*) real(gk(nw)), imag(gk(nw))
        enddo
        close(20)
        return
    end subroutine m_write_w
    
    subroutine m_write_k_tau_2(gk, filek, momindex)
        complex(kind=8), dimension(Lq, Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex
        integer :: nt
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%xk_v(momindex, 1), Latt%xk_v(momindex, 2)
        do nt = 1, Ltrot
            write(20,*) real(gk(momindex, nt)), imag(gk(momindex, nt))
        enddo
        close(20)
        return
    end subroutine m_write_k_tau_2
    
    subroutine m_write_k_tau_4(gk, filek, momindex, no1, no2)
        complex(kind=8), dimension(Lq, Nbond, Nbond, Ltrot), intent(in) :: gk
        character(len=*), intent(in) :: filek
        integer, intent(in) :: momindex, no1, no2
        integer :: nt
        open(unit=20, file=filek, status='unknown', action="write", position="append")
        write(20,*) Latt%xk_v(momindex, 1), Latt%xk_v(momindex, 2)
        do nt = 1, Ltrot
            write(20,*) real(gk(momindex, no1, no2, nt))
        enddo
        close(20)
        return
    end subroutine m_write_k_tau_4
    
    subroutine m_write_obs_equal(this, Obs)
        class(FourierTrans), intent(inout) :: this
        class(ObserEqual), intent(in) :: Obs
        complex(kind=8) :: correlation_up(Lq, Norb, Norb), correlation_do(Lq, Norb, Norb), correlation_updo(Lq)
        character(len=25) :: filek
        integer :: indexzero, no1, no2, site
        
        indexzero = Latt%inv_cell_list(1, 1)
            
        open(unit=80, file='density_up', status='unknown', action="write", position="append")
        write(80,*) Obs%density_up
        close(80)

        open(unit=80, file='density_do', status='unknown', action="write", position="append")
        write(80,*) Obs%density_do
        close(80)

        open(unit=80, file='density', status='unknown', action="write", position="append")
        write(80,*) Obs%density
        close(80)

        open(unit=80, file='density_total', status='unknown', action="write", position="append")
        write(80,*) Obs%density_total
        close(80)

        open(unit=80, file='density_site_total', status='unknown', action="write", position="append")
        do site = 1, Lq
            if (site .lt. Lq) then
                write(80, '(1X,E24.16)', advance='no') Obs%density_site_total(site)
            else
                write(80, '(1X,E24.16)') Obs%density_site_total(site)
            endif
        enddo
        close(80)

        open(unit=80, file='num_up', status='unknown', action="write", position="append")
        write(80,*) Obs%num_up
        close(80)

        open(unit=80, file='num_do', status='unknown', action="write", position="append")
        write(80,*) Obs%num_do
        close(80)

        open(unit=80, file='kinetic', status='unknown', action="write", position="append")
        write(80,*) Obs%kinetic
        close(80)

        open(unit=80, file='doubleOcc', status='unknown', action="write", position="append")
        write(80,*) Obs%doubleOcc
        close(80)

        open(unit=80, file='squareOcc', status='unknown', action="write", position="append")
        write(80,*) Obs%squareOcc
        close(80)

        open(unit=80, file='local_numsquare', status='unknown', action="write", position="append")
        write(80,*) Obs%local_numsquare
        close(80)

        open(unit=80, file='numsquare_up', status='unknown', action="write", position="append")
        write(80,*) Obs%numsquare_up
        close(80)

        open(unit=80, file='numsquare_do', status='unknown', action="write", position="append")
        write(80,*) Obs%numsquare_do
        close(80)

        open(unit=80, file='onsite_n2_up', status='unknown', action="write", position="append")
        write(80,*) Obs%onsite_n2_up
        close(80)

        open(unit=80, file='onsite_n2_do', status='unknown', action="write", position="append")
        write(80,*) Obs%onsite_n2_do
        close(80)

        open(unit=80, file='pair_equal', status='unknown', action="write", position="append")
        write(80,*) Obs%pair_equal
        close(80)

        open(unit=80, file='interaction_energy_density', status='unknown', action="write", position="append")
        write(80,*) Obs%interaction_energy_density
        close(80)

        open(unit=80, file='pairing_energy_density', status='unknown', action="write", position="append")
        write(80,*) Obs%pairing_energy_density
        close(80)

        open(unit=80, file='chemical_energy_density', status='unknown', action="write", position="append")
        write(80,*) Obs%chemical_energy_density
        close(80)

        open(unit=80, file='grand_energy_density', status='unknown', action="write", position="append")
        write(80,*) Obs%grand_energy_density
        close(80)

        open(unit=80, file='energy_density', status='unknown', action="write", position="append")
        write(80,*) Obs%energy_density
        close(80)

        open(unit=80, file='sf_K', status='unknown', action="write", position="append")
        write(80,*) real(Obs%sf_K), aimag(Obs%sf_K)
        close(80)

        open(unit=80, file='dw_K', status='unknown', action="write", position="append")
        write(80,*) real(Obs%dw_K), aimag(Obs%dw_K)
        close(80)

        open(unit=80, file='psf_Gamma', status='unknown', action="write", position="append")
        write(80,*) real(Obs%psf_Gamma), aimag(Obs%psf_Gamma)
        close(80)

        call Fourier_R_to_K(Obs%den_corr_up, correlation_up, Latt)
        call Fourier_R_to_K(Obs%den_corr_do, correlation_do, Latt)
        call Fourier_R_to_K(Obs%den_corr_updo, correlation_updo, Latt)

        do no1 = 1, Norb
            do no2 = 1, Norb
                write(filek, "('den_upup_sub',I0,'',I0)") no1, no2
                call this%write_k(correlation_up, filek, indexzero, no1, no2 )
                write(filek, "('den_dodo_sub',I0,'',I0)") no1, no2
                call this%write_k(correlation_do, filek, indexzero, no1, no2 )
            enddo
        enddo

        filek = 'den_updo'
        call this%write_k(correlation_updo, filek, indexzero )

        ! filek = 'green'
        ! call this%write_cmplx(Obs%single_corr, filek)

        return
    end subroutine m_write_obs_equal

    subroutine m_process_obs_equal(this, Obs)
!#define DEC
        include 'mpif.h'
! Arguments:
        class(FourierTrans), intent(inout) :: this
        class(ObserEqual), intent(inout) :: Obs
! Local:
!        complex(kind=8), dimension(Lq, Nbond, Nbond) :: Collect3
        real(kind=8), dimension(Lq) :: Collect1
        real(kind=8), dimension(Lq, Nbond) :: Collect2
        real(kind=8), dimension(Naux, Lq) :: Collect2prime
        real(kind=8) :: Collect0, Collect1prime(Nbond)
        complex(kind=8), dimension(Lq) :: Collect1cmplx
        complex(kind=8) :: Collect0cmplx
        complex(kind=8) :: Collect3(Lq, Norb, Norb)
        integer :: N
        
        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density_do = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%density_total, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density_total = Collect0/dble(ISIZE)

        Collect1 = 0.d0
        call MPI_REDUCE(Obs%density_site_total, Collect1, Lq, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%density_site_total = Collect1/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%kinetic, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%kinetic = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%doubleOcc, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%doubleOcc = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%squareOcc, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%squareOcc = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%local_numsquare, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%local_numsquare = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%pair_equal, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%pair_equal = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%num_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%num_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%num_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%num_do = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%numsquare_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%numsquare_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%numsquare_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%numsquare_do = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%onsite_n2_up, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%onsite_n2_up = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%onsite_n2_do, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%onsite_n2_do = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%interaction_energy_density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%interaction_energy_density = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%pairing_energy_density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%pairing_energy_density = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%chemical_energy_density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%chemical_energy_density = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%grand_energy_density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%grand_energy_density = Collect0/dble(ISIZE)

        Collect0 = 0.d0
        call MPI_REDUCE(Obs%energy_density, Collect0, 1, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%energy_density = Collect0/dble(ISIZE)

        Collect0cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%sf_K, Collect0cmplx, 1, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%sf_K = Collect0cmplx/dcmplx(dble(ISIZE),0.d0)

        Collect0cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%dw_K, Collect0cmplx, 1, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%dw_K = Collect0cmplx/dcmplx(dble(ISIZE),0.d0)

        Collect0cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%psf_Gamma, Collect0cmplx, 1, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%psf_Gamma = Collect0cmplx/dcmplx(dble(ISIZE),0.d0)

        N = Lq * Norb * Norb

        Collect3 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%den_corr_up, Collect3, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_corr_up = Collect3/dcmplx(dble(ISIZE),0.d0)

        Collect3 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%den_corr_do, Collect3, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_corr_do = Collect3/dcmplx(dble(ISIZE),0.d0)

        Collect3 = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%single_corr, Collect3, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%single_corr = Collect3/dcmplx(dble(ISIZE),0.d0)

        N = Lq

        Collect1cmplx = dcmplx(0.d0, 0.d0)
        call MPI_REDUCE(Obs%den_corr_updo, Collect1cmplx, N, MPI_complex16, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Obs%den_corr_updo = Collect1cmplx/dcmplx(dble(ISIZE),0.d0)

        if (IRANK == 0) call this%write_obs_equal(Obs)
        return
    end subroutine m_process_obs_equal
    
    subroutine m_write_obs_tau(this, Obs)
        class(FourierTrans), intent(inout) :: this
        class(ObserTau), intent(in) :: Obs
        character(len=25) :: filek
        return
    end subroutine m_write_obs_tau
    
    subroutine m_process_obs_tau(this, Obs)
!#define DEC
        include 'mpif.h'
! Arguments:
        class(FourierTrans), intent(inout) :: this
        class(ObserTau), intent(inout) :: Obs
! Local:
        complex(kind=8), dimension(Lq, Ltrot) :: Collect2
!        complex(kind=8), dimension(Lq, Nbond, Nbond, Ltrot) :: Collect4
        integer :: N

        N = Lq * Ltrot

        if (IRANK == 0) call this%write_obs_tau(Obs)
        return
    end subroutine m_process_obs_tau
end module FourierTrans_mod
