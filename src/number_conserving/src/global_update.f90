module GlobalUpdate_mod ! a combination of shift&Wolff update
    use GlobalK_mod
    use Stabilize_mod
    implicit none
    
    public
    private :: spin_reflect, action_dif, Wolff, phi_new
    
    type WolffContainer
        logical, dimension(:), allocatable :: is_marked
        integer, dimension(:), allocatable :: new_adjoined
        integer :: count ! counter of the stack
        real(kind=8) :: size_cluster
    contains
        procedure :: init => Wolff_init
        procedure :: reset => Wolff_reset
        procedure :: flip => Wolff_flip
        final :: Wolff_clear
    end type WolffContainer
    
    type :: GlobalUpdate
        type(Propagator), allocatable, private :: prop
        type(WrapList), allocatable, private :: wrlist
    contains
        procedure :: init => Global_init
        procedure :: clear => Global_clear
        procedure, private :: reset => Global_reset
        procedure, private :: flip => Global_flip
        procedure, private :: sweep_L => Global_sweep_L
        procedure, private :: sweep_R => Global_sweep_R
        procedure :: sweep => Global_sweep
        procedure, nopass :: ctrl_print => Global_control_print
    end type GlobalUpdate
    
    type(AccCounter) :: Acc_U_global
    type(WolffContainer), allocatable :: Wolff
    real(kind=8), dimension(:,:,:), allocatable :: phi_new
    
contains
    subroutine Wolff_init(this)
        class(WolffContainer), intent(inout) :: this
        allocate(this%is_marked(Lq*Ltrot))
        this%is_marked = .false.
        allocate(this%new_adjoined(Lq*Ltrot))
        this%new_adjoined = -1
        this%count = 0
        this%size_cluster = 0.d0
        return
    end subroutine Wolff_init
    
    subroutine Wolff_reset(this)
        class(WolffContainer), intent(inout) :: this
        this%is_marked = .false.
        this%new_adjoined = -1
        this%count = 0
        return
    end subroutine Wolff_reset
    
    subroutine Wolff_clear(this)
        type(WolffContainer), intent(inout) :: this
        deallocate(this%is_marked)
        deallocate(this%new_adjoined)
        return
    end subroutine Wolff_clear
    
    subroutine Wolff_flip(this, iseed, size_cluster)
! Arguments:
        class(WolffContainer), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
! Local:
        real(kind=8), external :: ranf
        real(kind=8) :: theta, xflip, random, ratio, dif
        real(kind=8), dimension(Naux) :: vec_old, vec_new, vec_j
        integer :: iit0, iit, jjt, ii, nti, jj, ntj, n, nf, eff_bond
! randomly choose the initial space-time site iit0 = (ii0, ntau0) of the cluster
        xflip = ranf(iseed)
        iit0 = nranf(iseed, Lq*Ltrot)
! randomly choose the unit vector theta to perform the reflection
        xflip = ranf(iseed)
        theta = xflip * 2 * PI
! initialize the stack: add the first site into the stack
        this%count = 1
        this%new_adjoined(this%count) = iit0
        size_cluster = 1
! The process stops if the stack is exhausted
        if (this%count > 0) then
! choose the first site in new_adjoined, and pop it up from new_adjoined
            iit = this%new_adjoined(1)
            do n = 1, this%count
                this%new_adjoined(n) = this%new_adjoined(n+1)
            enddo
            this%count = this%count - 1
! iit should not be marked for now in principle; otherwise there is an error
            if (.not. this%is_marked(iit)) then
! iit is now being processed. Any site processed should be marked and thus not to be processed again
                this%is_marked(iit) = .true.
! as a part of the cluster, it is now flipped
                ii = Latt%dimt_list(iit, 1)
                nti = Latt%dimt_list(iit, 2)
                vec_old(:) = Conf%phi_list(:, ii, nti)
                vec_new = spin_reflect(vec_old, theta)
                phi_new(:, ii, nti) = vec_new(:) ! output
! approach the adjacent sites if unmarked, and add them into the stack with probability.
                eff_bond = 0
                do nf = 1, 2*Nbond+2
                    jjt = Latt%LT_bonds(iit, nf)
                    if (.not. this%is_marked(jjt)) then
                        eff_bond = eff_bond + 1
                        jj = Latt%dimt_list(jjt, 1)
                        ntj = Latt%dimt_list(jjt, 2)
                        vec_j(:) = Conf%phi_list(:, jj, ntj)
                        if (nf .le. 2*Nbond) then ! space-adjacent
                            dif = action_dif(vec_new, vec_old, vec_j, .true.)
                            ratio = exp(-dif)
                        else ! temporal-adjacent
                            dif = action_dif(vec_new, vec_old, vec_j, .false.)
                            ratio = 1.d0 - exp(-dif)
                        endif
                        random = ranf(iseed)
                        if (ratio .gt. random) then
                            this%count = this%count + 1
                            this%new_adjoined(this%count) = jjt
                            size_cluster = size_cluster + 1
                        endif
                    endif
                enddo
                if (eff_bond == 0) write(6,*) "Wolff: All the adjacent sites of iit=", iit, " have been marked; size_cluster=", size_cluster
            else
                write(6,*) "ERROR: the site popped from the stack has already been marked"; stop
            endif
        endif
        return
    end subroutine Wolff_flip
    
    pure function action_dif(vec_new, vec_old, vec_j, is_space)
        real(kind=8) :: action_dif
        real(kind=8), dimension(Naux), intent(in) :: vec_new, vec_old, vec_j
        logical, intent(in) :: is_space
        real(kind=8), dimension(Naux) :: vec_tmp
        real(kind=8) :: action_new, action_old
! vec_new, vec_old for the marked site iit; vec_j for the unmarked site jjt
        if (is_space) then ! space-adjacent sites
            vec_tmp = vec_new - vec_j
            action_new = sqr_vec(vec_tmp)
            vec_tmp = vec_old - vec_j
            action_old = sqr_vec(vec_tmp)
        else ! temporal-adjacent sites
            vec_tmp = vec_new - vec_j
            vec_tmp = vec_tmp / (Dtau)
            action_new = sqr_vec(vec_tmp)
            vec_tmp = vec_old - vec_j
            vec_tmp = vec_tmp / (Dtau)
            action_old = sqr_vec(vec_tmp)
        endif
        action_dif = (action_new - action_old) * Dtau / 2.0
        return
    end function action_dif
    
    pure function spin_reflect(vec, theta)
        real(kind=8), dimension(Naux) :: spin_reflect
        real(kind=8), dimension(Naux), intent(in) :: vec
        real(kind=8), intent(in) :: theta
        real(kind=8) :: inner_prod
        inner_prod = vec(1)*cos(theta) + vec(2)*sin(theta)
        spin_reflect(1) = vec(1) - 2.0 * cos(theta) * inner_prod
        spin_reflect(2) = vec(2) - 2.0 * sin(theta) * inner_prod
        return
    end function spin_reflect
    
    subroutine Global_init(this)
        class(GlobalUpdate), intent(inout) :: this
        allocate(phi_new(Naux, Lq, Ltrot))
        phi_new = 0.d0
        allocate(this%prop)
        call this%prop%make()
        allocate(this%wrlist)
        call this%wrlist%make()
        call Acc_U_global%init()
        allocate(Wolff)
        call Wolff%init()
        return
    end subroutine Global_init
    
    subroutine Global_clear(this)
        class(GlobalUpdate), intent(inout) :: this
        deallocate(phi_new)
        deallocate(this%prop)
        deallocate(this%wrlist)
        deallocate(Wolff)
        return
    end subroutine Global_clear
    
    subroutine Global_reset(this, Prop, WrList)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(in) :: Prop
        class(WrapList), intent(in) :: WrList
        phi_new = Conf%phi_list
        call this%prop%asgn(Prop)
        call this%wrlist%asgn(WrList)
        call Acc_U_global%reset()
        return
    end subroutine Global_reset
    
    subroutine Global_flip(this, iseed, size_cluster)
        class(GlobalUpdate), intent(inout) :: this
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
! Local: 
        real(kind=8) :: xflip, X
        real(kind=8), external :: ranf
        integer :: ns
        
        call Wolff%reset()
        call Wolff%flip(iseed, size_cluster) ! update phi_new and size_cluster
        do ns = 1, Naux
            xflip = ranf(iseed)
            X = dble( (xflip - 0.5) * abs(shiftGlb(ns)) )
            phi_new(ns, 1:Lq, 1:Ltrot) = phi_new(ns, 1:Lq, 1:Ltrot) + X
        enddo
        return
    end subroutine Global_flip
    
    subroutine Global_sweep_L(this, Prop, WrList, iseed, size_cluster, is_beta)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
        logical, intent(inout) :: is_beta
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: ratio_boson, ratio_fermion, ratio_re, random
        integer :: nt
        
        ratio_fermion = 1.d0
        call this%flip(iseed, size_cluster)
        do nt = Ltrot, 1, -1
            if (mod(nt, Nwrap) == 0) call Wrap_L(this%prop, this%wrlist, nt)
            call propT_L(this%prop)
            if (U1 > Zero) call GlobalK_prop_L(this%prop, ratio_fermion, phi_new, nt)
        enddo
        call Wrap_L(this%prop, this%wrlist, 0)
        ratio_boson = Conf%bosonratio(phi_new)
        ratio_re = ratio_fermion * ratio_boson
        random = ranf(iseed)
        if (ratio_re .ge. random) then
            call Acc_U_global%count(.true.)
            Conf%phi_list = phi_new
            call Prop%asgn(this%prop)
            call WrList%asgn(this%wrlist)
            is_beta = .false.
        else
            call Acc_U_global%count(.false.)
            phi_new = Conf%phi_list
            call this%prop%asgn(Prop)
            call this%wrlist%asgn(WrList)
            is_beta = .true. 
        endif
        return
    end subroutine Global_sweep_L
    
    subroutine Global_sweep_R(this, Prop, WrList, iseed, size_cluster, is_beta)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: iseed
        integer, intent(out) :: size_cluster
        logical, intent(inout) :: is_beta
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: ratio_boson, ratio_fermion, ratio_re, random
        integer :: nt
        
        ratio_fermion = 1.d0
        call this%flip(iseed, size_cluster)
        call Wrap_R(this%prop, this%wrlist, 0)
        do nt = 1, Ltrot
            if (U1 > Zero) call GlobalK_prop_R(this%prop, ratio_fermion, phi_new, nt)
            call propT_R(this%prop)
            if (mod(nt, Nwrap) == 0) call Wrap_R(this%prop, this%wrlist, nt)
        enddo
        ratio_boson = Conf%bosonratio(phi_new)
        ratio_re = ratio_fermion * ratio_boson
        random = ranf(iseed)
        if (ratio_re .ge. random) then
            call Acc_U_global%count(.true.)
            Conf%phi_list = phi_new
            call Prop%asgn(this%prop)
            call WrList%asgn(this%wrlist)
            is_beta = .true.
        else
            call Acc_U_global%count(.false.)
            phi_new = Conf%phi_list
            call this%prop%asgn(Prop)
            call this%wrlist%asgn(WrList)
            is_beta = .false. 
        endif
        return
    end subroutine Global_sweep_R
    
    subroutine Global_sweep(this, Prop, WrList, iseed, is_beta)
        class(GlobalUpdate), intent(inout) :: this
        class(Propagator), intent(inout) :: Prop
        class(WrapList), intent(inout) :: WrList
        integer, intent(inout) :: iseed
        logical, intent(inout) :: is_beta
        integer :: ngl, size_count
        real(kind=8) :: size_cluster
        
        call this%reset(Prop, WrList)
        size_cluster = 0.d0
        do ngl = 1, Nglobal
            if (is_beta) then
                call this%sweep_L(Prop, WrList, iseed, size_count, is_beta)
            else
                call this%sweep_R(Prop, WrList, iseed, size_count, is_beta)
            endif
            size_cluster = size_cluster + dble(size_count)
        enddo
        call Acc_U_global%ratio()
        size_cluster = size_cluster / dble(Nglobal)
        Wolff%size_cluster = Wolff%size_cluster + size_cluster
        return
    end subroutine Global_sweep
    
    subroutine Global_control_print()
        include 'mpif.h'
        real(kind=8) :: collect
        collect = 0.d0
        call MPI_Reduce(Acc_U_global%acc, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Acc_U_global%acc = collect / dble(ISIZE * Nbin)
        collect = 0
        call MPI_Reduce(Wolff%size_cluster, collect, 1, MPI_Real8, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
        if (IRANK == 0) Wolff%size_cluster = collect / dble(ISIZE * Nbin)
        if (IRANK == 0) then
            write(50,*) 'Average size of Wolff cluster                  :', Wolff%size_cluster
            write(50,*) 'Accept_Kglobal                                 :', Acc_U_global%acc
        endif
        return
    end subroutine Global_control_print
end module GlobalUpdate_mod