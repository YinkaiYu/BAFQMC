module ProcessMatrix
    use MyLattice
    implicit none
    public
    
    type :: Propagator ! allocated in one spin-orbital sector
        complex(kind=8), dimension(:,:), allocatable :: UUR, VUR, UUL, VUL
        complex(kind=8), dimension(:), allocatable :: DUL, DUR
        complex(kind=8), dimension(:,:), allocatable :: Gr
        real(kind=8) :: Xmaxm, Xmeanm
    contains
        procedure :: make => Prop_make
        procedure :: asgn => Prop_assign
        final :: Prop_clear
    end type Propagator
    
    type :: PropGreen
        complex(kind=8), dimension(:,:), allocatable :: Gr00, Grtt, Grt0, Gr0t
        real(kind=8) :: Xmaxm(3), Xmeanm(3)
    contains
        procedure :: make => Propgr_make
        procedure :: reset => Propgr_reset
        final :: Propgr_clear
    end type PropGreen
    
    type, public :: WrapList
        complex(kind=8), dimension(:,:,:), allocatable :: URlist, VRlist, ULlist, VLlist
        complex(kind=8), dimension(:,:), allocatable :: DRlist, DLlist
    contains
        procedure :: make => Wrlist_make
        procedure :: asgn => Wrlist_assign
        final :: Wrlist_clear
    end type WrapList
    
contains
    subroutine Prop_make(this)
        class(Propagator), intent(inout) :: this
        allocate(this%UUR(Ndim, Ndim), this%VUR(Ndim, Ndim), this%UUL(Ndim, Ndim), this%VUL(Ndim, Ndim))
        allocate(this%DUL(Ndim), this%DUR(Ndim))
        allocate(this%Gr(Ndim, Ndim))
        this%UUR = ZKRON; this%VUR = ZKRON; this%DUR = dcmplx(1.d0, 0.d0)
        this%UUL = ZKRON; this%VUL = ZKRON; this%DUL = dcmplx(1.d0, 0.d0)
        this%Gr = ZKRON
        this%Xmaxm = 0.d0; this%Xmeanm = 0.d0
        return
    end subroutine Prop_make
    
    subroutine Prop_assign(this, that)
        class(Propagator), intent(inout) :: this
        class(Propagator), intent(in) :: that
        this%UUL = that%UUL; this%VUL = that%VUL; this%DUL = that%DUL
        this%UUR = that%UUR; this%VUR = that%VUR; this%DUR = that%DUR
        this%Gr = that%Gr
        this%Xmaxm = that%Xmaxm; this%Xmeanm = that%Xmeanm
        return
    end subroutine Prop_assign
    
    subroutine Prop_clear(this)
        type(Propagator), intent(inout) :: this
        deallocate(this%UUR, this%VUR, this%UUL, this%VUL)
        deallocate(this%DUL, this%DUR, this%Gr)
        return
    end subroutine Prop_clear
    
    subroutine Propgr_make(this)
        class(PropGreen), intent(inout) :: this
        allocate(this%Gr00(Ndim, Ndim), this%Grtt(Ndim, Ndim))
        allocate(this%Grt0(Ndim, Ndim), this%Gr0t(Ndim, Ndim))
        this%Gr00 = ZKRON
        this%Grtt = ZKRON
        this%Gr0t = dcmplx(1.d0, 0.d0)
        this%Grt0 = ZKRON
        this%Xmaxm = 0.d0
        this%Xmeanm = 0.d0
        return
    end subroutine Propgr_make

    subroutine Propgr_reset(this, Gr)
        class(PropGreen), intent(inout) :: this
        complex(kind=8), dimension(:,:), intent(in) :: Gr
        this%Gr00 = Gr
        this%Grtt = Gr
        this%Grt0 = Gr
        this%Gr0t = - ( ZKRON - Gr )
        return
    end subroutine Propgr_reset
    
    subroutine Propgr_clear(this)
        type(PropGreen), intent(inout) :: this
        deallocate(this%Gr00, this%Grtt, this%Grt0, this%Gr0t)
        return
    end subroutine Propgr_clear
    
    subroutine Wrlist_make(this)
        class(WrapList), intent(inout) :: this
        allocate(this%URlist(Ndim, Ndim, 0:Nst), this%ULlist(Ndim, Ndim, 0:Nst))
        allocate(this%VRlist(Ndim, Ndim, 0:Nst), this%VLlist(Ndim, Ndim, 0:Nst))
        allocate(this%DRlist(Ndim, 0:Nst), this%DLlist(Ndim, 0:Nst))
        this%URlist = dcmplx(0.d0, 0.d0); this%ULlist = dcmplx(0.d0, 0.d0); this%VRlist = dcmplx(0.d0, 0.d0)
        this%VLlist = dcmplx(0.d0, 0.d0); this%DRlist = dcmplx(0.d0, 0.d0); this%DLlist = dcmplx(0.d0, 0.d0)
        return
    end subroutine Wrlist_make

    subroutine Wrlist_assign(this, that)
        class(WrapList), intent(inout) :: this
        class(WrapList), intent(in) :: that
        this%URlist = that%URlist; this%ULlist = that%ULlist; this%VRlist = that%VRlist
        this%VLlist = that%VLlist; this%DRlist = that%DRlist; this%DLlist = that%DLlist
        return
    end subroutine Wrlist_assign
    
    subroutine Wrlist_clear(this)
        type(WrapList), intent(inout) :: this
        deallocate(this%URlist, this%VRlist, this%ULlist, this%VLlist)
        deallocate(this%DRlist, this%DLlist)
        return
    end subroutine Wrlist_clear
end module ProcessMatrix