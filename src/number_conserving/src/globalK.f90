module GlobalK_mod
    use Multiply_mod
    implicit none
    
    public
    private :: ratioK_fermion
    
contains
    real(kind=8) function ratioK_fermion(Gr, phi_new, ii, ntau) result(ratio_fermion)
        use MyMats
! Arguments:
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi_new    
        integer, intent(in) :: ii, ntau
! Local: 
        integer :: P(Norb), nr, nl, no, sign, j
        complex(kind=8), dimension(Norb, Norb) :: Prod, Prodinv, Gr_local, mat_tmp
        complex(kind=8) :: Proddet
        complex(kind=8) :: Vhlp(2, Ndim), Uhlp(Ndim, 2), temp(Ndim, 2), Diff(Ndim, Ndim)
        real(kind=8), dimension(Naux) :: vec_new, vec_old, vec_tmp
        real(kind=8) :: Xdif
        
        vec_new(:) = phi_new(:, ii, ntau)
        vec_old(:) = Conf%phi_list(:, ii, ntau)
        vec_tmp = vec_new - vec_old
        Xdif = sqr_vec(vec_tmp)
        if (Xdif > Zero) then
            do no = 1, Norb
                P(no) = Latt%inv_dim_list(ii, no)
            enddo
            sign = 1
            call Op_U1%get_delta(vec_old, vec_new, sign)
            Prod =  dcmplx(0.d0, 0.d0)
            do nr = 1, Norb
                do nl = 1, Norb
                    Gr_local(nl, nr) = ZKRON(nl, nr) - Gr(P(nl), P(nr))
                enddo
            enddo
            call mmult(mat_tmp, Op_U1%Delta, Gr_local) ! 2*2 matrix multiplication
            do nr = 1, Norb
                do nl = 1, Norb
                    Prod(nl, nr) = ZKRON(nl, nr) + mat_tmp(nl, nr)
                enddo
            enddo
            Proddet = Prod(1,1) * Prod(2,2) - Prod(1,2) * Prod(2,1)
            ratio_fermion = real(Proddet * dconjg(Proddet)) ! output
! Update Green's function
            Prodinv(1,1) = Prod(2,2)
            Prodinv(2,2) = Prod(1,1)
            Prodinv(1,2) = - Prod(1,2)
            Prodinv(2,1) = - Prod(2,1)
            Prodinv = Prodinv / Proddet
            Uhlp = dcmplx(0.d0, 0.d0); Vhlp = dcmplx(0.d0, 0.d0)
            temp = dcmplx(0.d0, 0.d0); Diff = dcmplx(0.d0, 0.d0)
! Vhlp(1:2, 1:Ndim) = Del(1:2) * (1 - Grup)(P(1):P(2), 1:Ndim); Uhlp(1:Ndim, 1:2) = Grup(1:Ndim, P(1):P(2))
            do no = 1, Norb
                Uhlp(:, no) = Gr(:, P(no))
                do j = 1, Ndim
                    Vhlp(no, j) = - Op_U1%Delta(no, 1) * Gr(P(1), j) - Op_U1%Delta(no, 2) * Gr(P(2), j)
                enddo
                Vhlp(no, P(1)) = Vhlp(no, P(1)) + Op_U1%Delta(no, 1)
                Vhlp(no, P(2)) = Vhlp(no, P(2)) + Op_U1%Delta(no, 2)
            enddo
            call mmult(temp, Uhlp, Prodinv)
            call mmult(Diff, temp, Vhlp)
            Gr = Gr - Diff ! output
        else
            ratio_fermion = 1.d0 ! output
        endif
        return
    end function ratioK_fermion
    
    subroutine GlobalK_prop_L(Prop, ratio_fermion, phi_new, nt) ! for shift & wolff flip
        class(Propagator), intent(inout) :: Prop
        real(kind=8), intent(inout) :: ratio_fermion
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi_new
        integer, intent(in) :: nt
        integer :: ii
        do ii = Lq, 1, -1
            ratio_fermion = ratio_fermion * ratioK_fermion(Prop%Gr, phi_new, ii, nt)
        enddo
        call Op_U1%mmult_L(Prop%Gr, Latt, phi_new, nt, 1)
        call Op_U1%mmult_R(Prop%Gr, Latt, phi_new, nt, -1)
        call Op_U1%mmult_L(Prop%UUL, Latt, phi_new, nt, 1)
        return
    end subroutine GlobalK_prop_L
    
    subroutine GlobalK_prop_R(Prop, ratio_fermion, phi_new, nt)
        class(Propagator), intent(inout) :: Prop
        real(kind=8), intent(inout) :: ratio_fermion
        real(kind=8), dimension(Naux, Lq, Ltrot), intent(in) :: phi_new
        integer, intent(in) :: nt
        integer :: ii
        call Op_U1%mmult_R(Prop%Gr, Latt, Conf%phi_list, nt, 1)
        call Op_U1%mmult_L(Prop%Gr, Latt, Conf%phi_list, nt, -1)
        do ii = 1, Lq
            ratio_fermion = ratio_fermion * ratioK_fermion(Prop%Gr, phi_new, ii, nt)
        enddo
        call Op_U1%mmult_R(Prop%UUR, Latt, phi_new, nt, 1)
        return
    end subroutine GlobalK_prop_R
end module GlobalK_mod