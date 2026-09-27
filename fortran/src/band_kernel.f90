!> Newman BAND block solver (Electrochemical Systems, Appendix C), library form.
!>
!> Solves, for dc(n,nj) with nj >= 3,
!>   j = 1:      B_1 dc_1 + D_1 dc_2 + X dc_3                  = G_1
!>   1<j<nj:     A_j dc_{j-1} + B_j dc_j + D_j dc_{j+1}         = G_j
!>   j = nj:     Y dc_{nj-2} + A_nj dc_{nj-1} + B_nj dc_nj      = G_nj
!> by the Appendix C forward elimination dc_j = E_j dc_{j+1} + e_j and back substitution.
!> Loop and accumulation order follow the archival BAND/MATINV so that PIVOT_LEGACY
!> reproduces the historical results; see docs/math.md for the known differences.
module bandsolver_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: BAND_OK = 0
    integer, parameter, public :: BAND_SINGULAR = 1
    integer, parameter, public :: BAND_INVALID_ARGUMENT = 2
    integer, parameter, public :: BAND_NOT_CONVERGED = 3
    integer, parameter, public :: BAND_NON_FINITE = 4
    integer, parameter, public :: BAND_CALLBACK_ERROR = 5

    integer, parameter, public :: PIVOT_PARTIAL = 0  !< row partial pivoting (default)
    integer, parameter, public :: PIVOT_LEGACY = 1   !< archival MATINV pivot heuristic

    public :: band_solve, block_solve

contains

    !> Solve the Appendix C block system. Inputs are not modified.
    !> status: BAND_OK, BAND_SINGULAR (fail_node set), BAND_NON_FINITE, BAND_INVALID_ARGUMENT.
    !> min_rel_pivot: smallest |pivot| / max|block entry| over all node factorizations.
    subroutine band_solve(n, nj, A, B, D, G, dc, status, X, Y, pivot, fail_node, min_rel_pivot)
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj)
        real(dp), intent(out) :: dc(n,nj)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: X(n,n), Y(n,n)
        integer, intent(in), optional :: pivot
        integer, intent(out), optional :: fail_node
        real(dp), intent(out), optional :: min_rel_pivot

        real(dp), allocatable :: E(:,:,:), R(:,:)
        real(dp) :: Xp(n,n), Yw(n,n), Am(n,n), Bm(n,n), Gm(n), rel
        integer :: piv, j, i, k, l, m, np1

        if (present(fail_node)) fail_node = 0
        if (present(min_rel_pivot)) min_rel_pivot = huge(1.0_dp)
        piv = PIVOT_PARTIAL
        if (present(pivot)) piv = pivot
        ! Validate before touching dc: with invalid sizes its extent is not trustworthy.
        if (n < 1 .or. nj < 3 .or. (piv /= PIVOT_PARTIAL .and. piv /= PIVOT_LEGACY)) then
            status = BAND_INVALID_ARGUMENT
            return
        end if
        dc = 0
        Xp = 0; Yw = 0
        if (present(X)) Xp = X
        if (present(Y)) Yw = Y
        if (.not. (all(ieee_is_finite(A)) .and. all(ieee_is_finite(B)) .and. all(ieee_is_finite(D)) &
                   .and. all(ieee_is_finite(G)) .and. all(ieee_is_finite(Xp)) .and. all(ieee_is_finite(Yw)))) then
            status = BAND_NON_FINITE
            return
        end if

        np1 = n + 1
        allocate(E(n,np1,nj), R(n,2*n+1))

        ! Node 1: B_1 [S_D | S_X | s_G] = [D_1 | X | G_1]
        Bm = B(:,:,1)
        R(:,1:n) = D(:,:,1)
        R(:,n+1:2*n) = Xp
        R(:,2*n+1) = G(:,1)
        call block_solve(n, 2*n+1, Bm, R, piv, status, rel)
        if (.not. node_ok(1)) return
        do k = 1, n
            E(k,np1,1) = R(k,2*n+1)
            do l = 1, n
                E(k,l,1) = -R(k,l)
                Xp(k,l) = -R(k,n+l)
            end do
        end do

        do j = 2, nj
            Am = A(:,:,j)
            Bm = B(:,:,j)
            Gm = G(:,j)
            R(:,1:n) = D(:,:,j)
            if (j == 2) then
                ! Fold node 1's second-neighbour term into the upper block: D_2 += A_2 X'
                do i = 1, n
                    do k = 1, n
                        do l = 1, n
                            R(i,k) = R(i,k) + Am(i,l)*Xp(l,k)
                        end do
                    end do
                end do
            end if
            if (j == nj) then
                ! Eliminate Y dc_{nj-2} using dc_{nj-2} = E dc_{nj-1} + e (+ X' dc_3 if nj-2 = 1).
                do i = 1, n
                    do l = 1, n
                        Gm(i) = Gm(i) - Yw(i,l)*E(l,np1,j-2)
                        do m = 1, n
                            Am(i,l) = Am(i,l) + Yw(i,m)*E(m,l,j-2)
                        end do
                    end do
                end do
                if (nj == 3) then
                    ! Not in the archival kernel: node 1's X' dc_3 term lands on dc_nj itself.
                    do i = 1, n
                        do k = 1, n
                            do l = 1, n
                                Bm(i,k) = Bm(i,k) + Yw(i,l)*Xp(l,k)
                            end do
                        end do
                    end do
                end if
            end if
            do i = 1, n
                R(i,np1) = -Gm(i)
                do l = 1, n
                    R(i,np1) = R(i,np1) + Am(i,l)*E(l,np1,j-1)
                    do k = 1, n
                        Bm(i,k) = Bm(i,k) + Am(i,l)*E(l,k,j-1)
                    end do
                end do
            end do
            call block_solve(n, np1, Bm, R(:,1:np1), piv, status, rel)
            if (.not. node_ok(j)) return
            do k = 1, n
                do m = 1, np1
                    E(k,m,j) = -R(k,m)
                end do
            end do
        end do

        ! Back substitution.
        dc(:,nj) = E(:,np1,nj)
        do j = nj - 1, 1, -1
            do k = 1, n
                dc(k,j) = E(k,np1,j)
                do l = 1, n
                    dc(k,j) = dc(k,j) + E(k,l,j)*dc(l,j+1)
                end do
            end do
        end do
        do l = 1, n
            do k = 1, n
                dc(k,1) = dc(k,1) + Xp(k,l)*dc(l,3)
            end do
        end do

        if (.not. all(ieee_is_finite(dc))) then
            status = BAND_NON_FINITE
        else
            status = BAND_OK
        end if

    contains

        logical function node_ok(jnode)
            integer, intent(in) :: jnode
            node_ok = status == BAND_OK
            if (node_ok) then
                if (present(min_rel_pivot)) min_rel_pivot = min(min_rel_pivot, rel)
            else
                if (present(fail_node)) fail_node = jnode
                dc = 0
            end if
        end function node_ok

    end subroutine band_solve

    !> Solve Bm * S = R in place (R becomes S; Bm is destroyed) for an n x n block
    !> with m right-hand sides. status = BAND_SINGULAR when a chosen pivot satisfies
    !> |pivot| <= n*epsilon*max|Bm| (numerically rank deficient). The archival MATINV has
    !> no such threshold; results are identical whenever the threshold is not triggered.
    !> rel = min |pivot| / max|Bm| (1 for a perfectly scaled diagonal block).
    subroutine block_solve(n, m, Bm, R, pivot, status, rel)
        integer, intent(in) :: n, m, pivot
        real(dp), intent(inout) :: Bm(n,n), R(n,m)
        integer, intent(out) :: status
        real(dp), intent(out) :: rel
        if (pivot == PIVOT_LEGACY) then
            call solve_legacy(n, m, Bm, R, status, rel)
        else
            call solve_partial(n, m, Bm, R, status, rel)
        end if
    end subroutine block_solve

    !> Gauss-Jordan elimination with row partial pivoting.
    subroutine solve_partial(n, m, Bm, R, status, rel)
        integer, intent(in) :: n, m
        real(dp), intent(inout) :: Bm(n,n), R(n,m)
        integer, intent(out) :: status
        real(dp), intent(out) :: rel
        real(dp) :: bscale, f, rowB(n), rowR(m)
        integer :: k, p, i

        bscale = maxval(abs(Bm))
        rel = huge(1.0_dp)
        status = BAND_SINGULAR
        if (bscale == 0) return
        do k = 1, n
            p = k - 1 + maxloc(abs(Bm(k:n,k)), dim=1)
            if (abs(Bm(p,k)) <= n*epsilon(1.0_dp)*bscale) return
            rel = min(rel, abs(Bm(p,k))/bscale)
            if (p /= k) then
                rowB = Bm(k,:); Bm(k,:) = Bm(p,:); Bm(p,:) = rowB
                rowR = R(k,:);  R(k,:) = R(p,:);   R(p,:) = rowR
            end if
            f = 1.0_dp/Bm(k,k)
            Bm(k,k:n) = Bm(k,k:n)*f
            R(k,:) = R(k,:)*f
            do i = 1, n
                if (i == k) cycle
                f = Bm(i,k)
                if (f == 0) cycle
                Bm(i,k:n) = Bm(i,k:n) - f*Bm(k,k:n)
                R(i,:) = R(i,:) - f*R(k,:)
            end do
        end do
        status = BAND_OK
    end subroutine solve_partial

    !> Archival MATINV (Appendix C.4): for each unused row, find its largest (btry) and
    !> second-largest (bnext) entries among unused columns; pivot on the row with the
    !> smallest bnext/btry ratio, swap that row into the pivot column's position, and
    !> Gauss-Jordan eliminate. Operation order is kept for bitwise agreement.
    !> Differences: reports BAND_SINGULAR (instead of printing DETERM=0 and continuing)
    !> when no nonzero entry remains among unused rows/columns or the pivot is below
    !> the relative threshold documented in block_solve.
    subroutine solve_legacy(n, m, Bm, R, status, rel)
        integer, intent(in) :: n, m
        real(dp), intent(inout) :: Bm(n,n), R(n,m)
        integer, intent(out) :: status
        real(dp), intent(out) :: rel
        ! 1.1 is a default-real literal in the archival source; keep its exact value.
        real(dp), parameter :: bmax0 = real(1.1, dp)
        logical :: used(n), found
        real(dp) :: bmax, bnext, btry, f, save, bscale
        integer :: nn, i, j, k, jc, irow, jcol

        bscale = maxval(abs(Bm))
        rel = huge(1.0_dp)
        status = BAND_SINGULAR
        used = .false.
        irow = 0; jcol = 0
        do nn = 1, n
            bmax = bmax0
            found = .false.
            do i = 1, n
                if (used(i)) cycle
                bnext = 0
                btry = 0
                do j = 1, n
                    if (used(j)) cycle
                    if (abs(Bm(i,j)) <= bnext) cycle
                    bnext = abs(Bm(i,j))
                    if (bnext <= btry) cycle
                    bnext = btry
                    btry = abs(Bm(i,j))
                    jc = j
                    found = .true.
                end do
                if (bnext >= bmax*btry) cycle
                bmax = bnext/btry
                irow = i
                jcol = jc
            end do
            if (.not. found) return
            used(jcol) = .true.
            if (jcol /= irow) then
                do j = 1, n
                    save = Bm(irow,j); Bm(irow,j) = Bm(jcol,j); Bm(jcol,j) = save
                end do
                do k = 1, m
                    save = R(irow,k); R(irow,k) = R(jcol,k); R(jcol,k) = save
                end do
            end if
            if (abs(Bm(jcol,jcol)) <= n*epsilon(1.0_dp)*bscale) return
            rel = min(rel, abs(Bm(jcol,jcol))/bscale)
            f = 1.0_dp/Bm(jcol,jcol)
            do j = 1, n
                Bm(jcol,j) = Bm(jcol,j)*f
            end do
            do k = 1, m
                R(jcol,k) = R(jcol,k)*f
            end do
            do i = 1, n
                if (i == jcol) cycle
                f = Bm(i,jcol)
                do j = 1, n
                    Bm(i,j) = Bm(i,j) - f*Bm(jcol,j)
                end do
                do k = 1, m
                    R(i,k) = R(i,k) - f*R(jcol,k)
                end do
            end do
        end do
        status = BAND_OK
    end subroutine solve_legacy

end module bandsolver_kernel
