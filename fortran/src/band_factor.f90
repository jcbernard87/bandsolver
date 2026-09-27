!> Factor once, solve many: the BAND elimination split into a factorization of the block
!> matrix (A, B, D, X, Y) and a cheap solve for any right-hand side G. Mirrors
!> cpp/src/factor.cpp. Per node it stores the partially pivoted LU of the reduced pivot
!> block B^_j = B_j + A'_j E_{j-1}, the elimination blocks E_j = -B^_j^{-1} D'_j and the
!> effective lower blocks A'_j (with the Y correction at the last node). A solve is a forward
!> sweep plus back substitution, O(nj n^2) instead of O(nj n^3). Loops run down columns.
module bandsolver_factor
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use bandsolver_kernel, only: BAND_OK, BAND_SINGULAR, BAND_INVALID_ARGUMENT, BAND_NON_FINITE
    implicit none
    private

    public :: band_factorization, band_factor, band_factor_solve

    type :: band_factorization
        integer :: n = 0, nj = 0
        integer :: status = BAND_INVALID_ARGUMENT
        integer :: fail_node = 0                 !< 1-based node of a singular block, else 0
        real(dp) :: min_rel_pivot = 0
        real(dp), allocatable :: lu(:,:,:)       !< (n,n,nj) LU of B^_j, unit lower L
        integer, allocatable :: perm(:,:)        !< (n,nj) row permutation of each LU
        real(dp), allocatable :: E(:,:,:)        !< (n,n,nj)
        real(dp), allocatable :: Aeff(:,:,:)     !< (n,n,nj)
        real(dp), allocatable :: Xp(:,:), Y(:,:) !< (n,n)
    end type band_factorization

contains

    !> In-place LU with row partial pivoting; same singularity rule as band_solve.
    subroutine lu_factor(n, a, perm, ok, rel)
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n,n)
        integer, intent(out) :: perm(n)
        logical, intent(out) :: ok
        real(dp), intent(out) :: rel
        real(dp) :: scale, t, inv
        integer :: i, k, p, c, itmp
        scale = 0
        do c = 1, n
            do i = 1, n
                scale = max(scale, abs(a(i,c)))
            end do
        end do
        rel = huge(1.0_dp)
        perm = [(i, i = 1, n)]
        ok = .false.
        if (scale == 0) return
        do k = 1, n
            p = k
            do i = k + 1, n
                if (abs(a(i,k)) > abs(a(p,k))) p = i
            end do
            if (abs(a(p,k)) <= n*epsilon(1.0_dp)*scale) return
            rel = min(rel, abs(a(p,k))/scale)
            if (p /= k) then
                do c = 1, n
                    t = a(p,c); a(p,c) = a(k,c); a(k,c) = t
                end do
                itmp = perm(p); perm(p) = perm(k); perm(k) = itmp
            end if
            inv = 1.0_dp/a(k,k)
            a(k+1:n,k) = a(k+1:n,k)*inv
            do c = k + 1, n
                a(k+1:n,c) = a(k+1:n,c) - a(k+1:n,k)*a(k,c)
            end do
        end do
        ok = .true.
    end subroutine lu_factor

    !> Solve (LU) x = b for one vector (permutation applied to b).
    subroutine lu_solve(n, lu, perm, b, x)
        integer, intent(in) :: n, perm(n)
        real(dp), intent(in) :: lu(n,n), b(n)
        real(dp), intent(out) :: x(n)
        integer :: i, k
        do i = 1, n
            x(i) = b(perm(i))
        end do
        do k = 1, n - 1                         ! forward, column-oriented
            x(k+1:n) = x(k+1:n) - lu(k+1:n,k)*x(k)
        end do
        do k = n, 1, -1                         ! backward, column-oriented
            x(k) = x(k)/lu(k,k)
            x(1:k-1) = x(1:k-1) - lu(1:k-1,k)*x(k)
        end do
    end subroutine lu_solve

    subroutine band_factor(n, nj, A, B, D, f, X, Y)
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: A(n,n,nj), B(n,n,nj), D(n,n,nj)
        type(band_factorization), intent(out) :: f
        real(dp), intent(in), optional :: X(n,n), Y(n,n)
        real(dp) :: Dm(n,n), col(n), rel
        logical :: ok
        integer :: j, c

        if (n < 1 .or. nj < 3) return                           ! status = invalid argument
        if (.not. (all(ieee_is_finite(A)) .and. all(ieee_is_finite(B)) .and. all(ieee_is_finite(D)))) then
            f%status = BAND_NON_FINITE
            return
        end if
        f%n = n; f%nj = nj
        allocate(f%lu(n,n,nj), f%perm(n,nj), f%E(n,n,nj), f%Aeff(n,n,nj), f%Xp(n,n), f%Y(n,n))
        f%E = 0; f%Aeff = 0; f%Xp = 0; f%Y = 0
        if (present(Y)) f%Y = Y
        if (present(X)) then
            if (.not. all(ieee_is_finite(X))) then
                f%status = BAND_NON_FINITE; return
            end if
        end if
        if (.not. all(ieee_is_finite(f%Y))) then
            f%status = BAND_NON_FINITE; return
        end if
        f%min_rel_pivot = huge(1.0_dp)

        f%lu(:,:,1) = B(:,:,1)
        call lu_factor(n, f%lu(:,:,1), f%perm(:,1), ok, rel)
        if (.not. ok) then
            f%status = BAND_SINGULAR; f%fail_node = 1; return
        end if
        f%min_rel_pivot = min(f%min_rel_pivot, rel)
        do c = 1, n
            call lu_solve(n, f%lu(:,:,1), f%perm(:,1), D(:,c,1), col)
            f%E(:,c,1) = -col
            if (present(X)) then
                call lu_solve(n, f%lu(:,:,1), f%perm(:,1), X(:,c), col)
                f%Xp(:,c) = -col
            end if
        end do

        do j = 2, nj
            f%Aeff(:,:,j) = A(:,:,j)
            f%lu(:,:,j) = B(:,:,j)
            Dm = D(:,:,j)
            if (j == 2) call gemm_add(n, f%Aeff(:,:,j), f%Xp, Dm)                 ! D_2 += A_2 X'
            if (j == nj) then
                call gemm_add(n, f%Y, f%E(:,:,nj-2), f%Aeff(:,:,j))               ! A += Y E_{nj-2}
                if (nj == 3) call gemm_add(n, f%Y, f%Xp, f%lu(:,:,j))           ! B += Y X'
            end if
            call gemm_add(n, f%Aeff(:,:,j), f%E(:,:,j-1), f%lu(:,:,j))           ! B^ = B + A' E
            call lu_factor(n, f%lu(:,:,j), f%perm(:,j), ok, rel)
            if (.not. ok) then
                f%status = BAND_SINGULAR; f%fail_node = j; return
            end if
            f%min_rel_pivot = min(f%min_rel_pivot, rel)
            if (j < nj) then
                do c = 1, n
                    call lu_solve(n, f%lu(:,:,j), f%perm(:,j), Dm(:,c), col)
                    f%E(:,c,j) = -col
                end do
            end if
        end do
        f%status = BAND_OK
    end subroutine band_factor

    !> C += A B (column-major, inner loop down columns).
    pure subroutine gemm_add(n, A, B, C)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n), B(n,n)
        real(dp), intent(inout) :: C(n,n)
        integer :: k, l
        do k = 1, n
            do l = 1, n
                C(:,k) = C(:,k) + A(:,l)*B(l,k)
            end do
        end do
    end subroutine gemm_add

    subroutine band_factor_solve(f, G, dc, status)
        type(band_factorization), intent(in) :: f
        real(dp), intent(in) :: G(:,:)
        real(dp), intent(out) :: dc(:,:)
        integer, intent(out) :: status
        real(dp) :: rhs(f%n)
        integer :: n, nj, j, l

        status = f%status
        if (status /= BAND_OK) return
        n = f%n; nj = f%nj
        if (size(G, 1) /= n .or. size(G, 2) /= nj .or. size(dc, 1) /= n .or. size(dc, 2) /= nj) then
            status = BAND_INVALID_ARGUMENT; return
        end if
        if (.not. all(ieee_is_finite(G))) then
            dc = 0; status = BAND_NON_FINITE; return
        end if
        call lu_solve(n, f%lu(:,:,1), f%perm(:,1), G(:,1), dc(:,1))
        do j = 2, nj
            rhs = G(:,j)
            if (j == nj) then
                do l = 1, n
                    rhs = rhs - f%Y(:,l)*dc(l,nj-2)
                end do
            end if
            do l = 1, n
                rhs = rhs - f%Aeff(:,l,j)*dc(l,j-1)
            end do
            call lu_solve(n, f%lu(:,:,j), f%perm(:,j), rhs, dc(:,j))
        end do
        do j = nj - 1, 1, -1
            do l = 1, n
                dc(:,j) = dc(:,j) + f%E(:,l,j)*dc(l,j+1)
            end do
        end do
        do l = 1, n
            dc(:,1) = dc(:,1) + f%Xp(:,l)*dc(l,3)
        end do
        status = merge(BAND_OK, BAND_NON_FINITE, all(ieee_is_finite(dc)))
    end subroutine band_factor_solve

end module bandsolver_factor
