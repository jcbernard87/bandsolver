! Shared helpers for Fortran tests: random Appendix C systems and an independent
! residual/backward-error evaluation that applies K directly (no solver involved).
module band_test_utils
    implicit none
    private
    public :: random_system, backward_error, check, n_failures, seed_rng
    integer, save :: n_failures = 0
contains
    subroutine seed_rng(s)
        integer, intent(in) :: s
        integer :: k, i
        integer, allocatable :: seed(:)
        call random_seed(size=k)
        allocate(seed(k))
        seed = [(s + 37*i, i = 1, k)]
        call random_seed(put=seed)
    end subroutine seed_rng

    subroutine random_system(n, nj, A, B, D, G, X, Y, with_xy)
        integer, intent(in) :: n, nj
        logical, intent(in) :: with_xy
        double precision, intent(out) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer :: j, i
        call random_number(A); A = 2*A - 1
        call random_number(B); B = 2*B - 1
        call random_number(D); D = 2*D - 1
        call random_number(G); G = 2*G - 1
        call random_number(X); X = 2*X - 1
        call random_number(Y); Y = 2*Y - 1
        if (.not. with_xy) then
            X = 0; Y = 0
        end if
        do j = 1, nj
            do i = 1, n
                B(i,i,j) = B(i,i,j) + 3.d0*n + 3.d0
            end do
        end do
    end subroutine random_system

    ! ||K dc - G||_inf / (||K||_inf ||dc||_inf + ||G||_inf)
    double precision function backward_error(n, nj, A, B, D, G, X, Y, dc) result(err)
        integer, intent(in) :: n, nj
        double precision, intent(in) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), dc(n,nj)
        double precision :: r(n), rowsum(n), knorm
        integer :: j
        err = 0; knorm = 0
        do j = 1, nj
            r = matmul(B(:,:,j), dc(:,j)) - G(:,j)
            rowsum = sum(abs(B(:,:,j)), dim=2)
            if (j > 1) then
                r = r + matmul(A(:,:,j), dc(:,j-1)); rowsum = rowsum + sum(abs(A(:,:,j)), dim=2)
            end if
            if (j < nj) then
                r = r + matmul(D(:,:,j), dc(:,j+1)); rowsum = rowsum + sum(abs(D(:,:,j)), dim=2)
            end if
            if (j == 1) then
                r = r + matmul(X, dc(:,3)); rowsum = rowsum + sum(abs(X), dim=2)
            end if
            if (j == nj) then
                r = r + matmul(Y, dc(:,nj-2)); rowsum = rowsum + sum(abs(Y), dim=2)
            end if
            err = max(err, maxval(abs(r)))
            knorm = max(knorm, maxval(rowsum))
        end do
        err = err / (knorm*maxval(abs(dc)) + maxval(abs(G)))
    end function backward_error

    subroutine check(cond, msg)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        if (cond) then
            print '(a,a)', 'PASS ', msg
        else
            print '(a,a)', 'FAIL ', msg
            n_failures = n_failures + 1
        end if
    end subroutine check
end module band_test_utils
