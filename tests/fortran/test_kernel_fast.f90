! P1: the column-major "fast" kernel must reproduce the reference loops bit for bit (same
! per-entry accumulation order), for solves and for Newton iterates.
program test_kernel_fast
    use band_test_utils
    use bandsolver_kernel
    implicit none
    integer, parameter :: ns(5) = [1, 2, 3, 7, 12], njs(4) = [3, 4, 50, 300]
    double precision, allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), X(:,:), Y(:,:), d1(:,:), d2(:,:)
    integer :: in, ij, xy, n, nj, s1, s2, f1, f2, mismatches, cases
    double precision :: r1, r2

    call seed_rng(77)
    mismatches = 0; cases = 0
    do in = 1, size(ns)
    do ij = 1, size(njs)
    do xy = 0, 1
        n = ns(in); nj = njs(ij)
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), d1(n,nj), d2(n,nj))
        call random_system(n, nj, A, B, D, G, X, Y, xy == 1)
        call band_solve(n, nj, A, B, D, G, d1, s1, X=X, Y=Y, kernel=KERNEL_REFERENCE, fail_node=f1, min_rel_pivot=r1)
        call band_solve(n, nj, A, B, D, G, d2, s2, X=X, Y=Y, kernel=KERNEL_FAST, fail_node=f2, min_rel_pivot=r2)
        cases = cases + 1
        if (s1 /= BAND_OK .or. s2 /= BAND_OK .or. any(d1 /= d2) .or. r1 /= r2) mismatches = mismatches + 1
        deallocate(A, B, D, G, X, Y, d1, d2)
    end do
    end do
    end do
    print '(a,i0,a,i0)', 'fast vs reference: bitwise mismatches ', mismatches, ' of ', cases
    call check(mismatches == 0, 'fast kernel is bit-identical to the reference loops (sweep incl. X/Y, nj=3)')

    ! Default kernel is fast; singular detection unchanged.
    n = 2; nj = 7
    allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), d1(n,nj), d2(n,nj))
    call random_system(n, nj, A, B, D, G, X, Y, .true.)
    call band_solve(n, nj, A, B, D, G, d1, s1, X=X, Y=Y)
    call band_solve(n, nj, A, B, D, G, d2, s2, X=X, Y=Y, kernel=KERNEL_FAST)
    call check(all(d1 == d2), 'default kernel is the fast kernel')
    A(:,:,4) = 0; B(2,:,4) = 2*B(1,:,4)
    call band_solve(n, nj, A, B, D, G, d1, s1, kernel=KERNEL_FAST, fail_node=f1)
    call band_solve(n, nj, A, B, D, G, d2, s2, kernel=KERNEL_REFERENCE, fail_node=f2)
    call check(s1 == BAND_SINGULAR .and. s2 == BAND_SINGULAR .and. f1 == 4 .and. f2 == 4, 'singular node reported by both kernels')
    call band_solve(n, nj, A, B, D, G, d1, s1, kernel=9)
    call check(s1 == BAND_INVALID_ARGUMENT, 'unknown kernel rejected')
    if (n_failures > 0) error stop 1
end program test_kernel_fast
