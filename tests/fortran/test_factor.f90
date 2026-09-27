! P3: Fortran factor/solve must match band_solve (sweep incl. X/Y and nj=3), serve repeated
! right-hand sides, and report singular blocks at factor time.
program test_factor
    use band_test_utils
    use bandsolver_kernel
    use bandsolver_factor
    implicit none
    integer, parameter :: ns(4) = [1, 2, 5, 12], njs(4) = [3, 4, 50, 300]
    double precision, allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), X(:,:), Y(:,:), d1(:,:), d2(:,:)
    type(band_factorization) :: f
    integer :: in, ij, xy, n, nj, st, r
    double precision :: worst
    logical :: ok

    call seed_rng(99)
    worst = 0; ok = .true.
    do in = 1, size(ns)
    do ij = 1, size(njs)
    do xy = 0, 1
        n = ns(in); nj = njs(ij)
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), d1(n,nj), d2(n,nj))
        call random_system(n, nj, A, B, D, G, X, Y, xy == 1)
        call band_solve(n, nj, A, B, D, G, d1, st, X=X, Y=Y)
        call band_factor(n, nj, A, B, D, f, X=X, Y=Y)
        ok = ok .and. f%status == BAND_OK
        call band_factor_solve(f, G, d2, st)
        ok = ok .and. st == BAND_OK
        worst = max(worst, maxval(abs(d2 - d1)) / maxval(abs(d1)))
        deallocate(A, B, D, G, X, Y, d1, d2)
    end do
    end do
    end do
    print '(a,es9.2)', 'max |factor+solve - band_solve| / max|band_solve| = ', worst
    call check(ok, 'factor + solve succeed on the whole sweep')
    call check(worst < 1.d-13, 'factor + solve matches band_solve (incl. X/Y, nj = 3)')

    n = 4; nj = 80
    allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), d1(n,nj), d2(n,nj))
    call random_system(n, nj, A, B, D, G, X, Y, .true.)
    call band_factor(n, nj, A, B, D, f, X=X, Y=Y)
    worst = 0
    do r = 1, 5
        call random_number(G)
        call band_solve(n, nj, A, B, D, G, d1, st, X=X, Y=Y)
        call band_factor_solve(f, G, d2, st)
        worst = max(worst, maxval(abs(d2 - d1)) / maxval(abs(d1)))
    end do
    call check(worst < 1.d-13, 'one factorization serves repeated right-hand sides')
    A(:,:,6) = 0; B(2,:,6) = 2*B(1,:,6)
    call band_factor(n, nj, A, B, D, f, X=X, Y=Y)
    call check(f%status == BAND_SINGULAR .and. f%fail_node == 6, 'singular block reported at factor time (1-based node 6)')
    call band_factor_solve(f, G, d2, st)
    call check(st == BAND_SINGULAR, 'solve with a failed factorization reports it')
    call band_factor(n, 2, A, B, D, f)
    call check(f%status == BAND_INVALID_ARGUMENT, 'nj < 3 rejected')
    if (n_failures > 0) error stop 1
end program test_factor
