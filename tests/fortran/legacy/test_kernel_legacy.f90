! Library kernel vs the frozen archival BAND/MATINV (built only when the private legacy
! oracle is available; see docs/provenance.md). Legacy pivot mode must match within
! 4 eps (bit-identical on the reference toolchain); partial pivoting within 1e-12.
program test_kernel_legacy
    use band_test_utils
    use legacy_oracle
    use bandsolver_kernel
    implicit none
    integer, parameter :: ns(4) = [1, 2, 5, 12], njs(4) = [3, 4, 50, 500]
    double precision, allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), X(:,:), Y(:,:), dc(:,:), dref(:,:)
    double precision :: ulp_max, rel_partial_max, scale
    integer :: in, ij, xy, n, nj, status

    call seed_rng(2024)
    ulp_max = 0; rel_partial_max = 0
    do in = 1, size(ns)
    do ij = 1, size(njs)
    do xy = 0, 1
        n = ns(in); nj = njs(ij)
        if (nj == 3 .and. xy == 1) cycle   ! known archival defect, see provenance notes
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), dc(n,nj), dref(n,nj))
        call random_system(n, nj, A, B, D, G, X, Y, xy == 1)
        call legacy_band_solve(n, nj, A, B, D, G, X, Y, dref)
        scale = maxval(abs(dref))
        call band_solve(n, nj, A, B, D, G, dc, status, X=X, Y=Y, pivot=PIVOT_LEGACY)
        ulp_max = max(ulp_max, maxval(abs(dc - dref)) / (epsilon(1.d0)*scale))
        call band_solve(n, nj, A, B, D, G, dc, status, X=X, Y=Y, pivot=PIVOT_PARTIAL)
        rel_partial_max = max(rel_partial_max, maxval(abs(dc - dref)) / scale)
        deallocate(A, B, D, G, X, Y, dc, dref)
    end do
    end do
    end do
    print '(a,f8.2)', 'max legacy-mode vs oracle diff [eps*max|dc|]: ', ulp_max
    print '(a,es10.2)', 'max partial vs oracle relative diff: ', rel_partial_max
    call check(ulp_max <= 4.d0, 'legacy pivot mode matches frozen oracle within 4 eps')
    call check(rel_partial_max < 1.d-12, 'partial pivot mode agrees with oracle')

    ! A nearly singular block: the archival MATINV divides by the tiny pivot. Legacy pivot with
    ! the exact singular rule must match it bit for bit; the relative rule stops instead.
    n = 2; nj = 6
    allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), dc(n,nj), dref(n,nj))
    call random_system(n, nj, A, B, D, G, X, Y, .false.)
    X = 0; Y = 0; D(:,:,1) = 0; A(:,:,2) = 0
    B(:,:,1) = reshape([1.d0, 2.d0, 0.5d0, 1.d0 + epsilon(1.d0)], [2, 2])
    call legacy_band_solve(n, nj, A, B, D, G, X, Y, dref)
    call band_solve(n, nj, A, B, D, G, dc, status, X=X, Y=Y, pivot=PIVOT_LEGACY, singular=SINGULAR_EXACT)
    call check(status == BAND_OK .and. all(dc == dref), 'nearly singular block: legacy + exact is bit-identical')
    call band_solve(n, nj, A, B, D, G, dc, status, X=X, Y=Y, pivot=PIVOT_LEGACY)
    call check(status == BAND_SINGULAR, 'nearly singular block: legacy + relative reports singular')
    deallocate(A, B, D, G, X, Y, dc, dref)
    if (n_failures > 0) error stop 1
end program test_kernel_legacy
