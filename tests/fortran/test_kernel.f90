! T2: Fortran kernel tests — backward error and failure paths. Agreement with the
! archival kernel is tested separately in legacy/test_kernel_legacy.f90 (optional).
program test_kernel
    use band_test_utils
    use bandsolver_kernel
    implicit none
    integer, parameter :: ns(4) = [1, 2, 5, 12], njs(4) = [3, 4, 50, 500]
    double precision, allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), X(:,:), Y(:,:)
    double precision, allocatable :: dc(:,:), dref(:,:), A0(:,:,:), B0(:,:,:), G0(:,:)
    double precision :: berr_max(0:1)
    integer :: in, ij, xy, piv, n, nj, status, node, fnode
    character(80) :: label

    call seed_rng(2024)
    berr_max = 0
    do in = 1, size(ns)
    do ij = 1, size(njs)
    do xy = 0, 1
        n = ns(in); nj = njs(ij)
        call alloc(n, nj)
        call random_system(n, nj, A, B, D, G, X, Y, xy == 1)
        A0 = A; B0 = B; G0 = G
        do piv = PIVOT_PARTIAL, PIVOT_LEGACY
            call band_solve(n, nj, A, B, D, G, dc, status, X=X, Y=Y, pivot=piv)
            write(label, '(a,i0,a,i0,a,i0,a,i0)') 'solve n=', n, ' nj=', nj, ' xy=', xy, ' pivot=', piv
            call check(status == BAND_OK, trim(label)//' status')
            berr_max(piv) = max(berr_max(piv), backward_error(n, nj, A, B, D, G, X, Y, dc))
        end do
        call check(all(A == A0) .and. all(B == B0) .and. all(G == G0), trim(label)//' inputs unchanged')
        call dealloc()
    end do
    end do
    end do
    print '(a,2es10.2)', 'max backward error (partial, legacy): ', berr_max
    call check(maxval(berr_max) < 1.d-13, 'backward error sweep < 1e-13 (incl. nj=3 with X and Y)')

    ! Singular pivot block detection with the failing node reported, both pivot modes.
    do piv = PIVOT_PARTIAL, PIVOT_LEGACY
        do node = 1, 3
            n = 2; nj = 7
            call alloc(n, nj)
            call random_system(n, nj, A, B, D, G, X, Y, .false.)
            select case (node)
            case (1); ij = 1
            case (2); ij = 4
            case default; ij = nj
            end select
            A(:,:,ij) = 0; B(2,:,ij) = 2*B(1,:,ij)   ! Bhat_j = B_j, rank 1
            call band_solve(n, nj, A, B, D, G, dc, status, X=X, Y=Y, pivot=piv, fail_node=fnode)
            write(label, '(a,i0,a,i0)') 'singular block at node ', ij, ' pivot=', piv
            call check(status == BAND_SINGULAR .and. fnode == ij, trim(label))
            call dealloc()
        end do
    end do

    ! Exactly singular 1x1 block.
    call alloc(1, 5)
    call random_system(1, 5, A, B, D, G, X, Y, .false.)
    A(:,:,3) = 0; B(:,:,3) = 0
    call band_solve(1, 5, A, B, D, G, dc, status, fail_node=fnode)
    call check(status == BAND_SINGULAR .and. fnode == 3, 'zero 1x1 block reported singular')
    call dealloc()

    ! Blocks that need an off-diagonal pivot: B = [[0,1],[1,0]].
    do piv = PIVOT_PARTIAL, PIVOT_LEGACY
        call alloc(2, 6)
        call random_system(2, 6, A, B, D, G, X, Y, .true.)
        A = 0.1d0*A; D = 0.1d0*D; X = 0.1d0*X; Y = 0.1d0*Y
        B = 0; B(1,2,:) = 1; B(2,1,:) = 1
        call band_solve(2, 6, A, B, D, G, dc, status, X=X, Y=Y, pivot=piv)
        call check(status == BAND_OK .and. backward_error(2, 6, A, B, D, G, X, Y, dc) < 1.d-14, &
            'zero-diagonal blocks solve (off-diagonal pivot required)')
        call dealloc()
    end do

    ! Non-finite input and invalid sizes.
    call alloc(2, 5)
    call random_system(2, 5, A, B, D, G, X, Y, .false.)
    G(1,2) = ieee_nan()
    call band_solve(2, 5, A, B, D, G, dc, status)
    call check(status == BAND_NON_FINITE, 'NaN input reported non-finite')
    call dealloc()
    call alloc(2, 2)
    call band_solve(2, 2, A, B, D, G, dc, status)
    call check(status == BAND_INVALID_ARGUMENT, 'nj < 3 rejected')
    call band_solve(2, 5, A, B, D, G, dc, status, pivot=7)
    call check(status == BAND_INVALID_ARGUMENT, 'unknown pivot rejected')
    call dealloc()

    if (n_failures > 0) error stop 1
contains
    subroutine alloc(n, nj)
        integer, intent(in) :: n, nj
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), dc(n,nj), dref(n,nj))
    end subroutine alloc
    subroutine dealloc()
        deallocate(A, B, D, G, X, Y, dc, dref)
    end subroutine dealloc
    double precision function ieee_nan()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        ieee_nan = ieee_value(0.d0, ieee_quiet_nan)
    end function ieee_nan
end program test_kernel
