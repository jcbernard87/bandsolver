! T1 smoke test: the frozen archival kernel, driven by the oracle wrapper, solves
! Appendix C systems (incl. X/Y endpoint blocks) to small backward error.
! Characterization: with nj=3 and both X and Y nonzero, the legacy BAND(NJ) step
! substitutes dc_1 = E_1 dc_2 + e_1 but omits node 1's X*dc_3 term, so the result is
! not a solution. The library kernel handles that case; the oracle must not be "fixed".
program test_legacy_oracle
    use band_test_utils
    use legacy_oracle
    implicit none
    integer, parameter :: cases(2,4) = reshape([1,3, 2,3, 3,10, 5,50], [2,4])
    double precision, allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), X(:,:), Y(:,:), dc(:,:)
    integer :: c, n, nj, xy
    character(64) :: label
    call seed_rng(12345)
    do c = 1, size(cases, 2)
        do xy = 0, 1
            n = cases(1,c); nj = cases(2,c)
            allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n), dc(n,nj))
            call random_system(n, nj, A, B, D, G, X, Y, xy == 1)
            call legacy_band_solve(n, nj, A, B, D, G, X, Y, dc)
            write(label, '(a,i0,a,i0,a,i0)') 'legacy oracle n=', n, ' nj=', nj, ' xy=', xy
            if (nj == 3 .and. xy == 1) then
                call check(backward_error(n, nj, A, B, D, G, X, Y, dc) > 1.d-8, &
                    trim(label)//' (known legacy limitation: X and Y at nj=3)')
            else
                call check(backward_error(n, nj, A, B, D, G, X, Y, dc) < 1.d-13, trim(label))
            end if
            deallocate(A, B, D, G, X, Y, dc)
        end do
    end do
    if (n_failures > 0) error stop 1
end program test_legacy_oracle
