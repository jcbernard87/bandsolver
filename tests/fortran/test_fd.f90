! F2: Fortran finite-difference Jacobians. Same nonlinear test problem as tests/cpp/test_fd.cpp:
! F(i,j) = sum_k Bm c(:,j) + Am c(:,j-1) + Dm c(:,j+1) + 0.5 sin(c(i,j)) c(i+1 mod n, j) - b(i) j
!          + [j=1] sum_k Xm(i,k) c(k,3)^2 + [j=nj] sum_k Ym(i,k) c(k,nj-2)^3
module fd_test_problem
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bandsolver_newton, only: band_problem
    use bandsolver_fd, only: band_residual_problem
    implicit none
    type, extends(band_residual_problem) :: nl_residual
        integer :: n = 0, nj = 0
        real(dp), allocatable :: Am(:,:), Bm(:,:), Dm(:,:), Xm(:,:), Ym(:,:), b(:)
        logical :: fail = .false.
    contains
        procedure :: residual => nl_residual_eval
    end type nl_residual
    type, extends(band_problem) :: nl_analytic
        type(nl_residual) :: p
        integer :: plant = 0   ! 0 none, 1 = D(2,3,5) += 0.5, 2 = X(3,1) := 0
    contains
        procedure :: fill => nl_analytic_fill
    end type nl_analytic
contains
    subroutine make_problem(p, n, nj)
        type(nl_residual), intent(out) :: p
        integer, intent(in) :: n, nj
        integer :: i
        p%n = n; p%nj = nj
        allocate(p%Am(n,n), p%Bm(n,n), p%Dm(n,n), p%Xm(n,n), p%Ym(n,n), p%b(n))
        call random_number(p%Am); p%Am = 0.3_dp*(2*p%Am - 1)
        call random_number(p%Bm); p%Bm = 0.3_dp*(2*p%Bm - 1)
        call random_number(p%Dm); p%Dm = 0.3_dp*(2*p%Dm - 1)
        call random_number(p%Xm); p%Xm = 0.3_dp*(2*p%Xm - 1)
        call random_number(p%Ym); p%Ym = 0.3_dp*(2*p%Ym - 1)
        call random_number(p%b); p%b = 2*p%b - 1
        do i = 1, n
            p%Bm(i,i) = p%Bm(i,i) + 3 + n
        end do
    end subroutine make_problem

    subroutine nl_residual_eval(self, n, nj, c, F, ierr)
        class(nl_residual), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(out) :: F(n,nj)
        integer, intent(out) :: ierr
        integer :: i, j
        ierr = 0
        if (self%fail) then
            ierr = 3
            return
        end if
        do j = 1, nj
            do i = 1, n
                F(i,j) = 0.5_dp*sin(c(i,j))*c(mod(i, n) + 1, j) - self%b(i)*j + dot_product(self%Bm(i,:), c(:,j))
                if (j > 1) F(i,j) = F(i,j) + dot_product(self%Am(i,:), c(:,j-1))
                if (j < nj) F(i,j) = F(i,j) + dot_product(self%Dm(i,:), c(:,j+1))
                if (j == 1) F(i,j) = F(i,j) + dot_product(self%Xm(i,:), c(:,3)**2)
                if (j == nj) F(i,j) = F(i,j) + dot_product(self%Ym(i,:), c(:,nj-2)**3)
            end do
        end do
    end subroutine nl_residual_eval

    subroutine nl_analytic_fill(self, n, nj, c, A, B, D, G, X, Y, ierr)
        class(nl_analytic), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(inout) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer, intent(out) :: ierr
        real(dp) :: F(n,nj)
        integer :: i, j, k, ip
        call self%p%residual(n, nj, c, F, ierr)
        G = -F
        do j = 1, nj
            B(:,:,j) = self%p%Bm
            if (j > 1) A(:,:,j) = self%p%Am
            if (j < nj) D(:,:,j) = self%p%Dm
            do i = 1, n
                ip = mod(i, n) + 1
                B(i,i,j) = B(i,i,j) + 0.5_dp*cos(c(i,j))*c(ip,j)
                B(i,ip,j) = B(i,ip,j) + 0.5_dp*sin(c(i,j))
            end do
        end do
        do k = 1, n
            X(:,k) = 2*self%p%Xm(:,k)*c(k,3)
            Y(:,k) = 3*self%p%Ym(:,k)*c(k,nj-2)**2
        end do
        if (self%plant == 1) D(2,3,5) = D(2,3,5) + 0.5_dp
        if (self%plant == 2) X(3,1) = 0
    end subroutine nl_analytic_fill
end module fd_test_problem

program test_fd
    use fd_test_problem
    use bandsolver_kernel
    use bandsolver_newton
    use bandsolver_fd
    use band_test_utils, only: check, n_failures, seed_rng
    implicit none
    integer :: n, nj, in, ij, status, evals, i
    integer, parameter :: ns(2) = [1, 3], njs(4) = [3, 4, 5, 10]
    type(nl_analytic) :: an
    type(newton_options) :: opts
    type(newton_result) :: ra, rf
    type(fd_options) :: fdo
    type(jacobian_check) :: chk
    real(dp), allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), c(:,:), cf(:,:)
    real(dp), allocatable :: Ae(:,:,:), Be(:,:,:), De(:,:,:), Ge(:,:)
    real(dp) :: X(3,3), Y(3,3), Xe(3,3), Ye(3,3), rel, scale
    character(100) :: label

    call seed_rng(31)
    do in = 1, size(ns)
        do ij = 1, size(njs)
            n = ns(in); nj = njs(ij)
            call make_problem(an%p, n, nj)
            allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), c(n,nj), Ae(n,n,nj), Be(n,n,nj), De(n,n,nj), Ge(n,nj))
            c = reshape([(0.3_dp + 0.01_dp*i, i = 0, n*nj - 1)], [n, nj])
            Ae = 0; Be = 0; De = 0; Xe = 0; Ye = 0
            call an%fill(n, nj, c, Ae, Be, De, Ge, Xe(1:n,1:n), Ye(1:n,1:n), status)
            call band_fd_jacobian(an%p, n, nj, c, A, B, D, G, X(1:n,1:n), Y(1:n,1:n), status, evaluations=evals)
            A(:,:,1) = 0; D(:,:,nj) = 0
            scale = max(maxval(abs(Be)), maxval(abs(Ae)), maxval(abs(De)))
            rel = max(maxval(abs(A - Ae)), maxval(abs(B - Be)), maxval(abs(D - De)), &
                      maxval(abs(X(1:n,1:n) - Xe(1:n,1:n))), maxval(abs(Y(1:n,1:n) - Ye(1:n,1:n))), &
                      maxval(abs(G - Ge)))/scale
            write(label, '(a,i0,a,i0,a,es8.1,a,i0,a)') 'fd_jacobian n=', n, ' nj=', nj, ' matches analytic (rel ', &
                rel, '), ', evals, ' evals'
            call check(status == BAND_OK .and. rel < 1.0e-6_dp .and. evals == 3*n + 1, trim(label))
            deallocate(A, B, D, G, c, Ae, Be, De, Ge, an%p%Am, an%p%Bm, an%p%Dm, an%p%Xm, an%p%Ym, an%p%b)
        end do
    end do

    ! newton_fd vs analytic Newton.
    n = 3; nj = 40
    call make_problem(an%p, n, nj)
    allocate(c(n,nj), cf(n,nj))
    c = reshape([(0.3_dp + 0.01_dp*i, i = 0, n*nj - 1)], [n, nj]); cf = c
    call band_newton(an, n, nj, c, opts, ra)
    call band_newton_fd(an%p, n, nj, cf, opts, rf)
    print '(a,i0,a,i0,a,i0,a,es9.2)', 'analytic: ', ra%iterations, ' iterations; fd: ', rf%iterations, &
        ' iterations, ', rf%residual_evaluations, ' residual evaluations; max diff ', maxval(abs(c - cf))
    call check(ra%converged .and. rf%converged .and. maxval(abs(c - cf)) < 1.0e-10_dp, 'band_newton_fd matches analytic Newton')
    call check(rf%residual_evaluations == (3*n + 1)*rf%iterations, 'residual evaluations = (3n+1) per iteration')

    ! Jacobian reuse with FD Jacobians: same root, fewer factorizations, and reuse iterations cost
    ! one residual evaluation instead of 3n+1.
    block
        type(newton_options) :: ro
        type(newton_result) :: rr
        real(dp), allocatable :: cr(:,:)
        allocate(cr(n,nj))
        cr = reshape([(0.3_dp + 0.01_dp*i, i = 0, n*nj - 1)], [n, nj])
        ro%jacobian_reuse = .true.
        call band_newton_fd(an%p, n, nj, cr, ro, rr)
        print '(a,i0,a,i0,a,i0,a,es9.2)', 'FD + reuse: ', rr%iterations, ' iterations, ', rr%factorizations, &
            ' factorizations, ', rr%residual_evaluations, ' residual evaluations; max diff ', maxval(abs(cr - c))
        call check(rr%converged .and. maxval(abs(cr - c)) < 1.0e-9_dp, 'band_newton_fd with reuse reaches the same root')
        call check(rr%factorizations < rr%iterations, 'reuse needs fewer factorizations than iterations')
        call check(rr%residual_evaluations == (3*n + 1)*rr%factorizations + (rr%iterations - rr%factorizations), &
            'reuse iterations cost one residual evaluation each')
    end block

    ! check_jacobian: clean, planted D error, missing X entry.
    c = reshape([(0.3_dp + 0.01_dp*i, i = 0, n*nj - 1)], [n, nj])
    an%plant = 0
    call band_check_jacobian(an, n, nj, c, chk, status)
    print '(a,es9.2)', 'correct fill: max error ', max_error(chk)
    call check(status == BAND_OK .and. max_error(chk) < 1.0e-4_dp, 'check_jacobian passes a correct fill')
    an%plant = 1
    call band_check_jacobian(an, n, nj, c, chk, status)
    print '(a,3(i0,1x),es9.2)', 'planted D error at node/row/col: ', chk%D%node, chk%D%row, chk%D%col, chk%D%error
    call check(chk%D%node == 5 .and. chk%D%row == 2 .and. chk%D%col == 3 .and. chk%D%error > 1.0e-2_dp, &
        'planted error located (1-based)')
    call check(max(chk%A%error, chk%B%error, chk%X%error, chk%Y%error) < 1.0e-4_dp, 'other blocks stay clean')
    an%plant = 2
    call band_check_jacobian(an, n, nj, c, chk, status)
    call check(chk%X%row == 3 .and. chk%X%col == 1 .and. chk%X%error > 1.0e-2_dp, 'missing X entry detected')
    an%plant = 0

    ! Errors.
    an%p%fail = .true.
    call band_newton_fd(an%p, n, nj, cf, opts, rf)
    call check(rf%status == BAND_CALLBACK_ERROR, 'residual error propagates through band_newton_fd')
    an%p%fail = .false.
    fdo%rel_step = 0
    allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj))
    call band_fd_jacobian(an%p, n, nj, c, A, B, D, G, X, Y, status, opts=fdo)
    call check(status == BAND_INVALID_ARGUMENT, 'rel_step <= 0 rejected')

    if (n_failures > 0) error stop 1
end program test_fd
