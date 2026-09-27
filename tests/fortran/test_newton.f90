! T3: Newton driver — steady 1-D diffusion -(D(c) c')' = f on [0,1], c(0)=c(1)=0,
! manufactured solution c* = sin(pi x); D = 1 (linear) or D = 1 + c^2 (nonlinear).
module newton_test_problem
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bandsolver_newton
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type, extends(band_problem) :: diffusion_bvp
        logical :: nonlinear = .true.
        logical :: fail = .false.
    contains
        procedure :: fill => diffusion_fill
    end type diffusion_bvp
contains
    pure real(dp) function dcoef(self, c)
        class(diffusion_bvp), intent(in) :: self
        real(dp), intent(in) :: c
        dcoef = 1
        if (self%nonlinear) dcoef = 1 + c*c
    end function dcoef
    pure real(dp) function dprime(self, c)
        class(diffusion_bvp), intent(in) :: self
        real(dp), intent(in) :: c
        dprime = 0
        if (self%nonlinear) dprime = 2*c
    end function dprime
    pure real(dp) function source(self, x)
        class(diffusion_bvp), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: c, cx, cxx
        c = sin(pi*x); cx = pi*cos(pi*x); cxx = -pi*pi*sin(pi*x)
        source = -(dprime(self, c)*cx*cx + dcoef(self, c)*cxx)
    end function source

    subroutine diffusion_fill(self, n, nj, c, A, B, D, G, X, Y, ierr)
        class(diffusion_bvp), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(inout) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer, intent(out) :: ierr
        real(dp) :: h, cm, cp, Dminus, Dplus, F
        integer :: j
        ierr = 0
        if (self%fail) then
            ierr = 42
            return
        end if
        h = 1.0_dp/(nj - 1)
        B(1,1,1) = 1;  G(1,1) = -c(1,1)
        B(1,1,nj) = 1; G(1,nj) = -c(1,nj)
        do j = 2, nj - 1
            cm = (c(1,j) + c(1,j-1))/2; cp = (c(1,j) + c(1,j+1))/2
            Dminus = dcoef(self, cm); Dplus = dcoef(self, cp)
            F = -(Dplus*(c(1,j+1) - c(1,j)) - Dminus*(c(1,j) - c(1,j-1)))/h**2 - source(self, (j-1)*h)
            G(1,j) = -F
            A(1,1,j) = (0.5_dp*dprime(self, cm)*(c(1,j) - c(1,j-1)) - Dminus)/h**2
            D(1,1,j) = -(Dplus + 0.5_dp*dprime(self, cp)*(c(1,j+1) - c(1,j)))/h**2
            B(1,1,j) = (Dplus + Dminus - 0.5_dp*dprime(self, cp)*(c(1,j+1) - c(1,j)) &
                        + 0.5_dp*dprime(self, cm)*(c(1,j) - c(1,j-1)))/h**2
        end do
        if (.false.) X = Y   ! endpoint blocks unused in this problem
    end subroutine diffusion_fill
end module newton_test_problem

program test_newton
    use newton_test_problem
    use bandsolver_kernel
    use band_test_utils, only: check, n_failures
    implicit none
    integer, parameter :: nj = 101
    type(diffusion_bvp) :: prob
    type(newton_options) :: opts
    type(newton_result) :: res
    real(dp) :: c(1,nj), cexact(1,nj), c_partial(1,nj)
    integer :: j, k
    logical :: quad

    cexact(1,:) = [(sin(pi*(j-1)/(nj-1.0_dp)), j = 1, nj)]

    ! Linear problem: exact after one solve, detected at iteration 2.
    prob%nonlinear = .false.
    c = 0
    call band_newton(prob, 1, nj, c, opts, res)
    print '(a,i0,a,2es10.2)', 'linear: iterations=', res%iterations, ' residual norms=', res%residual_norm
    call check(res%status == BAND_OK .and. res%converged .and. res%iterations == 2, 'linear converges at iteration 2')
    call check(res%residual_norm(2) < 1.0e-10_dp*res%residual_norm(1), 'linear residual eliminated by first solve')

    ! Nonlinear problem: quadratic convergence and O(h^2) discretization error.
    prob%nonlinear = .true.
    c = 0
    call band_newton(prob, 1, nj, c, opts, res)
    print '(a,i0)', 'nonlinear: iterations=', res%iterations
    print '(a,*(es10.2))', '  step norms: ', res%step_norm
    call check(res%status == BAND_OK .and. res%converged .and. res%iterations <= 8, 'nonlinear converges in <= 8 iterations')
    quad = .true.
    do k = 1, size(res%step_norm) - 1
        if (res%step_norm(k) < 0.1_dp .and. res%step_norm(k+1) > 1.0e-13_dp) &
            quad = quad .and. res%step_norm(k+1) <= 10*res%step_norm(k)**2
    end do
    call check(quad, 'quadratic convergence: s(k+1) <= 10 s(k)^2 once s(k) < 0.1')
    print '(a,es10.2)', '  max error vs sin(pi x): ', maxval(abs(c - cexact))
    call check(maxval(abs(c - cexact)) < 1.0e-3_dp, 'nonlinear solution within O(h^2) of exact')
    c_partial = c

    opts%pivot = PIVOT_LEGACY
    c = 0
    call band_newton(prob, 1, nj, c, opts, res)
    call check(res%converged .and. maxval(abs(c - c_partial)) < 1.0e-12_dp, 'legacy pivot Newton agrees')
    opts%pivot = PIVOT_PARTIAL

    ! Archival usage: one correction, no convergence requirement.
    opts%max_iter = 1; opts%require_convergence = .false.
    c = 0
    call band_newton(prob, 1, nj, c, opts, res)
    call check(res%status == BAND_OK .and. res%iterations == 1 .and. .not. res%converged, 'one-step legacy usage returns OK')
    opts%require_convergence = .true.
    c = 0
    call band_newton(prob, 1, nj, c, opts, res)
    call check(res%status == BAND_NOT_CONVERGED, 'max_iter exhausted -> NOT_CONVERGED')
    opts = newton_options()

    prob%fail = .true.
    call band_newton(prob, 1, nj, c, opts, res)
    call check(res%status == BAND_CALLBACK_ERROR, 'callback error propagates')
    prob%fail = .false.

    opts%damping = 0
    call band_newton(prob, 1, nj, c, opts, res)
    call check(res%status == BAND_INVALID_ARGUMENT, 'damping = 0 rejected')

    if (n_failures > 0) error stop 1
end program test_newton
