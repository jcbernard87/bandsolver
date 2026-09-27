! Quickstart: solve -c'' = 1 on [0,1], c(0) = c(1) = 0 (exact: x(1-x)/2, which the
! three-point stencil reproduces exactly), first as a linear solve, then via Newton.
module quickstart_problem
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bandsolver_newton, only: band_problem
    implicit none
    type, extends(band_problem) :: poisson
    contains
        procedure :: fill
    end type poisson
contains
    subroutine fill(self, n, nj, c, A, B, D, G, X, Y, ierr)
        class(poisson), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(inout) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer, intent(out) :: ierr
        real(dp) :: h
        integer :: j
        h = 1.0_dp/(nj - 1)
        B(1,1,1) = 1;  G(1,1) = -c(1,1)            ! G = -F (negative residual)
        B(1,1,nj) = 1; G(1,nj) = -c(1,nj)
        do j = 2, nj - 1
            A(1,1,j) = -1/h**2; B(1,1,j) = 2/h**2; D(1,1,j) = -1/h**2
            G(1,j) = (c(1,j+1) - 2*c(1,j) + c(1,j-1))/h**2 + 1
        end do
        ierr = 0
    end subroutine fill
end module quickstart_problem

program fortran_quickstart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bandsolver_kernel
    use bandsolver_newton
    use quickstart_problem
    implicit none
    integer, parameter :: n = 1, nj = 11
    real(dp) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), dc(n,nj), c(n,nj), x(nj)
    type(poisson) :: prob
    type(newton_options) :: opts
    type(newton_result) :: res
    integer :: status, j

    x = [((j - 1)/real(nj - 1, dp), j = 1, nj)]

    ! 1) Linear solve: assemble A, B, D, G directly.
    A = 0; B = 0; D = 0; G = 0
    B(1,1,1) = 1; B(1,1,nj) = 1
    do j = 2, nj - 1
        A(1,1,j) = -(nj - 1)**2; B(1,1,j) = 2*(nj - 1)**2; D(1,1,j) = -(nj - 1)**2; G(1,j) = 1
    end do
    call band_solve(n, nj, A, B, D, G, dc, status)
    print '(a,i0,a,es9.2)', 'band_solve status=', status, '  max error=', maxval(abs(dc(1,:) - x*(1 - x)/2))

    ! 2) Newton with a fill callback.
    c = 0
    call band_newton(prob, n, nj, c, opts, res)
    print '(a,i0,a,l1,a,i0,a,es9.2)', 'band_newton status=', res%status, '  converged=', res%converged, &
        '  iterations=', res%iterations, '  max error=', maxval(abs(c(1,:) - x*(1 - x)/2))
end program fortran_quickstart
