!> Finite-difference Jacobians for the Appendix C block system (the "AUTOBAND" idea).
!>
!> The user extends `band_residual_problem` and implements `residual`, which evaluates
!> F(c) (shape (n,nj)). Because F(:,j) depends only on c(:,j-1:j+1) (plus c(:,3) for
!> F(:,1) and c(:,nj-2) for F(:,nj)), perturbing unknown k at every third node at once
!> separates all derivatives: one Jacobian costs 3n + 1 residual evaluations for any nj.
!> Forward differences with step rel_step*max(|c|, typical). The residual MUST respect this
!> stencil; otherwise derivatives are silently mixed. Mirrors cpp/src/fd.cpp operation for
!> operation. Node/row/column indices reported here are 1-based.
module bandsolver_fd
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bandsolver_kernel, only: BAND_OK, BAND_INVALID_ARGUMENT, BAND_CALLBACK_ERROR
    use bandsolver_newton, only: band_problem, newton_options, newton_result, band_newton
    implicit none
    private

    public :: band_residual_problem, fd_options, jacobian_mismatch, jacobian_check
    public :: band_fd_jacobian, band_newton_fd, band_check_jacobian, max_error

    type, abstract :: band_residual_problem
    contains
        procedure(residual_iface), deferred :: residual
    end type band_residual_problem

    abstract interface
        !> Evaluate F(c). Return ierr /= 0 to signal failure (BAND_CALLBACK_ERROR).
        subroutine residual_iface(self, n, nj, c, F, ierr)
            import :: band_residual_problem, dp
            class(band_residual_problem), intent(inout) :: self
            integer, intent(in) :: n, nj
            real(dp), intent(in) :: c(n,nj)
            real(dp), intent(out) :: F(n,nj)
            integer, intent(out) :: ierr
        end subroutine residual_iface
    end interface

    type :: fd_options
        real(dp) :: rel_step = 1.4901161193847656e-08_dp   !< sqrt(machine epsilon)
        real(dp) :: typical = 1.0_dp                        !< step floor scale near zero
    end type fd_options

    !> error = |user - fd| / max(|user|, |fd|, 1e-3*rowscale); rowscale is the largest entry
    !> of the same equation row. X is reported at node 1, Y at node nj.
    type :: jacobian_mismatch
        real(dp) :: error = 0
        integer :: node = 0, row = 0, col = 0
        real(dp) :: user = 0, fd = 0
    end type jacobian_mismatch

    type :: jacobian_check
        type(jacobian_mismatch) :: A, B, D, X, Y
    end type jacobian_check

    !> Adapts a residual problem to the Newton driver's fill interface.
    type, extends(band_problem) :: fd_adapter
        class(band_residual_problem), pointer :: rp => null()
        type(fd_options) :: opts
        integer :: evaluations = 0
    contains
        procedure :: fill => fd_adapter_fill
    end type fd_adapter

    !> Residual F = -G taken from a user fill (for check_jacobian).
    type, extends(band_residual_problem) :: fill_residual
        class(band_problem), pointer :: fp => null()
    contains
        procedure :: residual => fill_residual_eval
    end type fill_residual

contains

    !> Blocks by finite differences and G = -F(c). status: BAND_OK, BAND_INVALID_ARGUMENT,
    !> or BAND_CALLBACK_ERROR. A(:,:,1) and D(:,:,nj) are returned as zero.
    subroutine band_fd_jacobian(problem, n, nj, c, A, B, D, G, X, Y, status, opts, evaluations)
        class(band_residual_problem), intent(inout) :: problem
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(out) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer, intent(out) :: status
        type(fd_options), intent(in), optional :: opts
        integer, intent(out), optional :: evaluations
        type(fd_options) :: o
        real(dp), allocatable :: F0(:,:), Fp(:,:), cp(:,:), h(:)
        real(dp) :: v, step
        integer :: r, k, m, i, ierr, evals

        if (present(opts)) o = opts
        if (present(evaluations)) evaluations = 0
        A = 0; B = 0; D = 0; G = 0; X = 0; Y = 0
        if (n < 1 .or. nj < 3 .or. .not. (o%rel_step > 0) .or. .not. (o%typical > 0)) then
            status = BAND_INVALID_ARGUMENT
            return
        end if
        allocate(F0(n,nj), Fp(n,nj), cp(n,nj), h(nj))
        cp = c
        status = BAND_CALLBACK_ERROR
        call problem%residual(n, nj, c, F0, ierr)
        evals = 1
        if (ierr /= 0) go to 90
        G = -F0
        do r = 1, 3
            do k = 1, n
                do m = r, nj, 3
                    v = c(k,m)
                    step = o%rel_step*max(abs(v), o%typical)*merge(-1.0_dp, 1.0_dp, v < 0)
                    cp(k,m) = v + step
                    h(m) = cp(k,m) - v          ! the step actually taken
                end do
                call problem%residual(n, nj, cp, Fp, ierr)
                evals = evals + 1
                if (ierr /= 0) go to 90
                do m = r, nj, 3
                    do i = 1, n
                        B(i,k,m) = (Fp(i,m) - F0(i,m))/h(m)
                        if (m >= 2) D(i,k,m-1) = (Fp(i,m-1) - F0(i,m-1))/h(m)
                        if (m + 1 <= nj) A(i,k,m+1) = (Fp(i,m+1) - F0(i,m+1))/h(m)
                        if (m == 3) X(i,k) = (Fp(i,1) - F0(i,1))/h(m)
                        if (m == nj - 2) Y(i,k) = (Fp(i,nj) - F0(i,nj))/h(m)
                    end do
                    cp(k,m) = c(k,m)
                end do
            end do
        end do
        status = BAND_OK
90      if (present(evaluations)) evaluations = evals
    end subroutine band_fd_jacobian

    subroutine fd_adapter_fill(self, n, nj, c, A, B, D, G, X, Y, ierr)
        class(fd_adapter), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(inout) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer, intent(out) :: ierr
        integer :: status, evals
        call band_fd_jacobian(self%rp, n, nj, c, A, B, D, G, X, Y, status, self%opts, evals)
        self%evaluations = self%evaluations + evals
        ierr = merge(0, 1, status == BAND_OK)
    end subroutine fd_adapter_fill

    !> Newton iteration with finite-difference Jacobians; sets res%residual_evaluations.
    subroutine band_newton_fd(problem, n, nj, c, opts, res, fd_opts)
        class(band_residual_problem), intent(inout), target :: problem
        integer, intent(in) :: n, nj
        real(dp), intent(inout) :: c(n,nj)
        type(newton_options), intent(in) :: opts
        type(newton_result), intent(out) :: res
        type(fd_options), intent(in), optional :: fd_opts
        type(fd_adapter) :: adapter
        adapter%rp => problem
        if (present(fd_opts)) adapter%opts = fd_opts
        if (.not. (adapter%opts%rel_step > 0) .or. .not. (adapter%opts%typical > 0)) then
            allocate(res%update_norm(0), res%step_norm(0), res%residual_norm(0))
            res%status = BAND_INVALID_ARGUMENT
            return
        end if
        call band_newton(adapter, n, nj, c, opts, res)
        res%residual_evaluations = adapter%evaluations
    end subroutine band_newton_fd

    subroutine fill_residual_eval(self, n, nj, c, F, ierr)
        class(fill_residual), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(out) :: F(n,nj)
        integer, intent(out) :: ierr
        real(dp), allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:)
        real(dp) :: X(n,n), Y(n,n)
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj))
        A = 0; B = 0; D = 0; G = 0; X = 0; Y = 0
        call self%fp%fill(n, nj, c, A, B, D, G, X, Y, ierr)
        F = -G
    end subroutine fill_residual_eval

    !> Compare a hand-written fill against finite differences of its own G (F = -G) at c.
    !> A correct Jacobian typically scores 1e-8 to 1e-5; above ~1e-3 indicates a bug.
    subroutine band_check_jacobian(problem, n, nj, c, check, status, opts)
        class(band_problem), intent(inout), target :: problem
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        type(jacobian_check), intent(out) :: check
        integer, intent(out) :: status
        type(fd_options), intent(in), optional :: opts
        type(fill_residual) :: fr
        real(dp), allocatable :: Au(:,:,:), Bu(:,:,:), Du(:,:,:), Gu(:,:), Af(:,:,:), Bf(:,:,:), Df(:,:,:), Gf(:,:)
        real(dp), allocatable :: rowscale(:,:)
        real(dp) :: Xu(n,n), Yu(n,n), Xf(n,n), Yf(n,n)
        integer :: i, j, k, ierr

        if (n < 1 .or. nj < 3) then
            status = BAND_INVALID_ARGUMENT
            return
        end if
        allocate(Au(n,n,nj), Bu(n,n,nj), Du(n,n,nj), Gu(n,nj), Af(n,n,nj), Bf(n,n,nj), Df(n,n,nj), Gf(n,nj))
        allocate(rowscale(n,nj))
        Au = 0; Bu = 0; Du = 0; Gu = 0; Xu = 0; Yu = 0
        call problem%fill(n, nj, c, Au, Bu, Du, Gu, Xu, Yu, ierr)
        if (ierr /= 0) then
            status = BAND_CALLBACK_ERROR
            return
        end if
        fr%fp => problem
        call band_fd_jacobian(fr, n, nj, c, Af, Bf, Df, Gf, Xf, Yf, status, opts)
        if (status /= BAND_OK) return

        rowscale = 0
        do j = 1, nj
            do i = 1, n
                do k = 1, n
                    if (j > 1) rowscale(i,j) = max(rowscale(i,j), abs(Au(i,k,j)), abs(Af(i,k,j)))
                    rowscale(i,j) = max(rowscale(i,j), abs(Bu(i,k,j)), abs(Bf(i,k,j)))
                    if (j < nj) rowscale(i,j) = max(rowscale(i,j), abs(Du(i,k,j)), abs(Df(i,k,j)))
                    if (j == 1) rowscale(i,j) = max(rowscale(i,j), abs(Xu(i,k)), abs(Xf(i,k)))
                    if (j == nj) rowscale(i,j) = max(rowscale(i,j), abs(Yu(i,k)), abs(Yf(i,k)))
                end do
            end do
        end do
        do j = 1, nj
            do i = 1, n
                do k = 1, n
                    if (j > 1) call consider(check%A, j, i, k, Au(i,k,j), Af(i,k,j))
                    call consider(check%B, j, i, k, Bu(i,k,j), Bf(i,k,j))
                    if (j < nj) call consider(check%D, j, i, k, Du(i,k,j), Df(i,k,j))
                end do
            end do
        end do
        do i = 1, n
            do k = 1, n
                call consider(check%X, 1, i, k, Xu(i,k), Xf(i,k))
                call consider(check%Y, nj, i, k, Yu(i,k), Yf(i,k))
            end do
        end do
    contains
        subroutine consider(mm, node, row, col, u, f)
            type(jacobian_mismatch), intent(inout) :: mm
            integer, intent(in) :: node, row, col
            real(dp), intent(in) :: u, f
            real(dp) :: e
            e = abs(u - f)/max(abs(u), abs(f), max(1.0e-3_dp*rowscale(row,node), tiny(1.0_dp)))
            if (mm%node == 0 .or. e > mm%error) mm = jacobian_mismatch(e, node, row, col, u, f)
        end subroutine consider
    end subroutine band_check_jacobian

    pure real(dp) function max_error(check)
        type(jacobian_check), intent(in) :: check
        max_error = max(check%A%error, check%B%error, check%D%error, check%X%error, check%Y%error)
    end function max_error

end module bandsolver_fd
