!> Newton iteration on the Appendix C block system.
!>
!> The user extends `band_problem` and implements `fill`, which evaluates at state c the
!> Jacobian blocks A, B, D (and endpoint X, Y) and G = -F(c), the negative residual.
!> Each iteration solves K dc = G and updates c <- c + damping*dc. The iteration has
!> converged when max_ij |dc_ij| / (atol + rtol*|c_ij|) <= 1 (evaluated after the update).
!> Setting max_iter = 1 and require_convergence = .false. reproduces the archival usage
!> of one linearized correction per call.
module bandsolver_newton
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use bandsolver_kernel
    use bandsolver_factor, only: band_factorization, band_factor, band_factor_solve
    implicit none
    private

    public :: band_problem, newton_options, newton_result, band_newton, RESIDUAL_NOT_PROVIDED

    !> Returned by the default band_problem%residual: the driver then falls back to fill.
    integer, parameter :: RESIDUAL_NOT_PROVIDED = -1

    type, abstract :: band_problem
    contains
        procedure(fill_iface), deferred :: fill
        !> Optional residual-only evaluation F(c) (used by Jacobian reuse to skip the blocks).
        !> Override it in your type; the default reports RESIDUAL_NOT_PROVIDED.
        procedure :: residual => band_problem_no_residual
    end type band_problem

    abstract interface
        !> Arrays arrive zeroed; set the entries that apply. A(:,:,1), D(:,:,nj) are unused.
        !> Return ierr /= 0 to abort the iteration with BAND_CALLBACK_ERROR.
        subroutine fill_iface(self, n, nj, c, A, B, D, G, X, Y, ierr)
            import :: band_problem, dp
            class(band_problem), intent(inout) :: self
            integer, intent(in) :: n, nj
            real(dp), intent(in) :: c(n,nj)
            real(dp), intent(inout) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
            integer, intent(out) :: ierr
        end subroutine fill_iface
    end interface

    type :: newton_options
        real(dp) :: rtol = 1.0e-10_dp
        real(dp) :: atol = 1.0e-12_dp
        real(dp) :: damping = 1.0_dp
        integer :: max_iter = 50
        integer :: pivot = PIVOT_PARTIAL
        logical :: require_convergence = .true.
        integer :: kernel = KERNEL_FAST             !< Fortran loop organisation (see bandsolver_kernel)
        !> Jacobian reuse (modified Newton): keep the factorization while updates contract
        !> (step_k <= reuse_contraction*step_{k-1}) and it has been used < reuse_max_iter times.
        logical :: jacobian_reuse = .false.
        integer :: reuse_max_iter = 5
        real(dp) :: reuse_contraction = 0.5_dp
    end type newton_options

    type :: newton_result
        integer :: status = BAND_OK
        integer :: iterations = 0
        integer :: fail_node = 0
        logical :: converged = .false.
        real(dp), allocatable :: update_norm(:)    !< scaled max |dc|/(atol+rtol|c|), per iteration
        real(dp), allocatable :: step_norm(:)      !< max |dc|, per iteration
        real(dp), allocatable :: residual_norm(:)  !< max |G| at the start of each iteration
        integer :: residual_evaluations = 0        !< residual-only calls (reuse) / band_newton_fd total
        integer :: jacobian_evaluations = 0        !< fill calls
        integer :: factorizations = 0              !< block factorizations (full Newton: one per iteration)
    end type newton_result

contains

    subroutine band_problem_no_residual(self, n, nj, c, F, ierr)
        class(band_problem), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(dp), intent(in) :: c(n,nj)
        real(dp), intent(out) :: F(n,nj)
        integer, intent(out) :: ierr
        F = 0 * c
        ierr = RESIDUAL_NOT_PROVIDED + 0 * storage_size(self)   ! (self unused by design)
    end subroutine band_problem_no_residual

    subroutine band_newton(problem, n, nj, c, opts, res)
        class(band_problem), intent(inout) :: problem
        integer, intent(in) :: n, nj
        real(dp), intent(inout) :: c(n,nj)
        type(newton_options), intent(in) :: opts
        type(newton_result), intent(out) :: res

        real(dp), allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), dc(:,:), F(:,:)
        real(dp) :: X(n,n), Y(n,n), snorm, prev_step
        type(band_factorization) :: fac
        logical :: refresh, need_jacobian, have_residual
        integer :: it, ierr, uses

        allocate(res%update_norm(0), res%step_norm(0), res%residual_norm(0))
        if (n < 1 .or. nj < 3 .or. opts%max_iter < 1 .or. .not. (opts%damping > 0 .and. opts%damping <= 1) &
            .or. opts%rtol < 0 .or. opts%atol < 0 .or. (opts%rtol == 0 .and. opts%atol == 0) &
            .or. (opts%jacobian_reuse .and. (opts%reuse_max_iter < 1 .or. .not. (opts%reuse_contraction > 0)))) then
            res%status = BAND_INVALID_ARGUMENT
            return
        end if
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), dc(n,nj), F(n,nj))
        refresh = .true.; have_residual = .true.; uses = 0; prev_step = 0

        do it = 1, opts%max_iter
            res%iterations = it
            need_jacobian = .not. opts%jacobian_reuse .or. refresh
            ierr = RESIDUAL_NOT_PROVIDED
            if (.not. need_jacobian .and. have_residual) then
                call problem%residual(n, nj, c, F, ierr)
                if (ierr == RESIDUAL_NOT_PROVIDED) then
                    have_residual = .false.
                else if (ierr /= 0) then
                    res%status = BAND_CALLBACK_ERROR
                    return
                else
                    res%residual_evaluations = res%residual_evaluations + 1
                    G = -F
                end if
            end if
            if (ierr == RESIDUAL_NOT_PROVIDED) then
                A = 0; B = 0; D = 0; G = 0; X = 0; Y = 0
                call problem%fill(n, nj, c, A, B, D, G, X, Y, ierr)
                if (ierr /= 0) then
                    res%status = BAND_CALLBACK_ERROR
                    return
                end if
                res%jacobian_evaluations = res%jacobian_evaluations + 1
            end if
            res%residual_norm = [res%residual_norm, maxval(abs(G))]

            if (.not. opts%jacobian_reuse) then
                call band_solve(n, nj, A, B, D, G, dc, res%status, X=X, Y=Y, pivot=opts%pivot, kernel=opts%kernel, &
                                fail_node=res%fail_node)
                res%factorizations = res%factorizations + 1
            else
                if (need_jacobian) then
                    call band_factor(n, nj, A, B, D, fac, X=X, Y=Y)
                    res%factorizations = res%factorizations + 1
                    uses = 0
                    refresh = .false.
                    if (fac%status /= BAND_OK) then
                        res%status = fac%status
                        res%fail_node = fac%fail_node
                        return
                    end if
                end if
                call band_factor_solve(fac, G, dc, res%status)
                uses = uses + 1
            end if
            if (res%status /= BAND_OK) return
            c = c + opts%damping*dc
            snorm = maxval(abs(dc))
            res%step_norm = [res%step_norm, snorm]
            res%update_norm = [res%update_norm, maxval(abs(dc)/(opts%atol + opts%rtol*abs(c)))]
            if (res%update_norm(it) <= 1) then
                res%converged = .true.
                exit
            end if
            if (opts%jacobian_reuse) &
                refresh = uses >= opts%reuse_max_iter .or. (uses > 1 .and. snorm > opts%reuse_contraction*prev_step)
            prev_step = snorm
        end do

        if (.not. res%converged .and. opts%require_convergence) res%status = BAND_NOT_CONVERGED
    end subroutine band_newton

end module bandsolver_newton
