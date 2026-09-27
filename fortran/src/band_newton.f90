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
    implicit none
    private

    public :: band_problem, newton_options, newton_result, band_newton

    type, abstract :: band_problem
    contains
        procedure(fill_iface), deferred :: fill
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
    end type newton_options

    type :: newton_result
        integer :: status = BAND_OK
        integer :: iterations = 0
        integer :: fail_node = 0
        logical :: converged = .false.
        real(dp), allocatable :: update_norm(:)    !< scaled max |dc|/(atol+rtol|c|), per iteration
        real(dp), allocatable :: step_norm(:)      !< max |dc|, per iteration
        real(dp), allocatable :: residual_norm(:)  !< max |G| at the start of each iteration
        integer :: residual_evaluations = 0        !< set by band_newton_fd only
    end type newton_result

contains

    subroutine band_newton(problem, n, nj, c, opts, res)
        class(band_problem), intent(inout) :: problem
        integer, intent(in) :: n, nj
        real(dp), intent(inout) :: c(n,nj)
        type(newton_options), intent(in) :: opts
        type(newton_result), intent(out) :: res

        real(dp), allocatable :: A(:,:,:), B(:,:,:), D(:,:,:), G(:,:), dc(:,:)
        real(dp) :: X(n,n), Y(n,n)
        integer :: it, ierr

        allocate(res%update_norm(0), res%step_norm(0), res%residual_norm(0))
        if (n < 1 .or. nj < 3 .or. opts%max_iter < 1 .or. .not. (opts%damping > 0 .and. opts%damping <= 1) &
            .or. opts%rtol < 0 .or. opts%atol < 0 .or. (opts%rtol == 0 .and. opts%atol == 0)) then
            res%status = BAND_INVALID_ARGUMENT
            return
        end if
        allocate(A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), dc(n,nj))

        do it = 1, opts%max_iter
            A = 0; B = 0; D = 0; G = 0; X = 0; Y = 0
            call problem%fill(n, nj, c, A, B, D, G, X, Y, ierr)
            res%iterations = it
            if (ierr /= 0) then
                res%status = BAND_CALLBACK_ERROR
                return
            end if
            res%residual_norm = [res%residual_norm, maxval(abs(G))]
            call band_solve(n, nj, A, B, D, G, dc, res%status, X=X, Y=Y, pivot=opts%pivot, &
                            fail_node=res%fail_node)
            if (res%status /= BAND_OK) return
            c = c + opts%damping*dc
            res%step_norm = [res%step_norm, maxval(abs(dc))]
            res%update_norm = [res%update_norm, maxval(abs(dc)/(opts%atol + opts%rtol*abs(c)))]
            if (res%update_norm(it) <= 1) then
                res%converged = .true.
                exit
            end if
        end do

        if (.not. res%converged .and. opts%require_convergence) res%status = BAND_NOT_CONVERGED
    end subroutine band_newton

end module bandsolver_newton
