!> bind(C) interface declared in include/bandsolver_f.h. C row-major blocks [nj][n][n]
!> are transposed per block on entry/exit; [nj][n] vectors share the Fortran (n,nj) layout.
module bandsolver_capi
    use, intrinsic :: iso_c_binding
    use bandsolver_kernel
    use bandsolver_newton
    implicit none
    private

    type, bind(c) :: c_newton_options
        real(c_double) :: rtol, atol, damping
        integer(c_int) :: max_iter, pivot, require_convergence
    end type c_newton_options

    type, bind(c) :: c_newton_result
        integer(c_int) :: status, iterations, converged, fail_node
        real(c_double) :: update_norm, step_norm, residual_norm
    end type c_newton_result

    abstract interface
        integer(c_int) function c_fill_iface(n, nj, c, A, B, D, G, X, Y, ctx) bind(c)
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: n, nj
            real(c_double), intent(in) :: c(*)
            real(c_double) :: A(*), B(*), D(*), G(*), X(*), Y(*)
            type(c_ptr), value :: ctx
        end function c_fill_iface
    end interface

    !> Adapts a C callback + context to the Fortran band_problem interface.
    type, extends(band_problem) :: c_problem
        type(c_funptr) :: cfill = c_null_funptr
        type(c_ptr) :: ctx = c_null_ptr
        real(c_double), allocatable :: At(:,:,:), Bt(:,:,:), Dt(:,:,:), Xt(:,:), Yt(:,:)
    contains
        procedure :: fill => c_problem_fill
    end type c_problem

contains

    pure subroutine to_fortran_blocks(n, nj, src, dst)
        integer, intent(in) :: n, nj
        real(c_double), intent(in) :: src(n,n,nj)   ! C row-major view: src(k,i,j) = block_j(i,k)
        real(c_double), intent(out) :: dst(n,n,nj)
        integer :: j
        do j = 1, nj
            dst(:,:,j) = transpose(src(:,:,j))
        end do
    end subroutine to_fortran_blocks

    subroutine bandsolver_f_default_options(opts) bind(c, name='bandsolver_f_default_options')
        type(c_newton_options), intent(out) :: opts
        type(newton_options) :: d
        opts = c_newton_options(d%rtol, d%atol, d%damping, d%max_iter, d%pivot, 1)
    end subroutine bandsolver_f_default_options

    integer(c_int) function bandsolver_f_solve(n, nj, A, B, D, G, X, Y, pivot, dc, fail_node, &
            min_rel_pivot) bind(c, name='bandsolver_f_solve') result(status)
        integer(c_int), value :: n, nj, pivot
        real(c_double), intent(in) :: A(*), B(*), D(*), G(*)
        type(c_ptr), value :: X, Y
        real(c_double), intent(out) :: dc(*)
        integer(c_int), intent(out) :: fail_node
        real(c_double), intent(out) :: min_rel_pivot
        real(c_double), allocatable :: Af(:,:,:), Bf(:,:,:), Df(:,:,:)
        real(c_double), pointer :: Xp(:,:), Yp(:,:)
        real(c_double), allocatable :: dcf(:,:)
        real(c_double) :: Xf(max(n,1),max(n,1)), Yf(max(n,1),max(n,1))
        integer :: fnode

        fail_node = 0
        min_rel_pivot = 0
        if (n < 1 .or. nj < 3) then
            status = BAND_INVALID_ARGUMENT
            return
        end if
        allocate(Af(n,n,nj), Bf(n,n,nj), Df(n,n,nj))
        call to_fortran_blocks(n, nj, A(1:n*n*nj), Af)
        call to_fortran_blocks(n, nj, B(1:n*n*nj), Bf)
        call to_fortran_blocks(n, nj, D(1:n*n*nj), Df)
        Xf = 0; Yf = 0
        if (c_associated(X)) then
            call c_f_pointer(X, Xp, [n, n]); Xf = transpose(Xp)
        end if
        if (c_associated(Y)) then
            call c_f_pointer(Y, Yp, [n, n]); Yf = transpose(Yp)
        end if
        allocate(dcf(n,nj))
        call band_solve(n, nj, Af, Bf, Df, reshape(G(1:n*nj), [n, nj]), dcf, status, &
                        X=Xf, Y=Yf, pivot=pivot, fail_node=fnode, min_rel_pivot=min_rel_pivot)
        dc(1:n*nj) = reshape(dcf, [n*nj])
        fail_node = fnode
    end function bandsolver_f_solve

    subroutine c_problem_fill(self, n, nj, c, A, B, D, G, X, Y, ierr)
        class(c_problem), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(c_double), intent(in) :: c(n,nj)
        real(c_double), intent(inout) :: A(n,n,nj), B(n,n,nj), D(n,n,nj), G(n,nj), X(n,n), Y(n,n)
        integer, intent(out) :: ierr
        procedure(c_fill_iface), pointer :: cfill
        call c_f_procpointer(self%cfill, cfill)
        self%At = 0; self%Bt = 0; self%Dt = 0; self%Xt = 0; self%Yt = 0
        ierr = cfill(int(n, c_int), int(nj, c_int), c, self%At, self%Bt, self%Dt, G, &
                          self%Xt, self%Yt, self%ctx)
        call to_fortran_blocks(n, nj, self%At, A)
        call to_fortran_blocks(n, nj, self%Bt, B)
        call to_fortran_blocks(n, nj, self%Dt, D)
        X = transpose(self%Xt)
        Y = transpose(self%Yt)
    end subroutine c_problem_fill

    integer(c_int) function bandsolver_f_newton(n, nj, fill, ctx, c, opts, res, update_history, &
            step_history, residual_history) bind(c, name='bandsolver_f_newton') result(status)
        integer(c_int), value :: n, nj
        type(c_funptr), value :: fill
        type(c_ptr), value :: ctx
        real(c_double), intent(inout), target :: c(*)
        type(c_newton_options), intent(in) :: opts
        type(c_newton_result), intent(out) :: res
        type(c_ptr), value :: update_history, step_history, residual_history
        type(c_problem) :: prob
        type(newton_options) :: fo
        type(newton_result) :: fr
        real(c_double), pointer :: cf(:,:)

        res = c_newton_result(BAND_INVALID_ARGUMENT, 0, 0, 0, 0, 0, 0)
        status = BAND_INVALID_ARGUMENT
        if (n < 1 .or. nj < 3 .or. .not. c_associated(fill)) return
        prob%cfill = fill
        prob%ctx = ctx
        allocate(prob%At(n,n,nj), prob%Bt(n,n,nj), prob%Dt(n,n,nj), prob%Xt(n,n), prob%Yt(n,n))
        fo = newton_options(opts%rtol, opts%atol, opts%damping, opts%max_iter, opts%pivot, &
                            opts%require_convergence /= 0)
        cf(1:n,1:nj) => c(1:n*nj)
        call band_newton(prob, n, nj, cf, fo, fr)

        status = fr%status
        res%status = fr%status
        res%iterations = fr%iterations
        res%converged = merge(1, 0, fr%converged)
        res%fail_node = fr%fail_node
        if (size(fr%update_norm) > 0) res%update_norm = fr%update_norm(size(fr%update_norm))
        if (size(fr%step_norm) > 0) res%step_norm = fr%step_norm(size(fr%step_norm))
        if (size(fr%residual_norm) > 0) res%residual_norm = fr%residual_norm(size(fr%residual_norm))
        call copy_history(fr%update_norm, update_history)
        call copy_history(fr%step_norm, step_history)
        call copy_history(fr%residual_norm, residual_history)
    end function bandsolver_f_newton

    subroutine copy_history(h, dst)
        real(c_double), intent(in) :: h(:)
        type(c_ptr), value :: dst
        real(c_double), pointer :: p(:)
        if (.not. c_associated(dst) .or. size(h) == 0) return
        call c_f_pointer(dst, p, [size(h)])
        p = h
    end subroutine copy_history

end module bandsolver_capi
