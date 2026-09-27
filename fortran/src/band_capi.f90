!> bind(C) interface declared in include/bandsolver_f.h. C row-major blocks [nj][n][n]
!> are transposed per block on entry/exit; [nj][n] vectors share the Fortran (n,nj) layout.
module bandsolver_capi
    use, intrinsic :: iso_c_binding
    use bandsolver_kernel
    use bandsolver_newton
    use bandsolver_fd
    use bandsolver_factor
    implicit none
    private

    type, bind(c) :: c_newton_options
        real(c_double) :: rtol, atol, damping
        integer(c_int) :: max_iter, pivot, require_convergence, kernel
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

    type, bind(c) :: c_fd_options
        real(c_double) :: rel_step, typical
    end type c_fd_options

    type, bind(c) :: c_mismatch
        real(c_double) :: error
        integer(c_int) :: node, row, col
        real(c_double) :: user, fd
    end type c_mismatch

    type, bind(c) :: c_jacobian_check
        type(c_mismatch) :: A, B, D, X, Y
    end type c_jacobian_check

    abstract interface
        integer(c_int) function c_residual_iface(n, nj, c, F, ctx) bind(c)
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: n, nj
            real(c_double), intent(in) :: c(*)
            real(c_double) :: F(*)
            type(c_ptr), value :: ctx
        end function c_residual_iface
    end interface

    !> Adapts a C residual callback + context to band_residual_problem.
    type, extends(band_residual_problem) :: c_residual_problem
        type(c_funptr) :: cres = c_null_funptr
        type(c_ptr) :: ctx = c_null_ptr
    contains
        procedure :: residual => c_residual_eval
    end type c_residual_problem

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
        opts = c_newton_options(d%rtol, d%atol, d%damping, d%max_iter, d%pivot, 1, d%kernel)
    end subroutine bandsolver_f_default_options

    integer(c_int) function bandsolver_f_solve(n, nj, A, B, D, G, X, Y, pivot, dc, fail_node, &
            min_rel_pivot) bind(c, name='bandsolver_f_solve') result(status)
        integer(c_int), value :: n, nj, pivot
        real(c_double), intent(in) :: A(*), B(*), D(*), G(*)
        type(c_ptr), value :: X, Y
        real(c_double), intent(out) :: dc(*)
        integer(c_int), intent(out) :: fail_node
        real(c_double), intent(out) :: min_rel_pivot
        status = bandsolver_f_solve_kernel(n, nj, A, B, D, G, X, Y, pivot, int(KERNEL_FAST, c_int), dc, &
                                           fail_node, min_rel_pivot)
    end function bandsolver_f_solve

    integer(c_int) function bandsolver_f_solve_kernel(n, nj, A, B, D, G, X, Y, pivot, kernel, dc, fail_node, &
            min_rel_pivot) bind(c, name='bandsolver_f_solve_kernel') result(status)
        integer(c_int), value :: n, nj, pivot, kernel
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
                        X=Xf, Y=Yf, pivot=pivot, fail_node=fnode, min_rel_pivot=min_rel_pivot, kernel=kernel)
        dc(1:n*nj) = reshape(dcf, [n*nj])
        fail_node = fnode
    end function bandsolver_f_solve_kernel

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
                            opts%require_convergence /= 0, opts%kernel)
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

    subroutine c_residual_eval(self, n, nj, c, F, ierr)
        class(c_residual_problem), intent(inout) :: self
        integer, intent(in) :: n, nj
        real(c_double), intent(in) :: c(n,nj)
        real(c_double), intent(out) :: F(n,nj)
        integer, intent(out) :: ierr
        procedure(c_residual_iface), pointer :: cres
        call c_f_procpointer(self%cres, cres)
        F = 0
        ierr = cres(int(n, c_int), int(nj, c_int), c, F, self%ctx)
    end subroutine c_residual_eval

    function get_fd_options(p) result(o)
        type(c_ptr), value :: p
        type(fd_options) :: o
        type(c_fd_options), pointer :: co
        if (c_associated(p)) then
            call c_f_pointer(p, co)
            o = fd_options(co%rel_step, co%typical)
        end if
    end function get_fd_options

    subroutine bandsolver_f_default_fd_options(opts) bind(c, name='bandsolver_f_default_fd_options')
        type(c_fd_options), intent(out) :: opts
        type(fd_options) :: d
        opts = c_fd_options(d%rel_step, d%typical)
    end subroutine bandsolver_f_default_fd_options

    integer(c_int) function bandsolver_f_fd_jacobian(n, nj, residual, ctx, c, fd_opts, A, B, D, G, X, Y, &
            evaluations) bind(c, name='bandsolver_f_fd_jacobian') result(status)
        integer(c_int), value :: n, nj
        type(c_funptr), value :: residual
        type(c_ptr), value :: ctx, fd_opts, evaluations
        real(c_double), intent(in) :: c(*)
        real(c_double), intent(out) :: A(*), B(*), D(*), G(*), X(*), Y(*)
        type(c_residual_problem) :: prob
        real(c_double), allocatable :: Af(:,:,:), Bf(:,:,:), Df(:,:,:), Gf(:,:), Xf(:,:), Yf(:,:)
        integer(c_long), pointer :: ev
        integer :: fstatus, evals

        status = BAND_INVALID_ARGUMENT
        if (n < 1 .or. nj < 3 .or. .not. c_associated(residual)) return
        prob%cres = residual
        prob%ctx = ctx
        allocate(Af(n,n,nj), Bf(n,n,nj), Df(n,n,nj), Gf(n,nj), Xf(n,n), Yf(n,n))
        call band_fd_jacobian(prob, n, nj, reshape(c(1:n*nj), [n, nj]), Af, Bf, Df, Gf, Xf, Yf, fstatus, &
                              get_fd_options(fd_opts), evals)
        call to_fortran_blocks(n, nj, Af, A(1:n*n*nj))
        call to_fortran_blocks(n, nj, Bf, B(1:n*n*nj))
        call to_fortran_blocks(n, nj, Df, D(1:n*n*nj))
        G(1:n*nj) = reshape(Gf, [n*nj])
        X(1:n*n) = reshape(transpose(Xf), [n*n])
        Y(1:n*n) = reshape(transpose(Yf), [n*n])
        if (c_associated(evaluations)) then
            call c_f_pointer(evaluations, ev)
            ev = evals
        end if
        status = fstatus
    end function bandsolver_f_fd_jacobian

    integer(c_int) function bandsolver_f_newton_fd(n, nj, residual, ctx, c, opts, fd_opts, res, update_history, &
            step_history, residual_history, evaluations) bind(c, name='bandsolver_f_newton_fd') result(status)
        integer(c_int), value :: n, nj
        type(c_funptr), value :: residual
        type(c_ptr), value :: ctx, fd_opts, update_history, step_history, residual_history, evaluations
        real(c_double), intent(inout), target :: c(*)
        type(c_newton_options), intent(in) :: opts
        type(c_newton_result), intent(out) :: res
        type(c_residual_problem), target :: prob
        type(newton_options) :: fo
        type(newton_result) :: fr
        real(c_double), pointer :: cf(:,:)
        integer(c_long), pointer :: ev

        res = c_newton_result(BAND_INVALID_ARGUMENT, 0, 0, 0, 0, 0, 0)
        status = BAND_INVALID_ARGUMENT
        if (n < 1 .or. nj < 3 .or. .not. c_associated(residual)) return
        prob%cres = residual
        prob%ctx = ctx
        fo = newton_options(opts%rtol, opts%atol, opts%damping, opts%max_iter, opts%pivot, &
                            opts%require_convergence /= 0, opts%kernel)
        cf(1:n,1:nj) => c(1:n*nj)
        call band_newton_fd(prob, n, nj, cf, fo, fr, get_fd_options(fd_opts))

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
        if (c_associated(evaluations)) then
            call c_f_pointer(evaluations, ev)
            ev = fr%residual_evaluations
        end if
    end function bandsolver_f_newton_fd

    integer(c_int) function bandsolver_f_check_jacobian(n, nj, fill, ctx, c, fd_opts, check) &
            bind(c, name='bandsolver_f_check_jacobian') result(status)
        integer(c_int), value :: n, nj
        type(c_funptr), value :: fill
        type(c_ptr), value :: ctx, fd_opts
        real(c_double), intent(in) :: c(*)
        type(c_jacobian_check), intent(out) :: check
        type(c_problem), target :: prob
        type(jacobian_check) :: fc
        integer :: fstatus

        check = c_jacobian_check(c_mismatch(0, 0, 0, 0, 0, 0), c_mismatch(0, 0, 0, 0, 0, 0), &
            c_mismatch(0, 0, 0, 0, 0, 0), c_mismatch(0, 0, 0, 0, 0, 0), c_mismatch(0, 0, 0, 0, 0, 0))
        status = BAND_INVALID_ARGUMENT
        if (n < 1 .or. nj < 3 .or. .not. c_associated(fill)) return
        prob%cfill = fill
        prob%ctx = ctx
        allocate(prob%At(n,n,nj), prob%Bt(n,n,nj), prob%Dt(n,n,nj), prob%Xt(n,n), prob%Yt(n,n))
        call band_check_jacobian(prob, n, nj, reshape(c(1:n*nj), [n, nj]), fc, fstatus, get_fd_options(fd_opts))
        check = c_jacobian_check(cm(fc%A), cm(fc%B), cm(fc%D), cm(fc%X), cm(fc%Y))
        status = fstatus
    contains
        type(c_mismatch) function cm(m)
            type(jacobian_mismatch), intent(in) :: m
            cm = c_mismatch(m%error, m%node, m%row, m%col, m%user, m%fd)
        end function cm
    end function bandsolver_f_check_jacobian

    !> Factor the block matrix; returns an opaque handle (free with bandsolver_f_factor_free).
    integer(c_int) function bandsolver_f_factor(n, nj, A, B, D, X, Y, handle, fail_node) &
            bind(c, name='bandsolver_f_factor') result(status)
        integer(c_int), value :: n, nj
        real(c_double), intent(in) :: A(*), B(*), D(*)
        type(c_ptr), value :: X, Y
        type(c_ptr), intent(out) :: handle
        integer(c_int), intent(out) :: fail_node
        type(band_factorization), pointer :: f
        real(c_double), allocatable :: Af(:,:,:), Bf(:,:,:), Df(:,:,:)
        real(c_double), pointer :: Xp(:,:), Yp(:,:)
        real(c_double) :: Xf(max(n,1),max(n,1)), Yf(max(n,1),max(n,1))

        handle = c_null_ptr
        fail_node = 0
        status = BAND_INVALID_ARGUMENT
        if (n < 1 .or. nj < 3) return
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
        allocate(f)
        call band_factor(n, nj, Af, Bf, Df, f, X=Xf, Y=Yf)
        status = f%status
        fail_node = f%fail_node
        handle = c_loc(f)
    end function bandsolver_f_factor

    integer(c_int) function bandsolver_f_factor_solve(handle, G, dc) bind(c, name='bandsolver_f_factor_solve') &
            result(status)
        type(c_ptr), value :: handle
        real(c_double), intent(in), target :: G(*)
        real(c_double), intent(out), target :: dc(*)
        type(band_factorization), pointer :: f
        real(c_double), pointer :: Gv(:,:), dv(:,:)
        integer :: fstatus
        status = BAND_INVALID_ARGUMENT
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, f)
        if (f%n < 1) then
            status = f%status; return
        end if
        Gv(1:f%n,1:f%nj) => G(1:f%n*f%nj)
        dv(1:f%n,1:f%nj) => dc(1:f%n*f%nj)
        call band_factor_solve(f, Gv, dv, fstatus)
        status = fstatus
    end function bandsolver_f_factor_solve

    subroutine bandsolver_f_factor_free(handle) bind(c, name='bandsolver_f_factor_free')
        type(c_ptr), value :: handle
        type(band_factorization), pointer :: f
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, f)
        deallocate(f)
    end subroutine bandsolver_f_factor_free

end module bandsolver_capi
