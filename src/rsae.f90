!###############################################################################
!
!   rsae: robust small area estimation
!   Copyright (C) 2024 Tobias Schoch
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
!   SUBJECT:        This file, rsae.f90, contains Fortran 90 code for
!                   robust estimation of the basic unit-level small area
!                   model by means of Huber M-estimation
!   AUTHOR:         Tobias Schoch (tobias.schoch@gmail.com)
!   DEPENDENCiES:   LAPACK
!                   BLAS
!                   Richard Brents zeroin.f (GPL; see file 'rsaeext.f90')
!
!===============================================================================
module blas_lapack
    implicit none
    public

    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            integer     :: m, n, nrhs, lda, ldb, lwork, info
            character   :: trans
            double precision, dimension( lda, * )   :: a
            double precision, dimension( ldb, * )   :: b
            double precision, dimension( * )        :: work
        end subroutine

        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            integer     :: m, n, lda, incx, incy
            character   :: trans
            double precision                        :: alpha, beta
            double precision, dimension( lda, * )   :: a
            double precision, dimension( * )        :: x
            double precision, dimension( * )        :: y
        end subroutine

        subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            integer     :: n, k, lda, ldc
            character                               :: uplo, trans
            double precision                        :: alpha, beta
            double precision, dimension( lda, * )   :: a
            double precision, dimension( ldc, * )   :: c
        end subroutine

        subroutine dpotri(uplo, n, a, lda, info)
            integer     :: n, lda, info
            character   uplo
            double precision, dimension( lda, * )   :: a
        end subroutine

        subroutine dpotrf (uplo, n, a, lda, info)
            integer     :: n, lda, info
            character   :: uplo
            double precision, dimension( lda, * )   :: a
        end subroutine

        subroutine dtrmm (side, uplo, transa, diag, m, n, alpha, a, lda, b, &
                          ldb)
            integer     :: m, n, lda, ldb
            character   :: side, uplo, transa, diag
            double precision                        :: alpha
            double precision, dimension( lda, * )   :: a
            double precision, dimension( ldb, * )   :: b
        end subroutine

    end interface
end module blas_lapack

!SUBROUTINE:   dconvumtofull
!DESCRIPTION:  convert a upper triagular matrix to a full matrix
!ON ENTRY:
!   INTEGER n(1)
!   REAL mat(n, n)
!ON RETURN
!   REAL mat(n, n)
!-------------------------------------------------------------------------------
subroutine dconvumtofull(n, mat)
    implicit none
    integer, intent(in):: n
    double precision, intent(inout) :: mat(n, n)
    !local declarations
    integer :: i, j
    do i = 2, n
        do j = 1, (i - 1)
            mat(i, j) = mat(j, i)
        end do
    end do
end subroutine
!===============================================================================
!SUBROUTINE:   drsaebeta
!DESCRIPTION:  compute regression coefficients
!DEPENDENCY:
!   dgemv (BLAS2 and LAPACK), dgels (LAPACK)
!   dhuberwgt, dsqrtinvva
!ON ENTRY:
!   INTEGER n(1), p(1), k(1), nsize(g), info(1), dec(1), decorr(1),
!       lwork_dgels(1)
!   REAL k(1), xmat(n,p) yvec(n), d(1), v(1)
!       beta(p), work_dgels(:)
!ON RETURN:
!   INTEGER info(1)
!   REAL beta(p), sumwgt(1)
!-------------------------------------------------------------------------------
subroutine drsaebeta(n, p, g, lwork_dgels, k, xmat, yvec, work_dgels, v, d, &
                     nsize, beta, sumwgt, info, dec, decorr)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, g, dec, decorr, lwork_dgels
    integer, intent(in) :: nsize(g)
    integer, intent(out) :: info
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(in) :: yvec(n)
    double precision, intent(in) :: d, v, k
    double precision, intent(inout) :: beta(p)
    double precision, intent(inout) :: work_dgels(lwork_dgels)
    double precision, intent(out) :: sumwgt
    !local declarations
    integer :: i, j
    double precision :: modyvec(n), res(n)
    double precision :: modxmat(n, p)
    res = yvec
    call dgemv("N", n, p, -1d0, xmat, n, beta, 1, 1d0, res, 1)
    call dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, decorr, res)
    call dhuberwgt(n, k, 1, res)
    modxmat = xmat
    modyvec = yvec
    call dsqrtinvva(n, p, g, nsize, d, v, 1, dec, 0, modxmat)
    call dsqrtinvva(n, 1, g, nsize, d, v, 1, dec, 0, modyvec)
    do j = 1, p
        sumwgt = 0d0
        do i = 1, n
            modxmat(i, j) = modxmat(i, j) * res(i)
            modyvec(i) = modyvec(i) * res(i)
            sumwgt = sumwgt + res(i)**2
        end do
    end do
    call dgels("N", n, p, 1, modxmat, n, modyvec, n, work_dgels, lwork_dgels, &
               info)
    if (info == 0) then
        beta = modyvec(1:p)
    else
        beta = 0
    end if
end subroutine
!===============================================================================
!SUBROUTINE:   drsaehubpredict
!DESCRIPTION:  robust prediction of random effects
!DEPENDENCY:   dsqrtinvva, huberpsi, dgemv(BLAS)
!ON ENTRY:
!   INTEGER n(1), p(1), g(1), nsize(g), dec(1)
!   REAL k(1), kappa(1), d(1), v(1), beta(p), yvec(n), xmat(n,p)
!ON RETURN
!   REAL predfe(g), predre(g)
!-------------------------------------------------------------------------------
subroutine drsaehubpredict(n, p, g, nsize, k, kappa, d, v, beta, yvec, &
                           xmat, predfe, predre, dec)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, g, dec
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: k, kappa, d, v
    double precision, intent(in) :: beta(p)
    double precision, intent(in) :: yvec(n)
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(out) :: predfe(g), predre(g)
    !local declarations
    integer :: i, j
    integer :: l(g)
    double precision :: workfe, workre, sigma2a
    double precision :: res(n), stdr(n), yhat(n)
    sigma2a = v * d
    call dgemv("N", n, p, 1d0, xmat, n, beta, 1, 0d0, yhat, 1)
    res = yvec
    call dgemv("N", n, p, -1d0, xmat, n, beta, 1, 1d0, res, 1)
    stdr = res
    call dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, 0, stdr)
    call dhuberpsi(n, k, stdr)
    call dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, 0, stdr)
    l(1) = 1
    do i = 2, g
        l(i) = l(i - 1) + nsize(i - 1)
    end do
    do i = 1, g
        workre = 0d0
        workfe = 0d0
        do j = 1, nsize(i)
            workre = workre + stdr(l(i) + j - 1)
            workfe = workfe + yhat(l(i) + j - 1)
        end do
        predre(i) = workre * sigma2a * (1 / kappa)
        predfe(i) = workfe / nsize(i)
    end do
end subroutine
!===============================================================================
!SUBROUTINE:   dsqrtinvva
!STATUS:       June 23, 2011; mod November 3, 2011
!DESCRIPTION:  pre-multiply a matrix (vector) by the square root of
!              the inverse covariance matrix with either
!              decomposition (dec):
!                   0: SVD
!                   1: Choleski
!
!              parametrization
!                   0: MLM
!                   1: Hartley-Rao
!                   (else): returns a zero matrix/vector
!
!              decorr (decorrelation; works only for amat of size 1; i.e,
!              y-variable)
!                   0: [as is]
!                   1: center by the median instead of the mean
!
!BENCHMARK:    (self-contained testing; approved June 23, 2011;
!              modifications: Nov 16, 2011)
!DEPENDENCY:   DPOTRF (LAPACK); dmedmad (sctrob)
!FORTRAN:      uses dynamic allocation (only v90, not v77)
!ON ENTRY:
!   INTEGER n(1), p(1), g(1), nsize(g), par(1), dec(1), decorr(1)
!   REAL d(1), v(1), amat(n,p)
!ON RETURN:
!   REAL amat(n,p)
!-------------------------------------------------------------------------------
subroutine dsqrtinvva(n, p, g, nsize, d, v, par, dec, decorr, amat)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, g, par, dec, decorr
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: d, v
    double precision, intent(inout) :: amat(n, p)
    !local declarations
    integer :: i, j, k, info
    integer :: l(g), u(g)
    double precision :: fd, sqrtv, med
    double precision :: m(p)
    double precision, allocatable :: winvv(:, :)
    l(1) = 1
    u(1) = nsize(1)
    do i = 2, g
        l(i) = l(i - 1) + nsize(i - 1)
        u(i) = u(i - 1) + nsize(i)
    end do
    if (dec == 0) then
        select case (par)
        case(1)  !sqrtinvv of type Hartley-Rao as default
            if(decorr == 1)  then
                !robust decorrelation (i.e., centering by median)
                do i = 1, g
                    fd = (1 / sqrt(1 + d * nsize(i))) - 1
                    call dmedmad(nsize(i), amat(l(i) : u(i), 1), 0, med)
                    do j = 1, nsize(i)
                        amat(l(i) + j - 1, 1) = (fd * med) + amat(l(i) &
                            + j - 1, 1)
                    end do
                end do
            !standard case (no robust centering)
            else
                do i = 1, g
                    fd = (1 / sqrt(1 + d * nsize(i))) - 1
                    m = fd * (sum(amat(l(i):u(i), :), 1) / nsize(i))
                    do k = 1, p
                        do j = 1, nsize(i)
                            amat(l(i) + j - 1, k) = m(k) + &
                                amat(l(i) + j - 1, k)
                        end do
                    end do
                end do
            end if
        case(0)  !default MLM sqrtinv * amat
            if (decorr == 1) then
                sqrtv = sqrt(v)
                do i = 1, g
                    fd = (1 / sqrt(1 + d * nsize(i))) - 1
                    call dmedmad(nsize(i), amat(l(i) : u(i), 1), 0, med)
                    do j = 1, nsize(i)
                        amat(l(i) + j - 1, 1) = ((fd / sqrtv ) * med) + &
                            amat(l(i) + j - 1, 1)
                    end do
                end do
            else
                sqrtv = sqrt(v)
                do i = 1, g
                    fd = (1 / sqrt(1 + d * nsize(i))) - 1
                    m = fd * (sum(amat(l(i):u(i), :), 1) / nsize(i))
                    do k = 1, p
                        do j = 1, nsize(i)
                            amat(l(i) + j - 1, k) = (1 / sqrtv) * (m(k) &
                                + amat(l(i) + j - 1, k))
                        end do
                    end do
                end do
            end if
        case default
            amat = 0d0
        end select
    else
        do i = 1, g
            allocate(winvv(nsize(i), nsize(i)))
            winvv =  (-d) / (1d0 + (d * nsize(i)))
            do j = 1, nsize(i)
                winvv(j, j) = winvv(j, j) + 1
            end do
            call dpotrf("U", nsize(i), winvv, nsize(i), info)
            call dtrmm("L", "U", "N", "N", nsize(i), p, 1d0, winvv, &
                       nsize(i), amat(l(i) : u(i), :), nsize(i))
            deallocate(winvv)
        end do
        if (par == 0) then
            do j = 1, p
                do i = 1, n
                    amat(i, j) = amat(i, j) * sqrt(v)
                end do
            end do
        end if
    end if
end subroutine
!
!===============================================================================
!FUNCTION:     is_converged
!DESCRIPTION:  computes the squared norm for two parameter vectors and
!              evaluates if their difference is smaller than the
!              reference acc; return 1 if true, 0 else.
!              this function computation the termination rule for
!              iterative algorithms (either for parameters or resids)
!DEPENDENCY:   none
!ON ENTRY:
!   INTEGER p(1)
!   REAL acc(1), oldvec(p), newvec(p)
!ON RETURN
!   INTEGER is_converged(1)
!-------------------------------------------------------------------------------
function is_converged(p, oldvec, newvec, acc)
    implicit none
    integer, intent(in) :: p
    integer :: is_converged
    double precision, intent(in) :: acc
    double precision, intent(in) :: oldvec(p), newvec(p)
    !local declarations
    double precision :: ratio
    double precision, parameter :: themin = 1.0d-15
    is_converged = 0
    ratio = sqrt(sum((oldvec - newvec)**2) / max(sum(oldvec**2), themin))
    if (ratio < acc) then
        is_converged = 1
    end if
end function is_converged
!===============================================================================
!SUBROUTINE:   drsaebetaiter
!DESCRIPTION:  fully iterated drsaebeta; info carries the # of iter
!DEPENDENCY:
!   drsaebeta
!   is_converged
!ON ENTRY:
!   INTEGER n(1), p(1), k(1), nsize(g), iter(1), dec(1), decorr(1),
!       lwork_dgels(1)
!   REAL k(1), xmat(n,p) yvec(n), v(1), d(1), acc(1)
!       beta(p), work_dgels(:)
!ON RETURN:
!   INTEGER converged(1), info(1)
!   REAL beta(p), sumwgt(1)
!-------------------------------------------------------------------------------
subroutine drsaebetaiter(n, p, g, lwork_dgels, k, xmat, yvec, work_dgels, v, &
                         d, nsize, acc, beta, iter, converged, sumwgt, info, &
                         dec, decorr)
    implicit none
    integer, intent(in) :: n, p, g, dec, decorr, lwork_dgels
    integer, intent(in) :: nsize(g)
    integer, intent(in) :: iter
    integer, intent(out) :: converged
    integer, intent(out) :: info
    double precision, intent(in) :: acc, v, d, k
    double precision, intent(in) :: yvec(n)
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(inout) :: beta(p), work_dgels(lwork_dgels)
    double precision, intent(out) :: sumwgt
    !local declarations
    integer :: i, coinfo, niter, is_converged
    double precision :: betaold(p)
    !iterative updating
    niter = 0
    do i = 1, iter
        betaold = beta
        call drsaebeta(n, p, g, lwork_dgels, k, xmat, yvec, work_dgels, v, d, &
                       nsize, beta, sumwgt, coinfo, dec, decorr)
        if (coinfo /= 0) then
            beta = 0
            exit
        end if
        niter = niter + 1
        converged = is_converged(p, betaold, beta, acc)
        if (converged == 1) then
            exit
        end if
    end do
    info = niter
end subroutine
!===============================================================================
!SUBROUTINE:   drsaehubvariance
!DESCRIPTION:  robust prediction of random effects
!DEPENDENCY:   dsqrtinvva, dsyrk(BLAS)
!ON ENTRY:
!  INTEGER
!  REAL
!ON RETURN
!  REAL
!-------------------------------------------------------------------------------
subroutine drsaehubvariance(n, p, g, nsize, v, d, xmat, vcovbeta, dec)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, g, dec
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: d, v
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(out) :: vcovbeta(p, p)
    !local declarations
    integer :: info
    double precision :: modx(n, p)
    double precision :: mxtmx(p, p), fmxtmx(p, p)
    vcovbeta = 0d0
    mxtmx = 0d0
    modx = xmat
    call dsqrtinvva(n, p, g, nsize, d, v, 0, dec, 0, modx)
    call dsyrk("U", "T", p, n, 1d0, modx, n, 0d0, mxtmx, p)
    fmxtmx = mxtmx
    call dconvumtofull(p, fmxtmx)
    call dpotrf("U", p, fmxtmx, p, info)
    if (info == 0) then
        call dpotri("U", p, fmxtmx, p, info)
        if (info == 0) then
            vcovbeta = fmxtmx
        else
            vcovbeta = info * 1d0
        end if
    else
        vcovbeta = info * 1d0
    end if
end subroutine
!===============================================================================
!SUBROUTINE:   dhuberpsi
!DESCRIPTION:  compute huber psi
!BENCHMARK:    robustbase:::huberPsi@psi (v 0.7-6),
!              approved June 18, 2011
!DEPENDENCY:   none
!ON ENTRY:
!   INTEGER n(1)
!   REAL k(1), vec(n)
!ON RETURN
!   REAL vec(n)
!-------------------------------------------------------------------------------
subroutine dhuberpsi(n, k, vec)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: k
    double precision, intent(inout) :: vec(n)
    !local declarations
    integer :: i
    double precision :: absvec(n)
    absvec = abs(vec)
    do i = 1, n
        if (absvec(i) >= k) then
            vec(i) = sign(k, vec(i))
        end if
    end do
end subroutine
!===============================================================================
!SUBROUTINE:   drsaehubdest
!DESCRIPTION:  evaluates the (non-linear) estimating of d; it is called
!              from the modified zeroin function; kappa is the
!              consistency correction term
!BENCHMARK:
!DEPENDENCY:
!   dhuberpsi, dsqrtinvva
!ON ENTRY:
!   INTEGER n(1), g(1), nsize(g), dec(1), decorr(1)
!   REAL d(1), v(1), k(1), kappa(1), res(n)
!ON RETURN: REAL eval(1)
!-------------------------------------------------------------------------------
subroutine drsaehubdest(n, g, nsize, d, v, k, kappa, res, eval, dec, &
                        decorr)
    implicit none
    integer, intent(in) :: n, g, dec, decorr
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: d, v, k, kappa
    double precision, intent(in) :: res(n)
    double precision, intent(out) :: eval
    !local declarations
    integer :: i, j
    integer :: l(g)
    double precision :: lhs, rhs, work
    double precision :: vec(n)
    !
    vec = res  !needed to create a new res-object because gfortran
    call dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, decorr, vec)
    call dhuberpsi(n, k, vec)
    l(1) = 1
    do i = 2, g
        l(i) = l(i - 1) + nsize(i - 1)
    end do
    lhs = 0d0
    rhs = 0d0
    do i = 1, g
        lhs = lhs + nsize(i) / (1 + d * nsize(i))
        work = 0d0
        do j = 1, nsize(i)
            work = work + vec(l(i) + j - 1) * sqrt(1 / (1 + d * nsize(i)))
        end do
        rhs = rhs + (work**2) / kappa
    end do
    eval = lhs - rhs
end subroutine
!===============================================================================
!SUBROUTINE:   drsaehubvest
!DESCRIPTION:  compute a Huber type M-scale^2 estimate; note that
!              stdres = v^(-1/2)*res, with res=y-x*beta and v is
!              parametrized in Hartley-Rao manner
!              notice, kappa denotes the consistency correction term
!              for Huber's Proposal 2 estimator.
!BENCHMARK:    MASS:::hubers(x, k, mu=0) (v 7.3-12),
!              approved on June 23, 2011
!DEPENDENCY:
!   dhuberwgt
!ON ENTRY:
!   INTEGER n(1), niter(1)
!   REAL k(1), acc(1), kappa(1), stdres(n)
!ON RETURN:
!   INTEGER info(1)
!   REAL v(1), sumwgt(1)
!-------------------------------------------------------------------------------
subroutine drsaehubvest(n, niter, v, k, acc, kappa, stdres, sumwgt, info)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: niter
    integer, intent(out) :: info
    double precision, intent(in) :: k, acc, kappa
    double precision, intent(in) :: stdres(n)
    double precision, intent(inout) :: v
    double precision, intent(out) :: sumwgt
    !local declarations
    integer :: i, iter
    double precision :: ssq, vold
    double precision :: workresid(n)
    vold = v
    do iter = 1, niter
        ssq = 0d0
        sumwgt = 0d0
        workresid = stdres / sqrt(vold)
        call dhuberwgt(n, k, 2, workresid)
        do i = 1, n
            ssq = ssq + workresid(i) * stdres(i)**2
            sumwgt = sumwgt + workresid(i)
        end do
        v = ssq / (n * kappa)
        if (abs(v / vold - 1d0) < acc) then
            exit
        else
            vold = v
        end if
    end do
    info = iter
end subroutine
!===============================================================================
!SUBROUTINE:   dhuberwgt
!DESCRIPTION:  compute huber psi-weight; NOTE:
!              typ = 1 for sqrt(wgt)
!              typ = 0 for wgt,
!              typ = 2 for wgt^2
!BENCHMARK:    robustbase:::huberPsi@wgt (v 0.7-6),
!              approved, June 19, 2011
!DEPENDENCY:   none
!ON ENTRY:
!   INTEGER n(1), typ(1)
!   REAL k(1), vec(n)
!ON RETURN:
!   REAL vec(n)
!-------------------------------------------------------------------------------
subroutine dhuberwgt(n, k, typ, vec)
    implicit none
    integer, intent(in) :: n, typ
    double precision, intent(in) :: k
    double precision, intent(inout) :: vec(n)
    !local declarations
    integer :: i
    double precision :: choice(n)
    choice = k / abs(vec)
    select case (typ)
        case(1) !take the square root of the weights
            do i = 1, n
                if (choice(i) < 1d0) then
                    vec(i) = sqrt(choice(i))
                else
                    vec(i) = 1
                end if
            end do
        case(0) !take the weights as they are
            do i = 1, n
                if (choice(i) < 1d0) then
                    vec(i) = choice(i)
                else
                    vec(i) = 1
                end if
            end do
        case(2) !the weights to the power of two
            do i = 1, n
                if (choice(i) < 1d0) then
                    vec(i) = choice(i)**2
                else
                    vec(i) = 1
                end if
            end do
        case default !an errorneous call returns one
            vec = 0
    end select
end subroutine
!===============================================================================
!SUBROUTINE:      drsaehub
!DESCRIPTION:     drsaehub is the main or workhorse function that
!                 can be called from R.
!                 NOTE that neither this nor any function called
!                 subsequently, i.e., from within rsaehub, does any
!                 tests (except standard data-type-conformity tests)
!                 whether the delivered args match or whether they
!                 are meaningful. This is the job of the calling
!                 instance, namely the R-function rsaehub.
!DEPENDENCY:      drsaebetaiter, dgemv(BLAS), dsqrtinvva, drsaehubvest
!                 drsaehubdestiter, is_converged
!ON ENTRY:
!   INTEGER  n(1), p(1), g(1), niter(1), nsize(g), iter(2), dec(1),
!            decorr(1)
!   REAL     k(2), allacc(1), epsd(1), acc(3), sumwgt(3), xmat(n, p),
!            yvec(n), kappa(2), iterrecord(niter,p+2),
!ON RETURN:
!   INTEGER: niter(1), allconverged(1)
!   REAL:    tau(p+2), taurecord(niter, p+2)
!-------------------------------------------------------------------------------
subroutine drsaehub(n, p, g, niter, nsize, iter, iterrecord, allacc, &
                    acc, sumwgt, xmat, yvec, k, kappa, epsd, tau, taurecord, &
                    allconverged, dec, decorr)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, g, dec, decorr
    integer, intent(in) :: niter
    integer, intent(in) :: nsize(g)
    integer, intent(in) :: iter(2)
    integer, intent(out):: allconverged
    double precision, intent(in) :: allacc, epsd
    double precision, intent(in) :: k(3)
    double precision, intent(in) :: kappa(2)
    double precision, intent(in) :: acc(3)
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(in) :: yvec(n)
    double precision, intent(out) :: iterrecord(niter, 3)
    double precision, intent(inout) :: tau(p + 2)
    double precision, intent(out) :: taurecord(niter, p + 2)
    double precision, intent(out) :: sumwgt(3)
    !local declarations
    integer :: i, j, info, lwork_dgels, convergedbeta, work, monitord, &
        is_converged
    double precision :: upper
    double precision :: oldtau(p + 2), prep_dgels(2)
    double precision :: res(n), stdres(n), sumwgtres(n)
    ! dynamically allocate work_array (used in dgels)
    double precision, allocatable :: work_dgels(:)
    lwork_dgels = -1
    call dgels("N", n, p, 1, xmat, n, yvec, n, prep_dgels, lwork_dgels, info)
    ! size of the dgels work array
    lwork_dgels = int(prep_dgels(1))
    ! allocate array and check
    if (info == 0) then
        allocate(work_dgels(lwork_dgels))
    end if
    if (.not. allocated(work_dgels)) then
        return
    end if
    ! initialize iteration
    iterrecord = 0
    allconverged = 0
    monitord = 0
    do i = 1, niter
        oldtau(1 : p) = tau(1 : p)
        oldtau(p + 1) = tau(p + 1)
        oldtau(p + 2) = tau(p + 2)
        ! compute regression coefficients
        call drsaebetaiter(n, p, g, lwork_dgels, k(1), xmat, yvec, work_dgels, &
                           tau(p + 1), tau(p + 2), nsize, acc(1), tau(1:p), &
                           iter(1), convergedbeta, sumwgt(1), work, dec, decorr)
        iterrecord(i, 1) = work
        if (convergedbeta /= 1) then
            iterrecord(1, i) = (-1) * iterrecord(1, i)
        end if
        ! residuals
        res = yvec
        call dgemv("N", n, p, -1d0, xmat, n, tau(1:p), 1, 1d0, res, 1)
        stdres = res
        ! std. residuals
        call dsqrtinvva(n, 1, g, nsize, tau(p + 2), tau(p + 1), 1, dec, &
                        decorr, stdres)
        ! variance component v (unit-level errors)
        call drsaehubvest(n, iter(2), tau(p + 1), k(2), acc(2), kappa(1), &
                          stdres, sumwgt(2), work)
        iterrecord(i, 2) = work
        if (monitord == 1) then
            tau(p + 2) = 0d0
            iterrecord(i, 3) = 0d0
        else
            upper = tau(p + 2) * 1d1
            ! ratio of variance components
            call drsaehubdestiter(n, g, nsize, tau(p + 1), k(3), kappa(2), &
                                  res, 0d0, upper, acc(3), tau(p + 2), work, &
                                  dec, decorr)
            iterrecord(i, 3) = work
            if (sum(taurecord(max(i - 2, 1):i, p + 2)) < 3 * epsd .and. &
                    i >= 3) then
                monitord = 1
            end if
        end if
        taurecord(i, :) = tau
        ! check for convergence
        allconverged = is_converged(p + 1, oldtau, tau, allacc)
        if (allconverged == 1) then
            exit
        end if
    end do
    sumwgtres = res
    call dsqrtinvva(n, 1, g, nsize, tau(p + 2), tau(p + 1), 0, dec, &
                    decorr, sumwgtres)
    call dhuberwgt(n, k(3), 0, sumwgtres)
    sumwgt(3) = 0d0
    do j = 1, n
        sumwgt(3) = sumwgt(3) + sumwgtres(j)
    end do
    deallocate(work_dgels)
end subroutine
!===============================================================================
!SUBROUTINE:   drsaeresid
!DESCRIPTION:  get the residuals and the huber wgt; i.e., res = e_ij
!              = y_ij - X_ij * beta - u_i; and stdres=V^(-1/2)*res;
!              the huber weight is w.r.t. to the model-psi (not the
!              prediction)
!BENCHMARK:
!DEPENDENCY:   dgemv(BLAS), dsqrtinvva, huberwgt
!ON ENTRY:
!   INTEGER n(1), p(1), g(1), nsize(g), dec(1)
!   REAL k(1), tau(p+2), xmat(n,p) yvec(n), u(g)
!ON RETURN:
!   REAL res(n), stdres(n), wgt(n)
!-------------------------------------------------------------------------------
subroutine drsaeresid(n, p, g, nsize, k, tau, u, xmat, yvec, res, stdres, &
                      wgt, dec)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, g, dec
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: k !robustness tuning constant
    double precision, intent(in) :: tau(p + 2) ! (beta, v, d)
    double precision, intent(in) :: u(g) !vector of random effects
    double precision, intent(in) :: yvec(n)
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(out) :: res(n), stdres(n), wgt(n)
    !local declarations
    integer :: i, j
    integer :: l(g)
    res = yvec
    call dgemv("N", n, p, -1d0, xmat, n, tau(1:p), 1, 1d0, res, 1)
    l(1) = 1
    do i = 2, g
        l(i) = l(i - 1) + nsize(i - 1)
    end do
    do i = 1, g
        do j = 1, nsize(i)
            res(l(i) + j - 1) = res(l(i) + j - 1) - u(i)
        end do
    end do
    stdres = res
    call dsqrtinvva(n, 1, g, nsize, tau(p + 2), tau(p + 1), 0, dec, 0, stdres)
    wgt = stdres
    call dhuberwgt(n, k, 0, wgt)
end subroutine
!===============================================================================
!SUBROUTINE:   dmedmed
!AUTHOR:       Tobias Schoch, November 19, 2011
!DESCRIPTION:  computes the median (type=2 in R; i.e., average over
!              discontinuities) and compute the mad
!              typ = 0; median
!              typ = 1; mad
!DEPENDENCY:   it uses qsort3 from R's utility functions (for C, it is
!              in the header file R_ext/Utils.h (included in R.h)
!ON ENTRY:
!   INTEGER n(1), typ(1)
!   REAL vx(n)
!ON RETURN
!   REAL res(1)
!-------------------------------------------------------------------------------
subroutine dmedmad(n, vx, typ, res)
    implicit none
    integer, intent(in):: n, typ
    double precision, intent(in) :: vx(n)
    double precision, intent(out) :: res
    !local declarations
    integer :: i
    double precision :: work(n)
    work = vx
    call qsort3(work, 1, n)
    if (modulo(n, 2) == 0) then
        res = (work(n / 2) + work((n / 2) + 1)) / 2d0
    else
        res = work(((n - 1) / 2) + 1)
    end if
    if (typ == 1) then
        do i = 1, n
            work(i) = abs(work(i) - res)
        end do
        call qsort3(work, 1, n)
        if (modulo(n, 2) == 0) then
            res = (1.4814 / 2d0) * (work(n / 2) + work((n / 2) + 1))
        else
            res = 1.4814 * work(((n - 1) / 2) + 1)
        end if
    end if
end subroutine
!===============================================================================
!SUBROUTINE:   drlm
!AUTHOR:       Tobias Schoch, November 20, 2011
!DESCRIPTION:  compute regression M-estimates, using the mad as the
!              preliminary scale estimate
!DEPENDENCY:
!   dgemv (BLAS2 and LAPACK), dgels (LAPACK)
!   dhuberwgt, dmedmad
!ON ENTRY:
!   INTEGER n(1), p(1),
!   REAL xmat(n, p), yvec(n), k(1), beta(p), acc(1)
!ON RETURN:
!   INTEGER info(1)
!   REAL beta(p), s(1)
!-------------------------------------------------------------------------------
subroutine drlm(n, p, xmat, yvec, k, beta, s, info, niter, acc)
    use blas_lapack
    implicit none
    integer, intent(in) :: n, p, niter
    integer, intent(out) :: info
    double precision, intent(in) :: k, acc
    double precision, intent(in) :: yvec(n)
    double precision, intent(in) :: xmat(n, p)
    double precision, intent(out) :: s
    double precision, intent(inout) :: beta(p)
    ! local declarations (most are used in dqrls)
    integer :: i, j, l, lwork_dgels, converged, is_converged
    double precision :: prep_dgels(2)
    double precision :: oldbeta(p)
    double precision :: modyvec(n), res(n)
    double precision :: modxmat(n, p)
    ! dynamically allocate work_array (used in dgels)
    double precision, allocatable :: work_dgels(:)
    lwork_dgels = -1
    call dgels("N", n, p, 1, xmat, n, yvec, n, prep_dgels, lwork_dgels, info)
    ! size of the dgels work array
    lwork_dgels = int(prep_dgels(1))
    ! allocate array and check
    if (info == 0) then
        allocate(work_dgels(lwork_dgels))
    end if
    if (.not. allocated(work_dgels)) then
        return
    end if
    ! start iteration
    do l = 1, niter
        oldbeta = beta
        res = yvec
        call dgemv("N", n, p, -1d0, xmat, n, oldbeta, 1, 1d0, res, 1)
        call dmedmad(n, res, 1, s)
        do i = 1, n
            res(i) = res(i) / s
        end do
        call dhuberwgt(n, k, 1, res)
        do j = 1, p
            do i = 1, n
                modyvec(i) = yvec(i) * res(i)
                modxmat(i, j) = xmat(i, j) * res(i)
            end do
        end do
        call dgels("N", n, p, 1, modxmat, n, modyvec, n, work_dgels, &
                   lwork_dgels, info)
        if (info == 0) then
            beta = modyvec(1:p)
        else
            beta = 0d0
        end if
        converged = is_converged(p, oldbeta, beta, acc)
        if (converged == 1) then
            info = l
            exit
        end if
    end do
    deallocate(work_dgels)
end subroutine
!===============================================================================
!SUBROUTINE:   drsaehubdestiter
!PART OF:      rsae
!DESCRIPTION:  root-finding device; evaluates drsaehubdest between
!              lower and upper; there must be a sign change in the
!              interval, otherwise the sub stops
!DEPENDENCY:   zero_rc
!ON ENTRY:
!   INTEGER n(1), g(1), nsize(g), dec(1)
!   REAL v(1), k(1), kappa(1), lower(1), upper(1), tol(1), res(n)
!ON RETURN
!   INTEGER info(1)
!   REAL zeroin(1)
!-------------------------------------------------------------------------------
subroutine drsaehubdestiter(n, g, nsize, v, k, kappa, res, lower, upper, &
                            tol, zeroin, info, dec, decorr)
    implicit none
    !zeroin declarations
    integer, intent(out) :: info
    double precision, intent(in) :: lower, upper, tol
    double precision, intent(out) :: zeroin
    integer, intent(in) :: n, g, dec, decorr
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: v, k, kappa
    double precision, intent(in) :: res(n)
    !local declarations
    integer, parameter :: ITMAX = 100
    integer :: iter, state
    double precision :: f_lower, f_upper, arg, f_value

    !evaluate the function at the interval borders (lower, upper)
    call drsaehubdest(n, g, nsize, lower, v, k, kappa, res, f_lower, dec, &
                      decorr)
    call drsaehubdest(n, g, nsize, upper, v, k, kappa, res, f_upper, dec, &
                      decorr)
    !stop if there is no sign change of the functions values over the interval
    if ((f_lower > 0d0 .and. f_upper > 0d0) .or. &
            (f_lower < 0d0 .and. f_upper < 0d0)) then
        info = -1
        zeroin = 0d0
        return
    end if

    !start with Brent's root finding algorithm
    info = 0
    state = 0               !this initializes the search with zero_rc
    arg = 0d0               !no need to initialize function argument
    f_value = 1d0           !matters only in the (iter + 1) call
    do iter = 1, ITMAX
        !check for zero
        call zero_rc(lower, upper, tol, arg, state, f_value)
        !termination criterion
        if (state == 0) then
            info = iter
            exit
        end if
        !evaluate function at arg
        call drsaehubdest(n, g, nsize, arg, v, k, kappa, res, f_value, &
                          dec, decorr)
    end do
    zeroin = arg
end subroutine
