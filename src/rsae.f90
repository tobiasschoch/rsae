!####################################################################
!
!     Copyright (c) 2011, Tobias Schoch
!     All rights reserved.
!
!     Redistribution and use in source and binary forms, with or
!     without modification, are permitted provided that the following
!     conditions are met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!
!     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
!     CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
!     INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
!     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!     DISCLAIMED. IN NO EVENT SHALL Tobias Schoch BE LIABLE FOR
!     ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!     CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
!     OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
!     OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!     LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
!     THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!     SUCH DAMAGE.
!
!####################################################################
!
!SUBJECT:      This file, rsae.f90, contains Fortran 90 code for
!              robust estimation of the basic unit-level small area
!              model by means of Huber M-estimation
!NOTES:        These subroutines are designed to be called from R; see
!              R Development Core Team (2011) for more details on R.
!              Do not call any routine by yourself unless you are
!              certain about the implementation details.
!COMPILER:     GNU gfortran v 4.5.1
!TESTS:        code tested on platform x86_64 SUSE Linux v 11.4
!FORTRAN SPEC: The layout of this file has been chosen to conform
!              with F77's limitation of 72 char per line
!              (otherwise gfortran must be called
!              with flag: -ffree-form; or with file ext. ".f90"
!AUTHOR:       Tobias Schoch, June 12, 2011
!LICENSE:      BSD 2 (i.e., modified, 2-clause Berkeley Software
!              Distribution License, aka Simplified BSD License, aka
!              FreeBSD License; see above)
!DEPENDENCY:   LAPACK
!              BLAS
!              Richard Brents zeroin.f (GPL; see file 'rsaeext.f90')
!
!####################################################################
!
!====================================================================
!SUBROUTINE:   dconvumtofull
!PART OF:      sct
!DESCRIPTION:  convert a upper triagular matrix to a full matrix
!ON ENTRY:
!  INTEGER n(1)
!  REAL mat(n, n)
!ON RETURN
!  REAL mat(n, n)
!--------------------------------------------------------------------
SUBROUTINE dconvumtofull(n, mat)
   IMPLICIT NONE
   INTEGER, INTENT(IN):: n
   DOUBLE PRECISION, INTENT(INOUT) :: mat(n, n)
   !local declarations
   INTEGER :: i, j
   !
   DO i = 2, n
      DO j = 1, (i - 1)
         mat(i, j) = mat(j, i)
      END DO
   END DO
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   drsaebeta
!PART OF:      rsaehuber
!DESCRIPTION:  compute the residual vector
!DEPENDENCY:
!  dgemv (BLAS2 and LAPACK), dgels (LAPACK)
!  dhuberwgt, dsqrtinvva
!ON ENTRY:
!  INTEGER n(1), p(1), k(1), nsize(g), info(1), dec(1), decorr(1)
!  REAL k(1), xmat(n,p) yvec(n), d(1), v(1)
!      beta(p)
!ON RETURN:
!  INTEGER info(1)
!  REAL beta(p), sumwgt(1)
!--------------------------------------------------------------------
SUBROUTINE drsaebeta(n, p, g, k, xmat, yvec, v, d, nsize, beta,&
      sumwgt, info, dec, decorr)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec, decorr
   INTEGER, INTENT(IN) :: nsize(g)
   INTEGER, INTENT(OUT) :: info
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: d, v, k
   DOUBLE PRECISION, INTENT(INOUT) :: beta(p)
   DOUBLE PRECISION, INTENT(OUT) :: sumwgt
   !local declarations (most are used in dqrls)
   INTEGER :: i, j, lwork
   INTEGER, PARAMETER :: lworkmax = 10000
   DOUBLE PRECISION :: work(lworkmax)
   DOUBLE PRECISION :: modyvec(n), res(n)
   DOUBLE PRECISION :: modxmat(n, p)
   !
   res = yvec
   CALL dgemv("N", n, p, -1D0, xmat, n, beta, 1, 1D0, res, 1)
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, decorr, res)
   CALL dhuberwgt(n, k, 1, res)
   modxmat = xmat
   modyvec = yvec
   CALL dsqrtinvva(n, p, g, nsize, d, v, 1, dec, 0, modxmat)
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 1, dec, 0, modyvec)
   DO j = 1, p
      sumwgt = 0D0
      DO i = 1, n
         modxmat(i, j) = modxmat(i, j) * res(i)
         modyvec(i) = modyvec(i) * res(i)
         sumwgt = sumwgt + res(i) ** 2
       END DO
   END DO
   lwork = -1
   CALL dgels("n", n, p, 1, modxmat, n, modyvec, n, work, lwork, info)
   lwork = MIN(lworkmax, INT(work(1)))
   CALL dgels("n", n, p, 1, modxmat, n, modyvec, n, work, lwork, info)
   IF (info == 0) THEN
      beta = modyvec(1:p)
   ELSE
      beta = 0
   END IF
END SUBROUTINE

!====================================================================
!SUBROUTINE:   drsaehubpredict
!PART OF:      rsaehub
!DESCRIPTION:  robust prediction of random effects
!DEPENDENCY:   dsqrtinvva, huberpsi, dgemv(BLAS)
!ON ENTRY:
!  INTEGER n(1), p(1), g(1), nsize(g), dec(1)
!  REAL k(1), kappa(1), d(1), v(1), beta(p), yvec(n), xmat(n,p)
!ON RETURN
!  REAL predfe(g), predre(g)
!--------------------------------------------------------------------
SUBROUTINE drsaehubpredict(n, p, g, nsize, k, kappa, d, v, beta, yvec, &
      xmat, predfe, predre, dec)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: k, kappa, d, v
   DOUBLE PRECISION, INTENT(IN) :: beta(p)
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(OUT) :: predfe(g), predre(g)
   !local declarations
   INTEGER :: i, j
   INTEGER :: l(g)
   DOUBLE PRECISION :: workfe, workre, sigma2a
   DOUBLE PRECISION :: res(n), stdr(n), yhat(n)
   sigma2a = v * d
   CALL dgemv("N", n, p, 1D0, xmat, n, beta, 1, 0D0, yhat, 1)
   res = yvec
   CALL dgemv("N", n, p, -1D0, xmat, n, beta, 1, 1D0, res, 1)
   stdr = res
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, 0, stdr)
   CALL dhuberpsi(n, k, stdr)
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, 0, stdr)
   l(1) = 1
   DO i = 2, g
      l(i) = l(i - 1) + nsize(i - 1)
   END DO
   DO i = 1, g
      workre = 0D0
      workfe = 0D0
      DO j = 1, nsize(i)
         workre = workre + stdr(l(i) + j - 1)
         workfe = workfe + yhat(l(i) + j - 1)
      END DO
      predre(i) = workre * sigma2a * (1 / kappa)
      predfe(i) = workfe / nsize(i)
   END DO
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   dsqrtinvva
!PART OF:      rsaehuber
!STATUS:       June 23, 2011; mod November 3, 2011
!DESCRIPTION:  pre-multiply a matrix (vector) by the square root of
!              the inverse covariance matrix with either
!              decomposition (dec):
!                 0: SVD
!                 1: Choleski
!
!              parametrization
!                 0: MLM
!                 1: Hartley-Rao
!                 (else): returns a zero matrix/vector
!
!              decorr (decorrelation; works only for amat of size 1;
!                 i.e, y-variable)
!                 0: [as is]
!                 1: center by the median instead of the mean
!
!BENCHMARK:    (self-contained testing; approved June 23, 2011;
!              modifications: Nov 16, 2011)
!DEPENDENCY:   DPOTRF (LAPACK); dmedmad (sctrob)
!FORTRAN:      uses dynamic allocation (only v90, not v77)
!ON ENTRY:
!  INTEGER n(1), p(1), g(1), nsize(g), par(1), dec(1), decorr(1)
!  REAL d(1), v(1), amat(n,p)
!ON RETURN:
!  REAL amat(n,p)
!--------------------------------------------------------------------
SUBROUTINE dsqrtinvva(n, p, g, nsize, d, v, par, dec, decorr, amat)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, par, dec, decorr
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: d, v
   DOUBLE PRECISION, INTENT(INOUT) :: amat(n, p)
   INTEGER :: i, j, k, info
   INTEGER :: l(g), u(g)
   DOUBLE PRECISION :: fd, sqrtv, med
   DOUBLE PRECISION :: m(p)
   DOUBLE PRECISION, ALLOCATABLE :: winvv(:, :)
   l(1) = 1
   u(1) = nsize(1)
   DO i = 2, g
      l(i) = l(i - 1) + nsize(i - 1)
      u(i) = u(i - 1) + nsize(i)
   END DO
   IF (dec == 0) THEN
      SELECT CASE (par)
      CASE(1)  !sqrtinvv of type Hartley-Rao as default
         IF(decorr == 1)  THEN
         !robust decorrelation (i.e., centering by median)
            DO i = 1, g
               fd = (1 / SQRT(1 + d * nsize(i))) - 1
               CALL dmedmad(nsize(i), amat( l(i) : u(i), 1), 0, med)
               DO j = 1, nsize(i)
                  amat(l(i) + j - 1, 1) = (fd * med) + amat(l(i) + j - 1, 1)
               END DO
            END DO
         !standard case (no robust centering)
         ELSE
            DO i = 1, g
               fd = (1 / SQRT(1 + d * nsize(i))) - 1
               m = fd * (SUM(amat(l(i):u(i), :), 1) / nsize(i))
               DO k = 1, p
                  DO j = 1, nsize(i)
                     amat(l(i) + j - 1, k) = m(k) + amat(l(i) + j - 1, k)
                  END DO
               END DO
            END DO
         END IF
      CASE(0)  !default MLM sqrtinv * amat
         IF (decorr == 1) THEN
            sqrtv = SQRT(v)
            DO i = 1, g
               fd = (1 / SQRT(1 + d * nsize(i))) - 1
               CALL dmedmad(nsize(i), amat( l(i) : u(i), 1), 0, med)
               DO j = 1, nsize(i)
                  amat(l(i) + j - 1, 1) = (( fd / sqrtv ) * med) + &
                     amat(l(i) + j - 1, 1)
               END DO
            END DO
         ELSE
            sqrtv = SQRT(v)
            DO i = 1, g
               fd = (1 / SQRT(1 + d * nsize(i))) - 1
               m = fd * (SUM(amat(l(i):u(i), :), 1) / nsize(i))
               DO k = 1, p
                  DO j = 1, nsize(i)
                     amat(l(i) + j - 1, k) = (1/sqrtv) * (m(k) &
                     + amat(l(i) + j - 1, k))
                  END DO
               END DO
            END DO
         END IF
      CASE DEFAULT
         amat = 0D0
      END SELECT
   ELSE
      DO i = 1, g
         ALLOCATE( winvv( nsize(i), nsize(i) ) )
         winvv =  (- d) / ( 1D0 + ( d * nsize(i) ) )
         DO j = 1, nsize(i)
            winvv(j, j) = winvv(j, j) + 1
         END DO
         CALL dpotrf("U", nsize(i), winvv, nsize(i), info)
         CALL dtrmm("L", "U", "N", "N", nsize(i), p, 1D0, winvv, &
            nsize(i), amat( l(i) : u(i), :), nsize(i))
         DEALLOCATE(winvv)
      END DO
      IF (par == 0) THEN
         DO j = 1, p
            DO i = 1, n
               amat(i, j) = amat(i, j) * SQRT(v)
            END DO
         END DO
      END IF
   END IF
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   ddelta
!PART OF:      sctrob
!DESCRIPTION:  computes the squared norm for two parameter vectors and
!              evaluates if their difference is smaller than the
!              reference acc; return 1 if true, 0 else.
!              this function computation the termination rule for
!              iterative algorithms (either for parameters or resids)
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER p(1), res(1)
!  REAL acc(1), oldvec(p), newvec(p)
!ON RETURN
!  INTEGER res(1)
!--------------------------------------------------------------------
SUBROUTINE ddelta(p, oldvec, newvec, acc, res)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: p
   INTEGER, INTENT(OUT) :: res
   DOUBLE PRECISION, INTENT(IN) :: acc
   DOUBLE PRECISION, INTENT(IN) :: oldvec(p), newvec(p)
   !local declaration
   DOUBLE PRECISION :: ratio
   DOUBLE PRECISION, PARAMETER :: themin = 1.0D-15
   !
   ratio = SQRT(SUM((oldvec - newvec)**2) / MAX(SUM(oldvec**2), themin))
   IF (ratio < acc) THEN
      res = 1
   ELSE
      res = 0
   END IF
END SUBROUTINE
!

!====================================================================
!SUBROUTINE:   drsaebetaiter
!PART OF:      rsaehuber
!DESCRIPTION:  fully iterated drsaebeta; info carries the # of iter
!DEPENDENCY:
!  drsaebeta
!  ddelta
!ON ENTRY:
!  INTEGER n(1), p(1), k(1), nsize(g), iter(1), dec(1), decorr(1)
!  REAL k(1), xmat(n,p) yvec(n), v(1), d(1), acc(1)
!      beta(p)
!ON RETURN:
!  INTEGER converged(1), info(1)
!  REAL beta(p), sumwgt(1)
!--------------------------------------------------------------------
SUBROUTINE drsaebetaiter(n, p, g, k, xmat, yvec, v, d, nsize, acc, &
      beta, iter, converged, sumwgt, info, dec, decorr)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec, decorr
   INTEGER, INTENT(IN) :: nsize(g)
   INTEGER, INTENT(IN) :: iter
   INTEGER, INTENT(OUT) :: converged
   INTEGER, INTENT(OUT) :: info
   DOUBLE PRECISION, INTENT(IN) :: acc, v, d, k
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(INOUT) :: beta(p)
   DOUBLE PRECISION, INTENT(OUT) :: sumwgt
   INTEGER :: i, coinfo, niter
   DOUBLE PRECISION :: betaold(p)
   !
   niter = 0
   DO i = 1, iter
      betaold = beta
      CALL drsaebeta(n, p, g, k, xmat, yvec, v, d, nsize,&
         beta, sumwgt, coinfo, dec, decorr)
      IF (coinfo /= 0) THEN
         beta = 0
         EXIT
      END IF
      niter = niter + 1
      CALL ddelta(p, betaold, beta, acc, converged)
      IF (converged == 1) THEN
         EXIT
      END IF
   END DO
   info = niter
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   drsaehubvariance
!PART OF:      rsaehub
!DESCRIPTION:  robust prediction of random effects
!DEPENDENCY:   dsqrtinvva, dsyrk(BLAS), dtrtri(LAPACK)
!ON ENTRY:
!  INTEGER
!  REAL
!ON RETURN
!  REAL
!--------------------------------------------------------------------
SUBROUTINE drsaehubvariance(n, p, g, nsize, kappa, v, d, xmat, &
      vcovbeta, dec)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: kappa, d, v
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(OUT) :: vcovbeta(p, p)
   !local declarations
   INTEGER :: info
   DOUBLE PRECISION :: modx(n, p)
   DOUBLE PRECISION :: mxtmx(p, p), fmxtmx(p, p)
   vcovbeta = 0D0
   mxtmx = 0D0
   modx = xmat
   CALL dsqrtinvva(n, p, g, nsize, d, v, 0, dec, 0, modx)
   CALL dsyrk("U", "T", p, n, 1D0, modx, n, 0D0, mxtmx, p)
   fmxtmx = mxtmx
   CALL dconvumtofull(p, fmxtmx)
   CALL dpotrf("U", p, fmxtmx, p, info)
   IF (info == 0) THEN
      CALL dpotri("U", p, fmxtmx, p, info)
      IF (info == 0) THEN
         vcovbeta = fmxtmx
      ELSE
         vcovbeta = info*1D0
      END IF
   ELSE
      vcovbeta = info*1D0
   END IF
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   dhuberpsi
!PART OF:      sctrob
!DESCRIPTION:  compute huber psi
!BENCHMARK:    robustbase:::huberPsi@psi (v 0.7-6),
!              approved June 18, 2011
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER n(1)
!  REAL k(1), vec(n)
!ON RETURN
!  REAL vec(n)
!--------------------------------------------------------------------
SUBROUTINE dhuberpsi(n, k, vec)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN) :: k
   DOUBLE PRECISION, INTENT(INOUT) :: vec(n)
   !local declarations
   INTEGER :: i
   DOUBLE PRECISION :: absvec(n)
   !
   absvec = ABS(vec)
   DO i = 1, n
      IF (absvec(i) >= k) THEN
         vec(i) = SIGN(k, vec(i))
      END IF
   END DO
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   drsaehubdest
!PART OF:      drsaehuber
!DESCRIPTION:  evaluates the (non-linear) estimating of d; it is called
!              from the modified zeroin function; kappa is the
!              consistency correction term
!BENCHMARK:
!DEPENDENCY:
!  dhuberpsi, dsqrtinvva
!ON ENTRY:
!  INTEGER n(1), g(1), nsize(g), dec(1), decorr(1)
!  REAL d(1), v(1), k(1), kappa(1), res(n)
!ON RETURN: REAL eval(1)
!--------------------------------------------------------------------
SUBROUTINE drsaehubdest(n, g, nsize, d, v, k, kappa, res, eval, dec, &
      decorr)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, g, dec, decorr
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: d, v, k, kappa
   DOUBLE PRECISION, INTENT(IN) :: res(n)
   DOUBLE PRECISION, INTENT(OUT) :: eval
   INTEGER :: i, j
   INTEGER :: l(g)
   DOUBLE PRECISION :: lhs, rhs, work
   DOUBLE PRECISION :: vec(n)
   !
   vec = res  !needed to create a new res-object because gfortran
   CALL dsqrtinvva(n, 1, g, nsize, d, v, 0, dec, decorr, vec)
   CALL dhuberpsi(n, k, vec)
   l(1) = 1
   DO i = 2, g
      l(i) = l(i - 1) + nsize(i - 1)
   END DO
   lhs = 0D0
   rhs = 0D0
   DO i = 1, g
      lhs = lhs + nsize(i) / (1 + d * nsize(i))
      work = 0D0
      DO j = 1, nsize(i)
         work = work + vec(l(i) + j - 1) * SQRT(1 / (1 + d * nsize(i)))
      END DO
      rhs = rhs + (work ** 2) / kappa
   END DO
   eval = lhs - rhs
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   drsaehubvest
!PART OF.      rsaehuber
!DESCRIPTION:  compute a Huber type M-scale^2 estimate; note that
!              stdres = v^(-1/2)*res, with res=y-x*beta and v is
!              parametrized in Hartley-Rao manner
!              notice, kappa denotes the consistency correction term
!              for Huber's Proposal 2 estimator.
!BENCHMARK:    MASS:::hubers(x, k, mu=0) (v 7.3-12),
!              approved on June 23, 2011
!DEPENDENCY:
!  dhuberwgt
!ON ENTRY:
!  INTEGER n(1), niter(1)
!  REAL k(1), acc(1), kappa(1), stdres(n)
!ON RETURN:
!  INTEGER info(1)
!  REAL v(1), sumwgt(1)
!--------------------------------------------------------------------
SUBROUTINE drsaehubvest(n, niter, v, k, acc, kappa, stdres, &
      sumwgt, info)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: niter
   INTEGER, INTENT(OUT) :: info
   DOUBLE PRECISION, INTENT(IN) :: k, acc, kappa
   DOUBLE PRECISION, INTENT(IN) :: stdres(n)
   DOUBLE PRECISION, INTENT(INOUT) :: v
   DOUBLE PRECISION, INTENT(OUT) :: sumwgt
   !local declarations
   INTEGER :: i, iter
   DOUBLE PRECISION :: ssq, vold
   DOUBLE PRECISION :: workresid(n)
   vold = v
   DO iter = 1, niter
      ssq = 0D0
      sumwgt = 0D0
      workresid = stdres / SQRT(vold)
      CALL dhuberwgt(n, k, 2, workresid)
      DO i = 1, n
         ssq = ssq + workresid(i) * stdres(i) ** 2
         sumwgt = sumwgt + workresid(i)
      END DO
      v = ssq / (n * kappa)
      IF (ABS(v/vold - 1D0) < acc) THEN
         EXIT
      ELSE
         vold = v
      END IF
   END DO
   info = iter
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   dhuberwgt
!PART OF:      sctrob
!DESCRIPTION:  compute huber psi-weight; NOTE:
!              typ = 1 for sqrt(wgt)
!              typ = 0 for wgt,
!              typ = 2 for wgt^2
!BENCHMARK:    robustbase:::huberPsi@wgt (v 0.7-6),
!              approved, June 19, 2011
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER n(1), typ(1)
!  REAL k(1), vec(n)
!ON RETURN:
!  REAL vec(n)
!--------------------------------------------------------------------
SUBROUTINE dhuberwgt(n, k, typ, vec)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, typ
   DOUBLE PRECISION, INTENT(IN) :: k
   DOUBLE PRECISION, INTENT(INOUT) :: vec(n)
   !local declarations
   INTEGER :: i
   DOUBLE PRECISION :: choice(n)
   !
   choice = k / ABS(vec)
   SELECT CASE (typ)
      CASE(1) !take the square root of the weights
         DO i = 1, n
            IF (choice(i) < 1D0) THEN
               vec(i) = SQRT(choice(i))
            ELSE
               vec(i) = 1
            END IF
         END DO
      CASE(0) !take the weights as they are
         DO i = 1, n
            IF (choice(i) < 1D0) THEN
               vec(i) = choice(i)
            ELSE
               vec(i) = 1
            END IF
         END DO
      CASE(2) !the weights to the power of two
         DO i = 1, n
            IF (choice(i) < 1D0) THEN
               vec(i) = choice(i) ** 2
            ELSE
               vec(i) = 1
            END IF
         END DO
      CASE DEFAULT !an errorneous call returns one
         vec = 0
   END SELECT
END SUBROUTINE
!
!====================================================================
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
!                 drsaehubdestiter, ddelta
!ON ENTRY:
!  INTEGER  n(1), p(1), g(1), niter(1), nsize(g), iter(2), dec(1),
!           decorr(1)
!  REAL     k(2), allacc(1), epsd(1), acc(3), sumwgt(3), xmat(n, p),
!           yvec(n), kappa(2), iterrecord(niter,p+2),
!ON RETURN:
!  INTEGER: niter(1), allconverged(1)
!  REAL:    tau(p+2), taurecord(niter, p+2)
!--------------------------------------------------------------------
SUBROUTINE drsaehub(n, p, g, niter, nsize, iter, iterrecord, allacc, &
      acc, sumwgt, xmat, yvec, k, kappa, epsd, tau, taurecord, &
      allconverged, dec, decorr)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec, decorr
   INTEGER, INTENT(IN) :: niter
   INTEGER, INTENT(IN) :: nsize(g)
   INTEGER, INTENT(IN) :: iter(2)
   INTEGER, INTENT(OUT):: allconverged
   DOUBLE PRECISION, INTENT(IN) :: allacc, epsd
   DOUBLE PRECISION, INTENT(IN) :: k(3)
   DOUBLE PRECISION, INTENT(IN) :: kappa(2)
   DOUBLE PRECISION, INTENT(IN) :: acc(3)
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(OUT) :: iterrecord(niter, 3)
   DOUBLE PRECISION, INTENT(INOUT) :: tau(p + 2)
   DOUBLE PRECISION, INTENT(OUT) :: taurecord(niter, p + 2)
   DOUBLE PRECISION, INTENT(OUT) :: sumwgt(3)
   INTEGER :: i, j, convergedbeta, work, monitord
   DOUBLE PRECISION :: upper
   DOUBLE PRECISION :: oldtau(p + 2)
   DOUBLE PRECISION :: res(n), stdres(n), sumwgtres(n)
   iterrecord = 0
   allconverged = 0
   monitord = 0
   DO i = 1, niter
      oldtau(1 : p) = tau(1 : p)
      oldtau(p + 1) = tau(p + 1)
      oldtau(p + 2) = tau(p + 2)
      CALL drsaebetaiter(n, p, g, k(1), xmat, yvec, tau(p+1), tau(p+2), &
         nsize, acc(1), tau(1:p), iter(1), convergedbeta, sumwgt(1), &
         work, dec, decorr)
      iterrecord(i, 1) = work
      IF (convergedbeta /= 1) THEN
         iterrecord(1, i) = (-1) * iterrecord(1, i)
      END IF
      res = yvec
      CALL dgemv("N", n, p, -1D0, xmat, n, tau(1:p), 1, 1D0, res, 1)
      stdres = res
      CALL dsqrtinvva(n, 1, g, nsize, tau(p+2), tau(p+1), 1, dec, &
         decorr, stdres)
      CALL drsaehubvest(n, iter(2), tau(p+1), k(2), acc(2), kappa(1), &
         stdres, sumwgt(2), work)
      iterrecord(i, 2) = work
      IF (monitord == 1) THEN
         tau(p+2) = 0D0
         iterrecord(i, 3) = 0D0
      ELSE
         upper = tau(p+2) * 1D1
         CALL drsaehubdestiter(n, g, nsize, tau(p+1), k(3), kappa(2), &
            res, 0D0, upper, acc(3), tau(p+2), work, dec, decorr)
         iterrecord(i, 3) = work
         IF (SUM(taurecord(MAX(i-2, 1):i, p+2)) < 3*epsd .AND. i >= 3) THEN
            monitord = 1
         END IF
      END IF
      taurecord(i, :) = tau
      CALL ddelta(p+1, oldtau, tau, allacc, allconverged)
      IF (allconverged == 1) THEN
         EXIT
      END IF
   END DO
   sumwgtres = res
   CALL dsqrtinvva(n, 1, g, nsize, tau(p+2), tau(p+1), 0, dec, &
      decorr, sumwgtres)
   CALL dhuberwgt(n, k(3), 0, sumwgtres)
   sumwgt(3) = 0D0
   DO j = 1, n
      sumwgt(3) = sumwgt(3) + sumwgtres(j)
   END DO
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   drsaeresid
!PART OF:      rsaehuber
!DESCRIPTION:  get the residuals and the huber wgt; i.e., res = e_ij
!              = y_ij - X_ij * beta - u_i; and stdres=V^(-1/2)*res;
!              the huber weight is w.r.t. to the model-psi (not the
!              prediction)
!BENCHMARK:
!DEPENDENCY:   dgemv(BLAS), dsqrtinvva, huberwgt
!ON ENTRY:
!  INTEGER n(1), p(1), g(1), nsize(g), dec(1)
!  REAL k(1), tau(p+2), xmat(n,p) yvec(n), u(g)
!ON RETURN:
!  REAL res(n), stdres(n), wgt(n)
!--------------------------------------------------------------------
SUBROUTINE drsaeresid(n, p, g, nsize, k, tau, u, xmat, yvec, res, &
      stdres, wgt, dec)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, g, dec
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: k !robustness tuning constant
   DOUBLE PRECISION, INTENT(IN) :: tau(p+2) ! (beta, v, d)
   DOUBLE PRECISION, INTENT(IN) :: u(g) !vector of random effects
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(OUT) :: res(n), stdres(n), wgt(n)
   !local declarations
   INTEGER :: i, j
   INTEGER :: l(g)
   res = yvec
   CALL dgemv("N", n, p, -1D0, xmat, n, tau(1:p), 1, 1D0, res, 1)
   l(1) = 1
   DO i = 2, g
      l(i) = l(i - 1) + nsize(i - 1)
   END DO
   !
   DO i = 1, g
      DO j = 1, nsize(i)
         res(l(i) + j - 1) = res(l(i) + j - 1) - u(i)
      END DO
   END DO
   stdres = res
   CALL dsqrtinvva(n, 1, g, nsize, tau(p+2), tau(p+1), 0, dec, 0, stdres)
   wgt = stdres
   CALL dhuberwgt(n, k, 0, wgt)
END SUBROUTINE
!
!====================================================================
!SUBROUTINE:   dmedmed
!AUTHOR:       Tobias Schoch, November 19, 2011
!PART OF:      sctrob
!DESCRIPTION:  computes the median (type=2 in R; i.e., average over
!              discontinuities) and compute the mad
!              typ = 0; median
!              typ = 1; mad
!DEPENDENCY:   it uses qsort3 from R's utility functions (for C, it is
!              in the header file R_ext/Utils.h (included in R.h)
!ON ENTRY:
!  INTEGER n(1), typ(1)
!  REAL vx(n)
!ON RETURN
!  REAL res(1)
!--------------------------------------------------------------------
SUBROUTINE dmedmad(n, vx, typ, res)
   IMPLICIT NONE
   INTEGER, INTENT(IN):: n, typ
   DOUBLE PRECISION, INTENT(IN) :: vx(n)
   DOUBLE PRECISION, INTENT(OUT) :: res
   INTEGER :: i
   DOUBLE PRECISION :: work(n)
   work = vx
   CALL qsort3(work, 1, n)
   IF (MODULO(n, 2) == 0) THEN
      res = ( work(n / 2) + work((n / 2) + 1) ) / 2D0
   ELSE
      res = work( ( (n - 1) / 2) + 1)
   END IF
   IF (typ == 1) THEN
      DO i = 1, n
         work(i) = ABS(work(i) - res)
      END DO
      CALL qsort3(work, 1, n)
      IF (MODULO(n, 2) == 0) THEN
         res = ( 1.4814 / 2D0 ) * ( work(n / 2) + work((n / 2) + 1) )
      ELSE
         res = 1.4814 * work( ( (n - 1) / 2) + 1)
      END IF
   END IF
END SUBROUTINE
!====================================================================
!SUBROUTINE:   drlm
!PART OF:      sctrob
!AUTHOR:       Tobias Schoch, November 20, 2011
!DESCRIPTION:  compute regression M-estimates, using the mad as the
!              preliminary scale estimate
!DEPENDENCY:
!  dgemv (BLAS2 and LAPACK), dgels (LAPACK)
!  dhuberwgt, dmedmad
!ON ENTRY:
!  INTEGER n(1), p(1),
!  REAL xmat(n, p), yvec(n), k(1), beta(p), acc(1)
!ON RETURN:
!  INTEGER info(1)
!  REAL beta(p), s(1)
!--------------------------------------------------------------------
SUBROUTINE drlm(n, p, xmat, yvec, k, beta, s, info, niter, acc)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, p, niter
   INTEGER, INTENT(OUT) :: info
   DOUBLE PRECISION, INTENT(IN) :: k, acc
   DOUBLE PRECISION, INTENT(IN) :: yvec(n)
   DOUBLE PRECISION, INTENT(IN) :: xmat(n, p)
   DOUBLE PRECISION, INTENT(OUT) :: s
   DOUBLE PRECISION, INTENT(INOUT) :: beta(p)
   !local declarations (most are used in dqrls)
   INTEGER :: i, j, l, lwork, converged
   INTEGER, PARAMETER :: lworkmax = 10000
   DOUBLE PRECISION :: oldbeta(p)
   DOUBLE PRECISION :: work(lworkmax)
   DOUBLE PRECISION :: modyvec(n), res(n)
   DOUBLE PRECISION :: modxmat(n, p)
   DO l = 1, niter
      oldbeta = beta
      res = yvec
      CALL dgemv("N", n, p, -1D0, xmat, n, oldbeta, 1, 1D0, res, 1)
      CALL dmedmad(n, res, 1, s)
      DO i = 1, n
         res(i) = res(i) / s
      END DO
      CALL dhuberwgt(n, k, 1, res)
      DO j = 1, p
         DO i = 1, n
            modyvec(i) = yvec(i) * res(i)
            modxmat(i, j) = xmat(i, j) * res(i)
         END DO
      END DO
      lwork = -1
      CALL dgels("n", n, p, 1, modxmat, n, modyvec, n, work, lwork, info)
      lwork = MIN(lworkmax, INT(work(1)))
      CALL dgels("n", n, p, 1, modxmat, n, modyvec, n, work, lwork, info)
      IF (info == 0) THEN
         beta = modyvec(1:p)
      ELSE
         beta = 0D0
      END IF
      CALL ddelta(p, oldbeta, beta, acc, converged)
      IF (converged == 1) THEN
         info = l
         EXIT
      END IF
   END DO
END SUBROUTINE
