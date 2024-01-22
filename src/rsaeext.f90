!####################################################################
!        The subroutine "drsaehubdestiter" derives code from Richard
!        Brent's Fortran77 function "zeroin"; cf., Brent, R.P.
!        1973, Algorithms for Minimization without Derivatives
!        (Englewood Cliffs, NJ: Prentice- Hall).
!
!        I obtained the code of "zeroin.f" from
!        http://www.netlib.org/go/zeroin.f,
!        on June 24, 2011.
!
!        Richard's original code is licensed under GNU General
!        Public License (cf.,
!        http//:maths.anu.edu.au/~brent/software.html; June 24, 2011).
!
!        I modified the original code in order that it meets
!        the gfortran f90 standards, and added some specific
!        blocks that I need for my computations.
!
!        I therefore license the subroutine "drsaehubdestiter"
!        also under GNU General Public License; GPL >=2.
!
!        You can find a copy of the GPL-2 license under
!        $RHOME/share/licenses/GPL-2
!        where $RHOME denotes the root directory of your R
!        (r-project.org) installation.
!
!####################################################################
!
!====================================================================
!SUBROUTINE:   drsaehubdestiter
!PART OF:      rsae
!DESCRIPTION:  root-finding device; evaluates drsaehubdest between
!              lower and upper; there must be a sign change in the
!              interval, otherwise the sub stops
!DEPENDENCY:   none
!ON ENTRY:
!  INTEGER n(1), g(1), nsize(g), dec(1)
!  REAL v(1), k(1), kappa(1), lower(1), upper(1), tol(1), res(n)
!ON RETURN
!  INTEGER info(1)
!  REAL zeroin(1)
!--------------------------------------------------------------------
!
SUBROUTINE drsaehubdestiter(n, g, nsize, v, k, kappa, res, lower, upper, &
      tol, zeroin, info, dec, decorr)
   IMPLICIT NONE
   !zeroin declarations
   INTEGER, INTENT(OUT) :: info  !info was not in Richard's original
   !defintion (since zeroin was a function)
   DOUBLE PRECISION, INTENT(IN) :: lower, upper, tol
   DOUBLE PRECISION, INTENT(OUT) :: zeroin
   !foo declarations
   INTEGER, INTENT(IN) :: n, g, dec, decorr
   INTEGER, INTENT(IN) :: nsize(g)
   DOUBLE PRECISION, INTENT(IN) :: v, k, kappa
   DOUBLE PRECISION, INTENT(IN) :: res(n)
   !local declarations
   INTEGER, PARAMETER :: ITMAX = 100
   INTEGER :: iter
   DOUBLE PRECISION :: a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s
   !define a precision constant, cf R .Machine$double.eps for double
   !data type, i.e., 64bit sized blocks
   DOUBLE PRECISION, PARAMETER :: EPS = 2.3E-16
   !initialize
   info = 0
   !start with checking whether there is a sign change, if not break
   a = lower
   b = upper
   !note that these function calls have been added to Richard's original code
   !compute fa (i.e., evaluation for a = lower)
   CALL drsaehubdest(n, g, nsize, a, v, k, kappa, res, fa, dec, decorr)
   !compute fb (i.e., evaluation for b = upper)
   CALL drsaehubdest(n, g, nsize, b, v, k, kappa, res, fb, dec, decorr)
   !check that f(lower) and f(upper) have different signs
   IF ((fa > 0D0 .and. fb > 0D0) .or. (fa < 0D0 .and. fb < 0D0)) THEN
      !return -1; this violates the assumption and therefore be detected
      info = -1
      zeroin = 0D0
   ELSE
      c = b
      fc = fb
      DO iter = 1,ITMAX
         IF ((fb > 0D0 .and. fc > 0D0) .or. &
               (fb < 0D0 .and. fc < 0D0)) THEN
            c = a
            fc = fa
            d = b - a
            e = d
         END IF
         IF (ABS(fc) < ABS(fb)) THEN
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         END IF
         tol1 = 2D0 * EPS * ABS(b)+ 5D-1 * tol
         xm = 5D-1 * (c-b)
         IF (ABS(xm) <= tol1 .or. fb == 0D0) THEN
            zeroin = b
            !if it converges, then info=# iterations
            info = iter
            EXIT
         END IF
         ! see if a bisection is forced
         IF (ABS(e) >= tol1 .and. ABS(fa) > ABS(fb)) THEN
            s = fb / fa
            IF (a == c) THEN
               !linear interpolation
               p = 2D0 * xm * s
               q = 1D0 - s
            ELSE
               !inverse quadratic interpolation
               q = fa / fc
               r = fb / fc
               p = s * (2D0 * xm * q * (q - r) - (b - a) * (r - 1D0))
               q = (q - 1D0) * (r - 1D0) * (s - 1D0)
            END IF
            IF (p > 0D0) THEN
               q = -q
            END IF
            p = ABS(p)
            IF (2D0 * p  <  MIN(3D0 * xm * q - ABS(tol1 * q), &
                  ABS(e * q))) THEN
               e = d
               d = p / q
            ELSE
               d = xm
               e = d
            END IF
         ELSE
            d = xm
            e = d
         END IF
         a = b
         fa = fb
         !note merge and sign work that way in gfortran (which may no
         !be true
         !for others; I did not check it)
         b = b + MERGE(d, SIGN(tol1, xm), ABS(d) > tol1)
         ! note this function call has been added to Richard's origial
         !code
         !eval function at b
         CALL drsaehubdest(n, g, nsize, b, v, k, kappa, res, fb, dec,&
            decorr)
      END DO
      !return b; however, it did not converge, therefore info=0, still
      zeroin = b
   END IF
END SUBROUTINE
!
!####################################################################
! EOF rsaeext.f90
!####################################################################
