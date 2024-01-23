!###############################################################################
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
!###############################################################################
!
!===============================================================================
!SUBROUTINE:   drsaehubdestiter
!PART OF:      rsae
!DESCRIPTION:  root-finding device; evaluates drsaehubdest between
!              lower and upper; there must be a sign change in the
!              interval, otherwise the sub stops
!DEPENDENCY:   none
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
    integer, intent(out) :: info  !info was not in richard's original
    !defintion (since zeroin was a function)
    double precision, intent(in) :: lower, upper, tol
    double precision, intent(out) :: zeroin
    !foo declarations
    integer, intent(in) :: n, g, dec, decorr
    integer, intent(in) :: nsize(g)
    double precision, intent(in) :: v, k, kappa
    double precision, intent(in) :: res(n)
    !local declarations
    integer, parameter :: ITMAX = 100
    integer :: iter
    double precision :: a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s
    !define a precision constant, cf R .Machine$double.eps for double
    !data type, i.e., 64bit sized blocks
    double precision, parameter :: EPS = 2.3e-16
    !initialize
    info = 0
    !start with checking whether there is a sign change, if not break
    a = lower
    b = upper
    !note that these function calls have been added to Richard's original code
    !compute fa (i.e., evaluation for a = lower)
    call drsaehubdest(n, g, nsize, a, v, k, kappa, res, fa, dec, decorr)
    !compute fb (i.e., evaluation for b = upper)
    call drsaehubdest(n, g, nsize, b, v, k, kappa, res, fb, dec, decorr)
    !check that f(lower) and f(upper) have different signs
    if ((fa > 0d0 .and. fb > 0d0) .or. (fa < 0d0 .and. fb < 0d0)) then
        !return -1; this violates the assumption and therefore be detected
        info = -1
        zeroin = 0d0
    else
        c = b
        fc = fb
        do iter = 1, ITMAX
            if ((fb > 0d0 .and. fc > 0d0) .or. (fb < 0d0 .and. fc < 0d0)) &
                    then
                c = a
                fc = fa
                d = b - a
                e = d
            end if
            if (abs(fc) < abs(fb)) then
                a = b
                b = c
                c = a
                fa = fb
                fb = fc
                fc = fa
            end if
            tol1 = 2d0 * EPS * abs(b) + 5d-1 * tol
            xm = 5d-1 * (c - b)
            if (abs(xm) <= tol1 .or. fb == 0d0) then
                zeroin = b
                !if it converges, then info=# iterations
                info = iter
                exit
            end if
            ! see if a bisection is forced
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                s = fb / fa
                if (a == c) then
                    !linear interpolation
                    p = 2d0 * xm * s
                    q = 1d0 - s
                else
                    !inverse quadratic interpolation
                    q = fa / fc
                    r = fb / fc
                    p = s * (2d0 * xm * q * (q - r) - (b - a) * (r - 1d0))
                    q = (q - 1d0) * (r - 1d0) * (s - 1d0)
                end if
                if (p > 0d0) then
                    q = -q
                end if
                p = abs(p)
                if (2d0 * p  <  min(3d0 * xm * q - abs(tol1 * q), &
                        abs(e * q))) then
                    e = d
                    d = p / q
                else
                    d = xm
                    e = d
                end if
            else
                d = xm
                e = d
            end if
            a = b
            fa = fb
            !note merge and sign work that way in gfortran
            b = b + merge(d, sign(tol1, xm), abs(d) > tol1)
            !note this function call has been added to Richard's origial
            !code eval function at b
            call drsaehubdest(n, g, nsize, b, v, k, kappa, res, fb, dec, &
                decorr)
        end do
        !return b; however, it did not converge, therefore info=0, still
        zeroin = b
    end if
end subroutine
!###############################################################################
! EOF rsaeext.f90
!###############################################################################
