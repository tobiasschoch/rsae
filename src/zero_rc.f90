!   ZERO_RC seeks the root of a function F(X) using reverse communication.
!   Copyright (C) 2013 John Burkardt
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Lesser General Public
!   License as published by the Free Software Foundation; either
!   version 2.1 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!   Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public
!   License along with this library; if not, write to the Free Software
!   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
!
!*******************************************************************************
!
!   I, the author of the R-package rsae (Tobias Schoch), obtained this
!   subroutine in the file Brent.f90 from John Burkardt's online repository,
!   https://people.math.sc.edu/Burkardt/f_src/brent/brent.html, on
!   January 23, 2024.
!
subroutine zero_rc ( a, b, t, arg, status, value )
!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent zero finder
!    algorithm, using reverse communication.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, double precission A, B, the endpoints of the change of
!    sign interval.
!
!    Input, double precission T, a positive error tolerance.
!
!    Output, double precission ARG, the currently considered point.  The user
!    does not need to initialize this value.  On return with STATUS positive,
!    the user is requested to evaluate the function at ARG, and return
!    the value in VALUE.  On return with STATUS zero, ARG is the routine's
!    estimate for the function's zero.
!
!    Input/output, integer ( kind = 4 ) STATUS, used to communicate between
!    the user and the routine.  The user only sets STATUS to zero on the first
!    call, to indicate that this is a startup call.  The routine returns STATUS
!    positive to request that the function be evaluated at ARG, or returns
!    STATUS as 0, to indicate that the iteration is complete and that
!    ARG is the estimated zero
!
!    Input, double precission VALUE, the function value at ARG, as requested
!    by the routine on the previous call.
!
  implicit none

  double precision a
  double precision arg
  double precision b
  double precision, save :: c
  double precision, save :: d
  double precision, save :: e
  double precision, save :: fa
  double precision, save :: fb
  double precision, save :: fc
  double precision m
  double precision, save :: machep
  double precision p
  double precision q
  double precision r
  double precision s
  double precision, save :: sa
  double precision, save :: sb
  integer status
  double precision t
  double precision tol
  double precision value
!
!  Input STATUS = 0.
!  Initialize, request F(A).
!
  if ( status == 0 ) then

    machep = epsilon ( a )

    sa = a
    sb = b
    e = sb - sa
    d = e

    status = 1
    arg = a
    return
!
!  Input STATUS = 1.
!  Receive F(A), request F(B).
!
  else if ( status == 1 ) then

    fa = value

    status = 2
    arg = sb
    return
!
!  Input STATUS = 2
!  Receive F(B).
!
  else if ( status == 2 ) then

    fb = value

    if ( 0.0D+00 < fa * fb ) then
      status = -1
      return
    end if

    c = sa
    fc = fa

  else

    fb = value

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end if
!
!  Compute the next point at which a function value is requested.
!
  if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

  end if

  tol = 2.0D+00 * machep * abs ( sb ) + t
  m = 0.5D+00 * ( c - sb )

  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
    status = 0
    arg = sb
    return
  end if

  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

  else

    s = fb / fa

    if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

    else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

    end if

    if ( 0.0D+00 < p ) then
      q = - q
    else
      p = - p
    end if

    s = e
    e = d

    if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
    else
      e = m
      d = e
    end if

  end if

  sa = sb
  fa = fb

  if ( tol < abs ( d ) ) then
    sb = sb + d
  else if ( 0.0D+00 < m ) then
    sb = sb + tol
  else
    sb = sb - tol
  end if

  arg = sb
  status = status + 1

  return
end
