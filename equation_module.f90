module equation_module

integer, public:: z,l
real(8), public:: energy
real(8), dimension(2), parameter:: uno=(/ 1.0D+00, 1.0D+00 /)

contains
!=======================================================================
!initially wanted to solve the equation on radial function R(r)
!where full wave function = R(r)/r * Y_lm
subroutine equations ( r, n, u, uprime )
! this subroutine evaluates the derivative uprime(n) given the time r and
!    solution vector u(n)
! r - independent variable
! n - number of equations in the system
! u(n) - solution vector (in our case u = (R, R') )
! uprime(n) - the derivative of u(n)

integer:: n
real(8):: r, u(n), uprime(n)

uprime(1) = u(2)
uprime(2) = - ( real(2*z)/r + 2*energy - real(l*(l+1))/r**2 ) * u(1)

end subroutine

!=============================================================
!then decided to do replacement
!=============================================================
!for small r we know the asymptotics R <--> r^(l+1)
!introduce new function f: R = r^l * f
!obtain new equation
subroutine eq_modern1( r, n, u, uprime )
! this subroutine evaluates the derivative uprime(n) given the time r and
!    solution vector u(n)
! r - independent variable
! n - number of equations in the system
! u(n) - solution vector (in our case u = (f, f') )
! uprime(n) - the derivative of u(n)

integer:: n
real(8):: r, u(n), uprime(n)

uprime(1) = u(2)
uprime(2) = -real(2*l)/r*u(2) - ( 2*energy - real(2*l)/r**2 + real(2*z)/r )*u(1)

end subroutine

!==============================================================
!for big r we know the asymptotics R <--> exp(-sqrt(-2*energy)*r)
!introduce new function f: R = exp(-sqrt(-2*energy)*r) * g
!obtain new equation
subroutine eq_modern2( r, n, u, uprime )
! this subroutine evaluates the derivative uprime(n) given the time r and
!    solution vector u(n)
! r - independent variable
! n - number of equations in the system
! u(n) - solution vector (in our case u = (g, g') )
! uprime(n) - the derivative of u(n)

integer:: n
real(8):: r, u(n), uprime(n)

uprime(1) = u(2)
uprime(2) = 2*sqrt(-2*energy)*u(2) - ( real(2*z)/r - real(l*(l+1))/r**2 )*u(1)

end subroutine

!===============================================================

end module equation_module
