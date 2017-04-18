program main_prog_hydrogen
use equation_module

implicit none
real(8) :: ri, rf, dr, dE, dEMax, energyMax, energyStepCoeff
real(8), allocatable :: wave_func_vec(:,:), coordinates(:)
integer :: N, Nr
integer :: i,k,j
integer :: controller, noChanges
real(8) :: closeness, closeness_mem
real(8) :: func_memL, func_memR
integer :: roots
real(8) :: u1_0(2), u2_0(2)
real(8) :: u1(2), u2(2)
character::rubbish

print*, 'Enter the nuclear charge Z = '
read*, z

print*, 'Enter the orbital quantum number l = '
read*, l

print*, 'Enter the main quantum number '
read *, Nr

if (Nr <= l) then
    print *, "Invalid Nr! It must be greater than l."
end if

energyMax = -real(Z)**2*1.17
energy = energyMax

if (energy .ge. 0.0D+00) then
    print*, 'the energy needs to be negative, type another number '
    read*, energy
endif

100 format(3(A15))
150 format(2(A15,f17.7))
200 format(6(e17.7))

open(unit=18, file='hydrogen_wavefunc.txt', status='replace')
!hydrogen_wavefunc.txt - file with radius, function and dericative values
open(unit=19, file='energy.txt', status='replace')
!energy.txt - file with energy and closeness values

write(19,100) 'energy', 'closeness'

!initialize cicle counter
controller=0

dEMax = abs(energy/500)
dE = dEMax
!initialize a big parameter to provoke the next step in the cicle
closeness_mem = 10000.0D+00

do
!do cicle - is looking for necessary energy value
!also is calculating wave function
!control the number of cicle steps
controller=controller+1
  ri = 1.0e-8
  !left start point
  rf = log(1.0D-04)/(-sqrt(-2*energy))

  do
    if ((rf**(l/2) * exp(-sqrt(-2*energy)*rf)) .le. 1.0D-04) then
    !if ((rf**l * exp(-sqrt(-2*energy) * rf))<= 1.0e-3) then
        exit
    else
        rf = rf*1.2D+00
        print *, rf
    end if
  end do

  !right start point (finish for the wave function)

u1_0 = (/ 0.0D+00, 1.0D+00 /)
!u1_0 - initial values for f and f' if R = r^l * f(r)
!u1 = (f, f')
u2_0 = (/ 1.0D+00, 0.0D+00 /)
!u2_0 - initial values for g and g' if R = exp( - sqrt( -2*energy ) ) * g(r)
!u2 - (g, g')

N = int( ( rf - ri )/5.0e-2 )*2
!2*N - number of small intervals on [ri,rf]

dr = ( rf - ri )/( 2*N )
!step size
!print*, 'dr = ', dr

allocate( wave_func_vec(2*N+1,2) )
!wave_func_vec - matrix, where we save wave function and its derivative values
allocate( coordinates(2*N+1) )
!coordunates - vector, where we save the radius points

wave_func_vec(1, 1) = u1_0(1) * ri**l
wave_func_vec(1, 2) = u1_0(2) * ri**l + l* ri**(l-1) * u1_0(1)

func_memL = 0.0D+00
func_memR = 0.0D+00

roots = 0

wave_func_vec(2*N+1, 1) = u2_0(1) * exp( - sqrt( -2 * energy)*rf )
wave_func_vec(2*N+1, 2) = u2_0(2) * exp( - sqrt( -2 * energy)*rf ) + &
( - sqrt( -2 * energy) ) *exp(- sqrt( -2 * energy)*rf ) *u2_0(1)

coordinates(1)=ri
coordinates(2*N+1)=rf

do k=1,N

call rk4vec( ri, 2, u1_0, dr, eq_modern1, u1 )
!runge-kutta for left part of wave function
call rk4vec( rf, 2, u2_0, -dr, eq_modern2, u2 )
!runge-kutta for right part of wave function

ri=ri+dr
rf=rf-dr

u1_0 = u1
u2_0 = u2

if (func_memL * u1(1) < 0) then
    roots = roots + 1
end if
if (func_memR * u2(1) < 0) then
    roots = roots + 1
end if


func_memL = u1(1)
func_memR = u2(1)

u1(2) = u1(2) * ri**l + l* ri**(l-1) * u1(1)
u1(1) = u1(1) * ri**l

wave_func_vec(k+1, 1:2) =   u1

u2(2) = u2(2) * exp(- sqrt( -2 * energy)*rf ) + ( - sqrt( -2 * energy) ) *exp(- sqrt( -2 * energy)*rf ) * u2(1)
u2(1) = u2(1) * exp(- sqrt( -2 * energy)*rf )

wave_func_vec(2*N-k+1, 1:2) = u2

coordinates(k+1)=ri
coordinates(2*N-k+1)=rf

end do
!print*, ri, rf
if ( abs( ri - rf ) .ge. 1.0e-3*dr ) then
    print*, ri, rf
    print*, 'Finish is not in central point!!!'
endif

wave_func_vec(N+1:2*N+1, 1:2) = wave_func_vec(N+1:2*N+1, 1:2) /u2(1)
wave_func_vec(1:N, 1:2) = wave_func_vec(1:N, 1:2) / u1(1)
!we need to scale one of the functions (right or left)
!because initial values are correct except for random constant factor
u2 = u2 / u2(1)
u1 = u1 / u1(1)
!measure of how close the functions are
!Has to be relative, since the more the n (main quantum number) of the radial function is,
!the smaller it is.

!closeness = ( u1(2)- u2(2) ) / max( abs(u1(2)), abs(u2(2)) )
!closeness = (u1(2) - u2(2)) / max( abs(u1(2)), abs(u2(2)) )
closeness = abs(u1(2) - u2(2))/max(abs(u1(2)),abs(u2(2)))
if (abs(closeness_mem) > 0.9D+04) then
    closeness_mem = closeness
    print *, "null closeness"
end if
!derivatives equality condition
energyStepCoeff = (energy/energyMax) ** (0.1)
if ( abs(closeness) .le. 1.0e-2 ) then
    print *, "Found state l=",l," N=", sqrt(-Z**2/(2*energy))," roots=",roots," energy steps=",controller
    if (roots == (Nr - l - 1)) then
        exit
    else
        print *,"Not the needed state, print something to continue or END to finish and write to file"
        read *, rubbish
        if (rubbish == "E") then
            exit
        end if
        dE = sign(dEMax, real(-roots + (Nr - l - 1),8))
        energy = energy + 50*dE
        closeness_mem = 1.0D+04
        controller = 0
    end if
elseif ((((closeness - closeness_mem)*dEMax/(abs(dE * energyStepCoeff))) .ge. 0.2D+00 ).and.( abs(closeness_mem) < 0.9D+03 )) then
    dE = - dE/2
    print * , "reverse, E= ", energy, " now dE=", dE, " closeness=", closeness, "n=", sqrt(-Z**2/(2*energy))
    noChanges = 0
    !exit
end if

if (noChanges > 9) then
    if (abs(dE)< abs(dEMax)) then
        dE = dE * 2
        noChanges = 0
    end if
end if

write(19,200) energy, closeness, u1(2), u2(2), real(roots,8)

energy = energy + dE*energyStepCoeff
if (energy >= 0) then
    print *, "Bad energy"
    exit
  end if
!if ((closeness - closeness_mem)/max(abs(closeness),abs(closeness_mem)) > 0.1D+00) then
!    exit
!end if
closeness_mem = closeness

!put the condition to avoid infinite looping
if (controller .ge. 3e3) then
    print*, 'Can not find the necessary energy point for l = ', l
    exit
endif

deallocate( wave_func_vec )
deallocate( coordinates )

end do !cicle for energy and wave function

write(18,*) '*********************************************************************'
write(18,150) 'E=', energy, 'theoretical n=', sqrt(-z**2/(2*energy))
write(18,*) '*********************************************************************'
write(18,100) 'r', 'wave func(r)', 'derivative'

do k=1,2*N+1
write(18,200) coordinates(k), wave_func_vec(k,1), wave_func_vec(k,2)
end do

print*, 'Number of energy steps = ', controller

deallocate( wave_func_vec )
deallocate( coordinates )

close(18)
close(19)

end program main_prog_hydrogen
