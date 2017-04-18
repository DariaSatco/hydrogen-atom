program main_prog_hydrogen
use equation_module

implicit none
real(8) :: ri, rf,rfl, dr, error, dE, dEMax, energyMax, energyStepCoeff
real(8), allocatable :: wave_func_vec(:,:), coordinates(:)
integer :: N1, N2, Nr
integer :: i,k,j
integer :: controller, noChanges, rootsNeeded
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

error = 1.0D-04

if (Nr <= l) then
    print *, "Invalid Nr! It must be greater than l."
end if

energyMax = -real(Z)**2*0.57
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

dEMax = abs(energy/100)
dE = dEMax

rootsNeeded = (Nr - l - 1)

!initialize a big parameter to provoke the next step in the cicle
closeness_mem = 10000.0D+00

do
!do cicle - is looking for necessary energy value
!also is calculating wave function
!control the number of cicle steps
controller=controller+1
  ri = 1.0e-8
  !left start point
  rfl = abs(energy/Z)
  rf = abs(log(1.0e-8))/(sqrt(-2*energy))

  rf = rfl + rf

  !right start point (finish for the wave function)

u1_0 = (/ 0.0D+00, 1.0D+00 /)
!u1_0 - initial values for f and f' if R = r^l * f(r)
!u1 = (f, f')
u2_0 = (/ 1.0D+00, 0.0D+00 /)
!u2_0 - initial values for g and g' if R = exp( - sqrt( -2*energy ) ) * g(r)
!u2 - (g, g')

dr=1.0e-3
N1 = int( ( rfl - ri )/dr )
N2 = int( (rf - rfl)/dr )
!2*N - number of small intervals on [ri,rf]

dr = ( rf - ri )/( (N1+N2) )
!step size
!print*, 'dr = ', dr

allocate( wave_func_vec(N1+N2+1,2) )
!wave_func_vec - matrix, where we save wave function and its derivative values
allocate( coordinates(N1+N2+1) )
!coordunates - vector, where we save the radius points

wave_func_vec(1, 1) = u1_0(1) * ri**l
wave_func_vec(1, 2) = u1_0(2) * ri**l + l* ri**(l-1) * u1_0(1)

func_memL = 0.0D+00
func_memR = 0.0D+00

roots = 0

wave_func_vec(N1+N2+1, 1) = u2_0(1) * exp( - sqrt( -2 * energy)*rf )
wave_func_vec(N1+N2+1, 2) = u2_0(2) * exp( - sqrt( -2 * energy)*rf ) + &
( - sqrt( -2 * energy) ) *exp(- sqrt( -2 * energy)*rf ) *u2_0(1)

coordinates(1)=ri
coordinates(N1+N2+1)=rf

do k=1,N1
    call rk4vec( ri, 2, u1_0, dr, eq_modern1, u1 )
    !runge-kutta for left part of wave function
    ri=ri+dr
    u1_0 = u1
    if (func_memL * u1(1) < 0) then
        roots = roots + 1
    end if
    func_memL = u1(1)

    u1(2) = u1(2) * ri**l + l* ri**(l-1) * u1(1)
    u1(1) = u1(1) * ri**l

    wave_func_vec(k+1, 1:2) =   u1
    coordinates(k+1)=ri
end do


do k=1,N2
    call rk4vec( rf, 2, u2_0, -dr, eq_modern2, u2 )
    !runge-kutta for right part of wave function
    rf=rf-dr
    u2_0 = u2

    if (func_memR * u2(1) < 0) then
        roots = roots + 1
    end if
    func_memR = u2(1)

    u2(2) = u2(2) * exp(- sqrt( -2 * energy)*rf ) + ( - sqrt( -2 * energy) ) *exp(- sqrt( -2 * energy)*rf ) * u2(1)
    u2(1) = u2(1) * exp(- sqrt( -2 * energy)*rf )

    wave_func_vec(N1+N2-k+1, 1:2) = u2
    coordinates(N1+N2-k+1)=rf
end do

!print*, ri, rf
if ( abs( ri - rf ) .ge. 1.0e-3*dr ) then
    print*, ri, rf
    print*, 'Finish is not in central point!!!'
endif

print *, "energy=",energy," step=",dE*energyStepCoeff
!print *, "rf=",rf," energy/Z=",abs(energy/Z), "diff=",rf-abs(energy/Z)

wave_func_vec(N1+1:N1+N2+1, 1:2) = wave_func_vec(N1+1:N1+N2+1, 1:2) /u2(1)
wave_func_vec(1:N1, 1:2) = wave_func_vec(1:N1, 1:2) / u1(1)
!we need to scale one of the functions (right or left)
!because initial values are correct except for random constant factor
u2 = u2 / u2(1)
u1 = u1 / u1(1)
!measure of how close the functions are
!Has to be relative, since the more the n (main quantum number) of the radial function is,
!the smaller it is.
closeness = u1(2) - u2(2)

energyStepCoeff = abs(energy/energyMax)

if (abs(dE*energyStepCoeff/energy) < error) then

    print*, "Found state L=",l," roots=",roots, " energy=", energy," n=", sqrt(-Z**2/(2*energy))

    if ( roots == rootsNeeded ) then
        exit
    else
        dE = -sign(dEMax, real(roots - rootsNeeded,8))
        energy = energy + 50.0D+00 * dE * energyStepCoeff
        closeness = 10000.0D+00
        controller = 0
    end if

elseif ( (closeness * closeness_mem) < 0.0D+00) then
    dE = -dE/2
    print * , "reverse, E= ", energy, " now dE=", dE, " closeness=", closeness, "n=", sqrt(-Z**2/(2*energy))

end if


write(19,200) energy, closeness, u1(2), u2(2), real(roots,8)

energy = energy + dE * energyStepCoeff
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

do k=1,N1+N2+1
write(18,200) coordinates(k), wave_func_vec(k,1), wave_func_vec(k,2)
end do

print*, 'Number of energy steps = ', controller

deallocate( wave_func_vec )
deallocate( coordinates )

close(18)
close(19)

end program main_prog_hydrogen
