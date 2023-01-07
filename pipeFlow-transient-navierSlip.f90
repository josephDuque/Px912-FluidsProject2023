program main
use iso_fortran_env, only : real64, int64
implicit none
!Number of nodes
integer, parameter :: n = 251, maxIter = 1000000, nPlots=20
integer :: info , i, j, totalTimeSteps, timeStep4Plot, t, col
integer, dimension ( n ) :: ipiv
! Constants G, nu and boundary length Wi
real(real64), parameter :: mu = (9.6_real64)*1.0e-4, rho = 730.0_real64, & ! viscosity and density of decane
  W = (22.0_real64)*1.0e-9, l=(41.0_real64)*1.0e-9, &
  totalSimTime = (1.0_real64)*1.0e-8, timeStep = (2.0_real64)*1.0e-15, pi=4.0_real64*ATAN(1.0_real64) 
real(real64) :: G, nu, eta, volumeFlux_ana, volumeFlux_num
! An array of real numbers of size nxn
real(real64), dimension (n , n ) :: A
! Note solution is now contained in f and r is spatial coordinates
real(real64), dimension ( n ) :: r, f, temp, source, f_ana, integrand, absErr, relErr
real(real64), allocatable, dimension(:,:) :: data2Plot
real(real64) :: h ! Mesh step 
character(len=10) :: fmt

!---------------------------------------------------------
! Initialsing kinematic viscosity and the source G, 
! calculated from the density and the pressure gradient, 
! estimated as the pressure difference/pipe-length 
! (length scales small enough to approximate as linear).
nu = mu/rho
G = ((7.1_real64)*1.0e3)/((78.0_real64)*1.0e-6) * 1.0_real64/rho
!---------------------------------------------------------

!---------------------------------------------------------
! Distance between nodes , note we make ’n ’ the integer into a real number using ’ real ’
h = W/(real(n, real64) - 1.0_real64)
!timeStep = (0.1*h**2)/nu
totalTimeSteps = int(totalSimTime/timeStep+1, int64)
!---------------------------------------------------------
allocate(data2Plot(n,nPlots))
timeStep4Plot = totalTimeSteps/nPlots
!---------------------------------------------------------
! Identify nodal positions x(i) and create analytic solution f_ana
do i =1 , n
  r(i) = (real(i, real64) - 1.0_real64) * h
  f_ana(i) = (G/(4.0_real64*nu))*(W**2 - r(i)**2 + 2.0_real64*W*l)
end do
!---------------------------------------------------------

!---------------------------------------------------------
! Initialise all values of A, f and source to be zero
A = 0.0_real64 ; f = 0.0_real64 ; source = 0.0_real64
!---------------------------------------------------------

!---------------------------------------------------------
! Fill the matrices A and b with known values by looping through the bulk / non-boundary nodes
do i = 2 ,n-1
  A(i, i-1) = (nu*timeStep)/(h**2) - (timeStep*nu)/(2.0_real64*h*r(i))
  A(i, i) = 1 - (2.0_real64*nu*timeStep)/(h**2)
  A(i, i+1) = (nu*timeStep)/(h**2) + (timeStep*nu)/(2.0_real64*h*r(i))
  source(i) = (timeStep*G) 
end do

!---------------------------------------------------------
! Insert Neumann boundary condition:

! Coefficients found from discretising Neumann boundary condition at the origin for a ghost node:
eta = (nu*timeStep)/(h**2)

! Using the Neumann b.c. to find an expression for a ghost node, 
! substituting that into the central difference scheme and rearranging
! we get the following coefficients:
A(1, 1) = 1.0_real64 - 2.0_real64*eta
A(1, 2) = 2.0_real64*eta
source(1) = timeStep*G
!---------------------------------------------------------

!---------------------------------------------------------
! Navier-slip b.c. at the boundary:
A(n, n-1) = (2.0_real64*nu*timeStep)/(h**2)
A(n, n) = 1.0_real64 - (nu*timeStep)/(l*r(n)) &
  - (2.0_real64*nu*timeStep)/(h**2) &
  - (2.0_real64*nu*timeStep)/(l*h)
source(n) = (timeStep*G)
!---------------------------------------------------------

!---------------------------------------------------------
! Printing for debugging
!print *, 'Mesh step:', h
!print *, 'Time step:', timeStep
!print *, 'Mesh:'
!do i=1, n
!  print *, r(i)
!end do
!print *, 'Source:'
!do i=1, n
!  print *, source(i)
!end do
!print *, 'Matrix:'
!do i=1, n
!  print *, A(i, :)
!end do
!print *, 'Matrix mult:'
!temp = matmul(A, f)
!do i=1, n
!  print *, temp(i)
!end do
!---------------------------------------------------------

!---------------------------------------------------------
! Initialise a column index:
col = 1
!---------------------------------------------------------

!---------------------------------------------------------
! Solve the time evolution explicitly: Forward Euler
! method for a number of iterations until it has converged 
! - If loop does not converge within a tolerance we can
! stop the program and denounce that the solution has
! diverged.
!do while (absDiff > tol .AND. loop <= maxIter) 
do t = 1, totalTimeSteps
  !temp = source + matmul(A, f)
  temp(1) = source(1) + A(1,1)*f(1) + A(1,2)*f(2)
  do i = 2, n-1
    temp(i) = source(i) + &
      A(i,i-1)*f(i-1) + A(i,i)*f(i) + A(i,i+1)*f(i+1)
  enddo
  temp(n) = source(n) + A(n, n-1)*f(n-1) + A(n, n)*f(n)
  !do i = 1, n
  !  if (abs(temp(i)) < 1.0e-9) then
  !    diffArr(i) = 0.0_real64
  !  else
  !    diffArr(i) = (temp(i) - f(i))/temp(i)
  !  end if
  !end do
  !absDiff = maxval(diffArr)
  ! Update the array
  f = temp
 
  ! Save flow for certain temporal positions:
  if (modulo(t,timeStep4Plot) == 0) then
    !----------------------------------------
    do i = 1, n
        data2Plot(i,col) = f(i)
    end do
    col = col + 1
  !---------------------------------------------------------
  end if
end do
!---------------------------------------------------------
! Calculate the absolute error and the relative error:
do i = 1, n
  absErr(i) = abs(f(i) - f_ana(i))
  if (abs(f_ana(i)) < 1.0e-9) then
    relErr(i) = 0.0_real64
  else
    relErr(i) = absErr(i)/f_ana(i)
  end if
end do
!---------------------------------------------------------

!---------------------------------------------------------
! Calculating the analyic volume flux and computing the 
! numeric volume flux:
volumeFlux_ana = (pi*G*(W**4))/(8.0_real64*nu)*(1.0_real64 + (4.0_real64*l)/W)
volumeFlux_num = 0.0_real64
do i = 1, n
  integrand(i) = f(i) * r(i)
end do
do i = 1, n-1
  volumeFlux_num = volumeFlux_num + (h/2.0_real64)*(integrand(i) + integrand(i+1))
end do
volumeFlux_num = 2.0_real64*pi*volumeFlux_num
print *, 'Analytic volume flux: ', volumeFlux_ana
print *, 'Numeric volume flux: ', volumeFlux_num
!---------------------------------------------------------

!---------------------------------------------------------
! Write the solution to a text file
open (9 , file = "pipeFlow-transient-navier.txt" , form = 'formatted')
23 FORMAT (5 ( ES23 .12 E3 ) )
! Choose format for text file ( real numbers )
do i =1 , n
  write (9 ,23) r(i) , f_ana(i) ,f(i), absErr(i), relErr(i)
end do
close (9)
!---------------------------------------------------------

!---------------------------------------------------------
open(unit=10, file='data2Plot-pipeFlow-navier.dat',status='unknown')
24 format(100(ES18.7,3X))
do i = 1, n
  write(10,24) r(i), (data2Plot(i,j), j = 1, nPlots)
enddo
close(10)
!---------------------------------------------------------

!if (loop > maxIter) then
!  print *, 'Iteration ended without convergence.'
!else 
!  print *, 'Loop iterated for ', loop, ' times.'
!  print *, 'Difference calculated: ', absDiff
!end if
end program main
