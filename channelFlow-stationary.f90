program main
use iso_fortran_env
implicit none
!Number of nodes
integer, parameter :: n = 51
integer :: info , i
integer, dimension ( n ) :: ipiv
! Constants G, nu and boundary length W
real(real64), parameter :: G = 2.0_real64, nu = 0.5_real64, W = 3.0_real64, pi=4.0_real64*ATAN(1.0_real64)
real(real64) :: volumeFlux_ana, volumeFlux_num
! An array of real numbers of size nxn
real(real64), dimension (n , n ) :: A
! Note solution is now contained in f and r is spatial coordinates
real(real64), dimension ( n ) :: y, b, f, f_ana, integrand, absErr, relErr


real(real64) :: h ! Mesh step 

!---------------------------------------------------------
! Distance between nodes , note we make ’n ’ the integer into a real number using ’ real ’
h = (2.0_real64*W)/(real(n, real64) - 1.0_real64)
!---------------------------------------------------------

!---------------------------------------------------------
! Identify nodal positions x(i) and create analytic solution f_ana
do i =1 , n
  y(i) = (real(i, real64) - 1.0_real64) * h - W
  f_ana(i) = G/(2.0_real64*nu)*(W**2 - y(i)**2)
end do
!---------------------------------------------------------

!---------------------------------------------------------
! Initialise all values of A and b to be zero
A = 0.0_real64 ; b = 0.0_real64
!---------------------------------------------------------

!---------------------------------------------------------
! Fill the matrices A and b with known values by looping through the bulk / non-boundary nodes
do i = 2 ,n-1
  A(i, i-1) = -1.0_real64/2.0_real64
  A(i, i) = 1.0_real64
  A(i, i+1) = -1.0_real64/2.0_real64
  b(i) = ((h**2)*G)/(2.0_real64*nu)
end do

!---------------------------------------------------------
! Insert Dirichlet boundary condition:
A(1, 1) = 1.0_real64
A(n ,n) = 1.0_real64
!---------------------------------------------------------

!---------------------------------------------------------
! Printing for debugging
!print *, 'Mesh step:', h
!print *, 'Mesh:'
!do i=1, n
!  print *, r(i)
!end do
!print *, 'Source:'
!do i=1, n
!  print *, b(i)
!end do
!print *, 'Matrix:'
!do i=1, n
!  print *, A(i, :)
!end do
!---------------------------------------------------------

!---------------------------------------------------------
! ‘ Call ’ the linear solver provided by Lapack
! Form of inputs and outputs can be found on Lapack website
call dgesv (n ,1 ,A ,n , ipiv ,b ,n , info )
! If dgesv subroutine fails , then the integer ’ info ’ will be non - zero .
! Then the if statement will be true and the write statement will occur
if (info .ne. 0) then
  write (* ,*) 'Inversion Failed'
else
  write (* ,*) 'Successful Inversion' ! info =0 so dgesv ran as planned
end if
! dgesv outputs the solution to b , which we assign to f for neatness
f=b
!---------------------------------------------------------

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
volumeFlux_ana = (2.0_real64*(W**3)*G)/(3.0_real64*nu)
volumeFlux_num = 0.0_real64
do i = 1, n
  integrand(i) = f(i)
end do
do i = 1, n-1
  volumeFlux_num = volumeFlux_num + (h/2.0_real64)*(integrand(i) + integrand(i+1))
end do
volumeFlux_num = volumeFlux_num
print *, 'Analytic volume flux: ', volumeFlux_ana
print *, 'Numeric volume flux: ', volumeFlux_num
!---------------------------------------------------------

!---------------------------------------------------------
! Write the solution to a text file
open (9 , file = "channelFlow-stationary.txt" , form = 'formatted')
23 FORMAT (5 ( ES23 .12 E3 ) )
! Choose format for text file ( real numbers )
do i =1 , n
  write (9 ,23) y(i) , f_ana(i) ,f(i), absErr(i), relErr(i)
end do
close (9)
!---------------------------------------------------------
end program main
