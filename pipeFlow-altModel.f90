program main
use iso_fortran_env
implicit none
!Number of nodes
!integer, parameter :: n = 51
integer :: info , i, n, fluid_int_index
integer, dimension ( : ), allocatable :: ipiv(:)
! Constants G, nu and boundary length W
real(real64), parameter :: mu_1 = (9.6_real64)*1.0e-4, rho = 730.0_real64, &
  W = (22.0_real64)*1.0e-9, pres_drop = (7.1_real64)*1.0e3, pipe_L = (78.0_real64)*1.0e-6, &
  l=(41.0_real64)*1.0e-9, alpha = (21.3_real64)*1.0e-9, pi=4.0_real64*ATAN(1.0_real64), &
  tol = (1.0_real64)*1.0e-15
real(real64) :: p_z, mu_2, enhance_fac, volumeFlux_ana, volumeFlux_num, eta, zeta
! An array of real numbers of size nxn
real(real64), dimension (: , :), allocatable :: A(:, :)
! Note solution is now contained in f and r is spatial coordinates
real(real64), dimension (:), allocatable :: r(:), b(:), f(:), f_ana(:), integrand(:), absErr(:), relErr(:)
real(real64) :: h ! Mesh step 

!---------------------------------------------------------
! Initialsing linearised pressure gradient, flow rate 
! enhancement factor, and the second viscosity:
p_z = -pres_drop/pipe_L
enhance_fac = 1.0_real64 + (4.0_real64*l)/W
mu_2 = mu_1*((W**4 - alpha**4)/(enhance_fac*(W**4) - alpha**4))
print *, 'enhancement factor: ', enhance_fac
print *, 'mu_1 = ', mu_1, 'and mu_2 = ', mu_2
!---------------------------------------------------------


!---------------------------------------------------------
! Distance between nodes , note we make ’n ’ the integer into a real number using ’ real ’
!h = W/(real(n, real64) - 1.0_real64)
! Or, hardcode the mesh-step with a certain value. Here, this is done to guarantee a node is
! present at the fluid-fluid interface of the model:
h = (1.0_real64)*1.0e-10
print *, 'W = ', W, 'and h = ', h
n = nint(W/h+1.0_real64)
print *, 'number of nodes: ', n
!---------------------------------------------------------

allocate(ipiv(1:n))
allocate(A(1:n, 1:n))
allocate(r(1:n))
allocate(b(1:n))
allocate(f(1:n))
allocate(f_ana(1:n))
allocate(integrand(1:n))
allocate(absErr(1:n))
allocate(relErr(1:n))

!---------------------------------------------------------
! Identify nodal positions x(i) and create analytic solution f_ana
do i = 1, n
  r(i) = (real(i, real64) - 1.0_real64) * h
  ! Identify what node is at the fluid boundary:
  if (abs(r(i) - alpha) < tol) then
    fluid_int_index = i
  end if
end do
do i = 1, n
  if (i <= fluid_int_index) then 
    f_ana(i) = p_z/(4.0_real64*mu_1) * (r(i)**2 - alpha**2) - p_z/(4.0_real64*mu_2) * (W**2 - alpha**2)
  else if (i > fluid_int_index .and. i <= n) then
    f_ana(i) = p_z/(4.0_real64*mu_2) * (r(i)**2 - W**2)
  end if
end do
!---------------------------------------------------------

!---------------------------------------------------------
! Initialise all values of A and b to be zero
A = 0.0_real64 ; b = 0.0_real64
!---------------------------------------------------------

eta = (2.0_real64*(h**2)*r(fluid_int_index))/(2.0_real64*r(fluid_int_index)*mu_2 - h*mu_2)
zeta = (2.0_real64*(h**2)*r(fluid_int_index))/(2.0_real64*r(fluid_int_index)*mu_1 + h*mu_1)

!---------------------------------------------------------
! Fill the matrices A and b with known values by looping through the bulk / non-boundary nodes
do i = 2 ,n-1
  if (i < fluid_int_index) then 
    A(i, i-1) = mu_1/h**2 - mu_1/(2.0_real64*h*r(i))
    A(i, i) = -(2.0_real64*mu_1)/h**2
    A(i, i+1) = mu_1/h**2 + mu_1/(2.0_real64*h*r(i))
    b(i) = p_z
  else if (i == fluid_int_index) then
    A(i, i-1) = 1.0_real64 - zeta*(mu_1/(h**2) - mu_1/(2.0_real64*h*r(i)))
    A(i, i) = ((2.0_real64*mu_2)/mu_1 - 2.0_real64) - (2.0_real64*(mu_2**2)*eta)/(mu_1*(h**2)) + (2.0_real64*mu_1*zeta)/(h**2)
    A(i, i+1) = mu_2/mu_1 * (eta*(mu_2/(h**2) + mu_2/(2.0_real64*h*r(i))) - 1.0_real64)
    b(i) = ((mu_2*eta)/mu_1 - zeta)*p_z
  else if (i > fluid_int_index .and. i <= n) then
    A(i, i-1) = mu_2/h**2 - mu_2/(2.0_real64*h*r(i))
    A(i, i) = -(2.0_real64*mu_2)/h**2
    A(i, i+1) = mu_2/h**2 + mu_2/(2.0_real64*h*r(i))
    b(i) = p_z
  end if
end do

!---------------------------------------------------------
! Insert Neumann, smoothness b.c. at origin, expressing 
! ghost node in terms of physical nodes and rearranging
! equation:
A(1, 1) = -(mu_1)/h
A(1, 2) = mu_1/h
b(1) = h/2.0_real64 * p_z
!---------------------------------------------------------

!---------------------------------------------------------
! Insert Dirichlet boundary condition:
A(n ,n) = 1.0_real64
b(n) = 0.0_real64
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
volumeFlux_ana = -(pi*p_z*(W**4))/(8.0_real64*mu_1)*enhance_fac
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
open (9 , file = "pipeFlow-altModel.txt" , form = 'formatted')
23 FORMAT (5 ( ES23 .12 E3 ) )
! Choose format for text file ( real numbers )
do i =1 , n
  write (9 ,23) r(i) , f_ana(i) ,f(i), absErr(i), relErr(i)
end do
close (9)
!---------------------------------------------------------
deallocate(ipiv); deallocate(A); deallocate(r); deallocate(b); deallocate(f);
deallocate(f_ana); deallocate(integrand); deallocate(absErr); deallocate(relErr);
end program main
