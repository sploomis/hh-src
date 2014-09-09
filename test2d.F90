!--------- Tests the code for a 1D chemical transport model:
! 
!		du                                                
!		  1         2                  1
!		---  =  d ∇  u   -  k  u   +  - k  u 
!		dt            1      1  1     2  2  2
!
!		du                                       
!		  2         2                         
!		---  =  d ∇  u   +  2 k  u   -  k  u 
!		dt            2        1  1      2  2
! 
!           with initial conditions
!
!                       { 0.25 (1 - cos(4*PI*x)) (1 - cos(4*PI*y))	0 < x < 0.5, 0 < y < 0.5
!               u (x) = { 
!                1      { 0						otherwise
! 
!                       { 0.25 (1 - cos(4*PI*x)) (1 - cos(4*PI*y))	0 < x < 0.5, 0.5 < y < 1
!               u (x) = {
!                2      { 0						otherwise
! 
!           over the region 0 < x < 1, 0 < y < 1 and time 0 < t < 10.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test2d
  use mpi
  use homer             !<-- Imports homer procedures
  use hesiod            !<-- Imports hesiod types
  implicit none
  type(mesh) :: base
  type(fields) :: u
  type(equation) :: eqn
  type(fn_ptr), dimension(2,1) :: operators
  real(kind(0.0d0)), dimension(:,:), allocatable :: values
  integer :: mpisize, el_per_proc, num_proc, my_id, ierr, ntsteps, nprint, i, p, m,m1, status(MPI_STATUS_SIZE)
  integer, dimension(:,:), allocatable :: map
  real(kind(0.0d0)) :: PI, courant,tf, dt, dx, x, y
  character(len=50) :: filename, mstring

  PI = acos(-1.0d0)
  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,ierr)
  num_proc = mpisize - 1              !<-- total processors = number of elemental processors + root processor
  el_per_proc = 5
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

  operators(1,1)%fp => RHS1           !<-- Stores the function RHS1 in operators(1,1)
  operators(2,1)%fp => RHS2           !<-- Stores the function RHS2 in operators(2,1)

  call eqn%addOps(1,operators,2,2,16) !<-- Inserts equation information into eqn
  call eqn%addTheta(0.0d0)            !<-- Inserts implicit parameter. Only works on 0 for now.

  m = 25                              !<-- Total number of elements in the square
  m1 = 5                              !<-- Number of elements on one side of the square

  tf = 10.0d0                         !<-- Gives us our final time.

  if (my_id .eq. 0) then
     dx = 1.0d0/m1
     dt = 0.0001
     ntsteps = 1.0d7
     nprint = 1.0d6
     write(*,*) ntsteps, nprint
     allocate(map(m,1))

     !--- Here we load up the base mesh from file.
     call base%load("testmesh2D")        !<--- This is different from test1d! Make sure that there is a mesh file named 'testmesh2D' and make sure it is consistent with the number of elements m.
     call u%setMesh(base,2)

     !--- Here we define the initial values of our field.
     allocate(values(2,base%num_pt))
     do p = 1, base%num_pt
        x = base%coords(p,1)
        y = base%coords(p,2)
        if ((x .le. 0.5) .and. (y .le. 0.5)) then
           values(1,p) = 0.25d0*(1.0d0-cos(4.0d0*PI*x))*(1.0d0-cos(4.0d0*PI*y))
        else
           values(1,p) = 0.0d0
        endif
        if ((x .le. 0.5) .and. (y .ge. 0.5)) then
           values(2,p) = 0.25d0*(1.0d0-cos(4.0d0*PI*x))*(1.0d0-cos(4.0d0*PI*y))
        else
           values(2,p) = 0.0d0
        endif
     enddo
     call u%init(values)
     call u%getElemValues

     !--- Here we define the map array. Since el_per_proc = 5, each row of the map contains five elements.
     do i = 1, m/el_per_proc
        do j = 1, el_per_proc
           map(i,j) = (i-1)*m1 + j
        enddo
        !--- We also send the number of timesteps and timestep size out to the elemental processors.
        call MPI_SEND(ntsteps,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,i,0,MPI_COMM_WORLD,ierr)
     enddo

     !--- Here we initialize the integration
     call homer_init(m,1,map,u,eqn)
     write(filename,*) "out"
     !--- Here we begin the integration in root
     call homer_integrate_root(m,1,map,u,ntsteps,nprint,filename)
     write(*,*) "finished"
  else
     !--- Here we receive the timestep number and size from root
     call MPI_RECV(ntsteps,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)

     !--- Here we begin the integration on the elemental processor
     call MPI_RECV(dt,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,status,ierr)
     call homer_integrate_elem(1,eqn,0.0d0,dt,ntsteps,my_id)
  endif
  call MPI_FINALIZE(ierr)
contains

  !--- These procedures store the equation information in functions matching the format required in hesiod.F90

  function RHS1(dim, num_fields, num_quad, diff,grad,Jw, u,q, t) result(ans)
     integer, intent(in) :: dim, num_fields, num_quad, q
     real(kind(0.0d0)), dimension(num_fields,num_quad), intent(in) :: u
     real(kind(0.0d0)), dimension(dim,num_quad,num_quad), intent(in) :: diff,grad
     real(kind(0.0d0)), dimension(num_quad), intent(in) :: Jw
     real(kind(0.0d0)) :: ans
     real(kind(0.0d0)), intent(in) :: t 
     real(kind(0.0d0)) :: kay1, kay2, d, diffusion, R, convection, vx, vy
     real(kind(0.0d0)), dimension(num_quad,num_quad) :: lap
     integer :: a, i, j, k

     kay1 = 1.0d0
     kay2 = 2.0d0
     d = 6.25d-3
     vx = 0.0d-1
     vy = 0.0d0

     do i = 1, num_quad
        do j = 1, num_quad
           diffusion = 0.0d0
           do k = 1, num_quad
              do a = 1, dim
                 diffusion = diffusion - Jw(k)*grad(a,k,i)*grad(a,k,j)
              enddo
           enddo
           lap(i,j) = diffusion/Jw(i)
        enddo
     enddo
       
     R = -kay1*u(1,q) + 0.5d0*kay2*u(2,q)**2.0d0
     diffusion = 0.0d0
     convection = 0.0d0
     do j = 1, num_quad
        diffusion = diffusion + d*lap(q,j)*u(1,j)
        convection = convection - (vx*grad(1,q,j) + vy*grad(2,q,j))*u(1,j)
     enddo
     ans = diffusion + convection + R
 end function

  function RHS2(dim, num_fields, num_quad, diff,grad,Jw, u,q, t) result(ans)
     integer, intent(in) :: dim, num_fields, num_quad,q
     real(kind(0.0d0)), dimension(num_fields,num_quad), intent(in) :: u
     real(kind(0.0d0)), dimension(dim,num_quad,num_quad), intent(in) :: diff,grad
     real(kind(0.0d0)), dimension(num_quad), intent(in) :: Jw 
     real(kind(0.0d0)) :: ans
     real(kind(0.0d0)), intent(in) :: t  
     real(kind(0.0d0)) :: kay1, kay2, d, diffusion, convection, R, vx, vy
     real(kind(0.0d0)), dimension(num_quad,num_quad) :: lap
     integer :: a, i, j, k
 
     kay1 = 1.0d0
     kay2 = 2.0d0
     d = 6.25d-3
     vx = 0.0d-1
     vy = 0.0d0

     do i = 1, num_quad
        do j = 1, num_quad
           diffusion = 0.0d0
           do k = 1, num_quad
              do a = 1, dim
                 diffusion = diffusion - Jw(k)*grad(a,k,i)*grad(a,k,j)
              enddo
           enddo
           lap(i,j) = diffusion/Jw(i)
        enddo
     enddo

     R = 2.0d0*kay1*u(1,q) - kay2*u(2,q)**2.0d0
     diffusion = 0.0d0
     do j = 1, num_quad
        diffusion = diffusion + d*lap(q,j)*u(2,j)
     enddo
     convection = 0.0d0
     do j = 1, num_quad
        convection = convection - (vx*grad(1,q,j) + vy*grad(2,q,j))*u(2,j)
     enddo
     ans = diffusion + convection + R
  end function
end program test2d
