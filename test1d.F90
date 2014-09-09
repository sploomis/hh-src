!--------- Tests the code for a 1D chemical transport model:
!
!		           2                        
!		du        d u                       
!		  1          1               1      
!		---  =  d ----  -  k  u   +  - k  u 
!		dt           2      1  1     2  2  2
!		           dx                       
!
!		           2                        
!		du        d u                       
!		  2          2                      
!		---  =  d ----  +  2 k  u   -  k  u 
!		dt           2        1  1      2  2
!		           dx
!
!           with initial conditions
!
!               u (x) = 1 + cos(PI*x)
!                1
!
!		u (x) = 0
!		 2
!
!           over the region 0 < x < 1 and time 0 < t < 1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test1d
  use mpi
  use homer               !<-- Imports homer procedures
  use hesiod              !<-- Imports hesiod types
  implicit none
  type(mesh) :: base
  type(fields) :: u
  type(equation) :: eqn
  type(fn_ptr), dimension(2,1) :: operators
  real(kind(0.0d0)), dimension(:,:), allocatable :: values
  integer :: mpisize, el_per_proc, num_proc, my_id, ierr, ntsteps, nprint, i, p, m, status(MPI_STATUS_SIZE)
  integer, dimension(:,:), allocatable :: map
  real(kind(0.0d0)) :: PI, courant,tf, dt, dx
  character(len=50) :: filename, mstring

  PI = acos(-1.0d0)
  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,ierr)
  num_proc = mpisize - 1              !<-- total processors = number of elemental processors + root processor
  el_per_proc = 1
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

  m = 40                  !<-- Setting number of elements. Must have generated a mesh with these elements through genMesh1d.

  operators(1,1)%fp => RHS1   !<-- Stores the function RHS1 in operators(1,1)
  operators(2,1)%fp => RHS2   !<-- Stores the function RHS2 in operators(2,1)

  call eqn%addOps(1,operators,1,2,4)  !<-- Inserts equation information into eqn
  call eqn%addTheta(0.0d0)            !<-- Inserts implicit parameter. Only works on 0 for now.

  tf = 1.0d0                          !<-- Gives us our final time.

  courant = 0.02d0

  if (my_id .eq. 0) then
     dx = 1.0d0/m
     dt = courant*dx**2.0d0           !<-- enforces Courant factor for our second-order equation.
     ntsteps = nint(tf/dt)            !<-- gets number of timesteps
     nprint = nint((1.0d0*ntsteps)/10.0d0)  !<-- We want to print 10 files

     !--- Here we load up the base mesh from file
     write(mstring,*) m
     allocate(map(m,1))
     call base%load("meshelem"//trim(adjustl(mstring)))
     call u%setMesh(base,2)

     !--- Here we define the initial values of our field.
     allocate(values(2,base%num_pt))
     do p = 1, base%num_pt
        values(1,p) = cos(PI*base%coords(p,1))+1
        values(2,p) = 0.0d0
     enddo
     call u%init(values)
     call u%getElemValues

     !--- Here we define the map array. Since el_per_proc = 1, each row of map contains one element.
     do i = 1, m
        map(i,1) = i
        !--- We also send the number of timesteps and timestep size out to the elemental processors.
        call MPI_SEND(ntsteps,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,i,0,MPI_COMM_WORLD,ierr)
     enddo

     !--- Here we initialize the integration, sending field data out to the relevant processors.
     call homer_init(m,1,map,u,eqn)
     write(filename,*) "out"
     !--- Here we begin the integration in root, which averages boundary values and prints them to file.
     call homer_integrate_root(m,1,map,u,ntsteps,nprint,filename)
     write(*,*) "finished"
  else
     !--- Here we receive the timestep number and size from root
     call MPI_RECV(ntsteps,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(dt,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,status,ierr)

     !--- Here we begin the integration on the elemental processor
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
     real(kind(0.0d0)) :: kay1, kay2, d, diffusion, R
     real(kind(0.0d0)), dimension(num_quad,num_quad) :: lap
     integer :: a, i, j, k

     kay1 = 1.0d0
     kay2 = 0.0d-1
     d = 0.1d0

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
       
     R = -kay1*u(1,q) + 0.5d0*kay2*u(2,q)
     diffusion = 0.0d0
     do j = 1, num_quad
        diffusion = diffusion + d*lap(q,j)*u(1,j)
     enddo
     ans = diffusion + R
  end function

  function RHS2(dim, num_fields, num_quad, diff,grad,Jw, u,q, t) result(ans)
     integer, intent(in) :: dim, num_fields, num_quad,q
     real(kind(0.0d0)), dimension(num_fields,num_quad), intent(in) :: u
     real(kind(0.0d0)), dimension(dim,num_quad,num_quad), intent(in) :: diff,grad
     real(kind(0.0d0)), dimension(num_quad), intent(in) :: Jw 
     real(kind(0.0d0)) :: ans
     real(kind(0.0d0)), intent(in) :: t  
     real(kind(0.0d0)) :: kay1, kay2, d, diffusion, R
     real(kind(0.0d0)), dimension(num_quad,num_quad) :: lap
     integer :: a, i, j, k
 
     kay1 = 1.0d0
     kay2 = 0.0d-1
     d = 0.1d0

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

     R = 2.0d0*kay1*u(1,q) - kay2*u(2,q)
     diffusion = 0.0d0
     do j = 1, num_quad
        diffusion = diffusion + d*lap(q,j)*u(2,j)
     enddo
     ans = diffusion + R
  end function
end program test1d
