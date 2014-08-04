!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------HOMER v1.0---------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module homer
  use mpi
  use hesiod
  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HOMER_INIT
  !
  ! The subroutine called to initialize integration on each processor. Called from root.
  !
  ! TASKS:	Send out map, mesh and initial data to each processor.
  !
  ! ARGUMENTS:	num_proc	integer		in	Number of elemental processes (MPI_COMM_SIZE - 1)
  !		el_per_proc	integer		in	Elements per processor
  !		map		int array	in	Table which stores which elements correspond to each processors
  !		fld		fields		in	Stores data about field values
  !		eqn		equation	in	Stores data about equation being integrated
  !
  ! OTHER VARIABLES:	i,j,m		integer		Dummy loop variables
  !			ierr		integer		MPI error return variable
  !			val		double array	Packages field values for sending
  !			grad		double array	Packages grad matrix for sending
  !			diff		double array    Packages diff matrix for sending
  !			Jw		double array	Packages mass matrix for sending

  subroutine homer_init(num_proc,el_per_proc,map,fld,eqn)
     integer, intent(in) :: num_proc, el_per_proc
     integer, dimension(num_proc,el_per_proc), intent(in) :: map
     type(fields), intent(in) :: fld
     type(equation), intent(in) :: eqn

     integer :: i, j, m

     integer :: ierr

     real(kind(0.0d0)), dimension(el_per_proc,eqn%num_fields,eqn%num_quad) :: val
     real(kind(0.0d0)), dimension(el_per_proc,eqn%dim,eqn%num_quad,eqn%num_quad) :: grad
     real(kind(0.0d0)), dimension(el_per_proc,eqn%num_quad) :: Jw

     do i = 1, num_proc
        call MPI_SEND(map(i,:),el_per_proc,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
        call MPI_SEND(fld%base%diff(:,:,:),size(fld%base%diff),MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,ierr)

        do j = 1, el_per_proc
           m = map(i,j)
           if (m .ne. -1) then
              val(j,:,:) = fld%elem_values(:,m,:)
              grad(j,:,:,:) = fld%base%grad(m,:,:,:)
              Jw(j,:) = fld%base%Jw(m,:)
           endif
        enddo
        call MPI_SEND(val,size(val),MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,ierr)
        call MPI_SEND(grad,size(grad),MPI_DOUBLE_PRECISION,i,3,MPI_COMM_WORLD,ierr)
        call MPI_SEND(Jw,size(Jw),MPI_DOUBLE_PRECISION,i,4,MPI_COMM_WORLD,ierr)
     enddo
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HOMER_INTEGRATE_ELEM
  !
  ! The subroutine called on the elemental processors for time integration.
  !
  ! TASKS: 	Receive relevant initial and mesh data, integrate elements, send
  !	   	discontinuous result to root, and receive averaged result back.
  !
  ! ARGUMENTS:  el_per_proc	integer		in	Number of elements per processor
  !		eqn		equation	in	Stores data about equation being integrated
  !		t0		double		in	Initial time
  !		dt		double		in	Size of timestep
  !		ntsteps		integer		in	Number of timesteps
  !		my_id		integer		in	ID number of processor
  !
  ! OTHER VARIABLES:	elem		int array		Stores map data for processor
  !			diff		double array		Stores diff matrix for processor's elements
  !			grad		double array		Stores grad matrix for processor's elements
  !			Jw		double array		Stores mass matrix for processor's elements
  !			val		double array		Stores field values for processor's elements
  !			ierr		integer			MPI error return variable
  !			status		int array		MPI status array
  !			m,j,n		integer			Dummy loop variables
  !			t		double			Time at each integration step

  subroutine homer_integrate_elem(el_per_proc,eqn,t0,dt,ntsteps,my_id)
     type(equation), intent(in) :: eqn
     integer, intent(in) :: el_per_proc, ntsteps, my_id
     real(kind(0.0d0)), intent(in) :: t0, dt

     integer, dimension(el_per_proc) :: elem
     real(kind(0.0d0)), dimension(eqn%dim,eqn%num_quad,eqn%num_quad) :: diff
     real(kind(0.0d0)), dimension(el_per_proc,eqn%dim,eqn%num_quad,eqn%num_quad) :: grad
     real(kind(0.0d0)), dimension(el_per_proc,eqn%num_quad) :: Jw
     real(kind(0.0d0)), dimension(el_per_proc,eqn%num_fields,eqn%num_quad) :: val

     integer :: ierr, status(MPI_STATUS_SIZE), m, j, n, f
     character(len=20) :: idstring, nstring
     real(kind(0.0d0)) :: t

     ! Receive data from homer_init

     call MPI_RECV(elem,el_per_proc,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(diff,size(diff),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(val,size(val),MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(grad,size(grad),MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(Jw,size(Jw),MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,status,ierr)

     ! Integrate data, communicating with homer_project

     t = t0

     do n = 1, ntsteps
        do j = 1, el_per_proc
           m = elem(j)
           if (m .ne. -1) then
              call homer_split(eqn,diff,grad(j,:,:,:),Jw(j,:),val(j,:,:),t,dt)
           endif
        enddo
        t = t + dt
        call MPI_SEND(val,size(val),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
        call MPI_RECV(val,size(val),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,status,ierr)
     enddo
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HOMER_INTEGRATE_ROOT
  !
  ! The subroutine called on the root processor for time integration.
  !
  ! TASKS: 	Print out data, receive discontinuous solutions from processors,
  !        	and perform weighted averages at elemental boundaries to make
  !        	the solution continuous.
  !
  ! ARGUMENTS:	num_proc	integer		in	Number of elemental processes (size(MPI_COMM_WORLD) - 1)
  !		el_per_proc	integer		in	Number of elements per processor
  !		map		int array	in	Table which stores which elements correspond to which processors
  !		fld		fields		inout	Stores data about the field values at each timestep
  !		ntsteps		integer		in	Number of integration steps
  !		nprint		integer		in	How often to print data (e.g. nprint = 2 means print at n=0,2,4,6,...)
  !										(nprint <= 0 results in no print)
  !		filename	string		in	Preface name of print files (e.g. filename = "out" means print files are
  !										    "out-0", "out-1", etc.)
  !
  ! OTHER VARIABLES:	val			double array	Dummy array for receiving field values
  !			i,j,k,n,m,p,a,f		integer		Dummy loop variables
  !			ierr			integer		MPI error return variable
  !			status			int array	MPI status array
  ! 			nstring			string		Dummy string for file names
                                                                                     
  subroutine homer_integrate_root(num_proc,el_per_proc,map,fld,ntsteps,nprint,filename)
     integer, intent(in) :: num_proc,el_per_proc,ntsteps, nprint
     integer, dimension(num_proc,el_per_proc), intent(in) :: map
     type(fields), intent(inout) :: fld
     character(len=50), intent(in) :: filename

     real(kind(0.0d0)), dimension(el_per_proc,fld%num_fields,fld%base%num_quad) :: val
     integer :: i, j, n, m, status(MPI_STATUS_SIZE), ierr, p, a, f
     character(len=20) :: nstring

     do n = 1, ntsteps

        ! Write most recent data to file

        if (nprint .ge. 1) then
           if (mod(n-1,nprint) .eq. 0) then
              write(nstring,*) n-1
              nstring = adjustl(nstring)
              open(1,file=trim(adjustl(filename))//"-"//trim(nstring),recl=500)
              do p = 1, fld%base%num_pt
                 write(1,*) p, (fld%base%coords(p,a),a=1,fld%base%dim), (fld%values(f,p),f=1,fld%num_fields)
              enddo
              write(*,*) "Timestep n="//trim(nstring)//" written to file."
              close(1)
           endif
        endif

        ! Receive elemental values of fields

        do i = 1, num_proc
           call MPI_RECV(val,size(val),MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,status,ierr)
           do j = 1, el_per_proc
              m = map(i,j)
              if (m .ne. -1) then
                 fld%elem_values(:,m,:) = val(j,:,:)
              endif
           enddo
        enddo

        ! Perform weighted average on element boundaries

        call fld%projectElemValues

        ! Send corrected data to elemental processors

        do i = 1, num_proc
           do j = 1, el_per_proc
              m = map(i,j)
              if (m .ne. -1) then
                 val(j,:,:) = fld%elem_values(:,m,:)
              endif
           enddo
           call MPI_SEND(val,size(val),MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,status,ierr)
        enddo
     enddo

     ! Print final data

     if (nprint .ge. 1) then
        if (mod(ntsteps,nprint) .eq. 0) then
           write(nstring,*) ntsteps
           nstring = adjustl(nstring)
           open(1,file=trim(adjustl(filename))//"-"//trim(nstring),recl=500)
           do p = 1, fld%base%num_pt
              write(1,*) p, (fld%base%coords(p,a),a=1,fld%base%dim), (fld%values(f,p),f=1,fld%num_fields)
           enddo
           write(*,*) "Timestep n="//trim(nstring)//" written to file."
           close(1)
        endif
     endif
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HOMER_SPLIT
  !
  ! TASK: 	Eventually, evolve data forwards with splitting methods. 
  !		Currently only works with trivial case (no splitting).
  !
  ! ARGUMENTS:	eqn	equation	Stores data about equation being integrated
  !		diff	double array	Stores diff matrix for current element
  !		grad	double array	Stores grad matrix for current element
  !		Jw	double array	Stores mass matrix for current element
  !		val	double array	Stores field values for current element
  !		t	double		Current time
  !		dt	double		Size of time step

  subroutine homer_split(eqn,diff,grad,Jw,val,t,dt)
     type(equation), intent(in) :: eqn
     real(kind(0.0d0)), dimension(eqn%dim,eqn%num_quad,eqn%num_quad), intent(in) :: diff, grad
     real(kind(0.0d0)), dimension(eqn%num_quad) :: Jw
     real(kind(0.0d0)), dimension(eqn%num_fields,eqn%num_quad) :: val
     real(kind(0.0d0)) :: t, dt

     if (eqn%num_split .eq. 1) then
        call homer_euler(eqn,1,diff,grad,Jw,val,t,dt)
     else
        write(*,*) "Not ready for higher-level splits yet! Sorry. - sploomis"
     endif
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HOMER_EULER
  !
  ! TASK:	Eventually, perform implicit Euler integration for a given theta parameter.
  !		Currently only works with explicit Euler (theta = 0).
  !
  ! ARGUMENTS:	eqn	equation	Stores data about equation being integrated
  !		s	integer		Splitting step
  !		diff	double array	Stores diff matrix for current element
  !		grad	double array	Stores grad matrix for current element
  !		Jw	double array	Stores mass matrix for current element
  !		val	double array	Stores field values for current element
  !		t	double		Current time
  !		dt	double		Size of time step

  subroutine homer_euler(eqn,s,diff,grad,Jw,val,t,dt)
     type(equation), intent(in) :: eqn
     integer, intent(in) :: s
     real(kind(0.0d0)), dimension(eqn%dim,eqn%num_quad,eqn%num_quad), intent(in) :: diff, grad
     real(kind(0.0d0)), dimension(eqn%num_quad), intent(in) :: Jw
     real(kind(0.0d0)), dimension(eqn%num_fields,eqn%num_quad), intent(inout) :: val
     real(kind(0.0d0)), intent(in) :: t, dt

     real(kind(0.0d0)), dimension(eqn%num_fields,eqn%num_quad) :: val_tmp
     real(kind(0.0d0)), dimension(eqn%num_quad) :: dval
     integer :: i,j, f
     

     character(len=20) tstring
 
     do f = 1, eqn%num_fields
        do i = 1, eqn%num_quad
           dval(i) = eqn%operators(f,s)%fp(eqn%dim,eqn%num_fields,eqn%num_quad,diff,grad,Jw,val,i,t)*dt
           val_tmp(f,i) = val(f,i) + dval(i)
        enddo
     enddo

     do f = 1, eqn%num_fields
        val(f,:) = val_tmp(f,:)
     enddo
  end subroutine
end module homer
