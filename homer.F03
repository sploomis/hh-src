!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------HOMER v1.0---------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module homer
  use mpi
  use hesiod
  implicit none

contains

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
              val(j,:,:) = fld%elem_values(m,:,:)
              grad(j,:,:,:) = fld%base%grad(m,:,:,:)
              Jw(j,:) = fld%base%Jw(m,:)
           endif
        enddo
        call MPI_SEND(val,size(val),MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,ierr)
        call MPI_SEND(grad,size(grad),MPI_DOUBLE_PRECISION,i,3,MPI_COMM_WORLD,ierr)
        call MPI_SEND(Jw,size(Jw),MPI_DOUBLE_PRECISION,i,4,MPI_COMM_WORLD,ierr)
     enddo
  end subroutine

  subroutine homer_integrate(el_per_proc,eqn,t0,dt,ntsteps)
     type(equation), intent(in) :: eqn
     integer, intent(in) :: el_per_proc, ntsteps
     real(kind(0.0d0)), intent(in) :: t0, dt

     integer, dimension(el_per_proc) :: elem
     real(kind(0.0d0)), dimension(eqn%dim,eqn%num_quad,eqn%num_quad) :: diff
     real(kind(0.0d0)), dimension(el_per_proc,eqn%dim,eqn%num_quad,eqn%num_quad) :: grad
     real(kind(0.0d0)), dimension(el_per_proc,eqn%num_quad) :: Jw
     real(kind(0.0d0)), dimension(el_per_proc,eqn%num_fields,eqn%num_quad) :: val

     integer :: ierr, status(MPI_STATUS_SIZE), m, j, n
     real(kind(0.0d0)) :: t

     ! Receive data from homer_init

     call MPI_RECV(elem,el_per_proc,MPI_INTEGER,0,0,MPI_COMM_WORLD,status,ierr)

     call MPI_RECV(diff,size(diff),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(val,size(val),MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(grad,size(grad),MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,status,ierr)
     call MPI_RECV(Jw,size(Jw),MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,status,ierr)

     ! Integrate data, communicating with homer_project

     write(*,*) "Got this far!"

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

  subroutine homer_project(num_proc,el_per_proc,map,fld,ntsteps)
     integer, intent(in) :: num_proc,el_per_proc,ntsteps
     integer, dimension(num_proc,el_per_proc), intent(in) :: map
     type(fields), intent(inout) :: fld

     real(kind(0.0d0)), dimension(el_per_proc,fld%num_fields,fld%base%num_quad) :: val
     integer :: i, j, n, m, status(MPI_STATUS_SIZE), ierr
     
     do n = 1, ntsteps
        do i = 1, num_proc
           call MPI_RECV(val,size(val),MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,status,ierr)
           do j = 1, el_per_proc
              m = map(i,j)
              if (m .ne. -1) then
                 fld%elem_values(m,:,:) = val(j,:,:)
              endif
           enddo
        enddo

        call fld%projectElemValues

        do i = 1, num_proc
           do j = 1, el_per_proc
              m = map(i,j)
              if (m .ne. -1) then
                 val(j,:,:) = fld%elem_values(m,:,:)
              endif
           enddo
           call MPI_SEND(val,size(val),MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,status,ierr)
        enddo
     enddo
  end subroutine

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

  subroutine homer_euler(eqn,s,diff,grad,Jw,val,t,dt)
     type(equation), intent(in) :: eqn
     integer, intent(in) :: s
     real(kind(0.0d0)), dimension(eqn%dim,eqn%num_quad,eqn%num_quad), intent(in) :: diff, grad
     real(kind(0.0d0)), dimension(eqn%num_quad), intent(in) :: Jw
     real(kind(0.0d0)), dimension(eqn%num_fields,eqn%num_quad), intent(inout) :: val
     real(kind(0.0d0)), intent(in) :: t, dt

     real(kind(0.0d0)), dimension(eqn%num_fields,eqn%num_quad) :: val_tmp
     integer :: i, f

     do f = 1, eqn%num_fields
        val_tmp(f,:) = val(f,:) + eqn%operators(f,s)%fp(eqn%dim,eqn%num_fields,eqn%num_quad,diff,grad,Jw,val(f,:),t)
     enddo
     do f = 1, eqn%num_fields
        val(f,:) = val_tmp(f,:)
     enddo
  end subroutine
end module homer
