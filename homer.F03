!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------HOMER v1.0---------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module homer
  use hesiod
  implicit none

contains
  subroutine integrate(fld, eqn, t0, tf, dt, info)
     type(fields), intent(inout) :: fld
     type(equation(fld%base%num_dim, fld%num_fields, fld%base%num_quad)), intent(in) :: eqn
     real(kind(0.0d0)), intent(in) :: t0, tf, dt
     integer, intent(out) :: info

     integer :: num_tsteps, n
     real(kind(0.0d0)), dimension(fld%num_fields, fld%base%num_quad) :: local_val

     num_tsteps = int((tf-t0)/dt)

     call MPI_INIT(ierr)
     if (ierr .ne. MPI_SUCCESS) then
        write(*,*) "Error starting MPI program. Terminating."
        call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
     endif

     call MPI_COMM_SIZE(MPI_COMM_WORLD, num, ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
     if (rank = 0) then
        do i = 1, num
           m = i
           do while (m .le. fld%base%num_elem)
              call MPI_SEND(fld%values(m,:),num_quad, MPI_DOUBLE_PRECISION, i, m, MPI_COMM_WORLD)
              m = m + num
           enddo
        enddo
        do n = 1, num_tsteps
           do i = 1, num
              m = i
              do while (m .le. fld%base%num_elem)
                 call MPI_RECV(fld%values(m,:), num_quad, MPI_DOUBLE_PRECISION, i, m, MPI_COMM_WORLD, status, ierr)
                 m = m + num
              enddo
           enddo
           call fld%projectElemValues
           do i = 1, num
              m = i
              do while (m .le. fld%base%num_elem)
                 call MPI_SEND(fld%values(m,:),num_quad,MPI_DOUBLE_PRECISION,i,m,MPI_COMM_WORLD)
                 m = m + num
              enddo
           enddo
        enddo
     else
        m = rank

        if (ierr .ne. MPI_SUCCESS) then
           write(*,*) "Error getting processor information. Terminating."
           call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        endif

        local_val = fld%elem_values(m,:,:)
        diff = fld%base%diff
        grad = fld%base%grad(m,:,:)
        t = t0
        do n = 1, num_tsteps
           do while (m .le. fld%base%num_elem)
              call MPI_RECEIVE(local_val, num_quad, MPI_DOUBLE_PRECISION, 0, m, MPI_COMM_WORLD, status, ierr)
              call euler(eqn,local_val,diff,grad,t,dt)
              call MPI_SEND(local_val, num_quad, MPI_DOUBLE_PRECISION, 0, m, MPI_COMM_WORLD)
              m = m + num
           enddo
           t = t + dt
        enddo
     endif
  end subroutine

  subroutine euler(eqn, val, diff, grad, t, dt)
     type(equation) :: eqn
     real(kind(0.0d0)), dimension(eqn%num_fields, eqn%num_quad) :: val
     double precision :: t, dt

     if ((eqn%num_split .ne. 1) .or. (eqn%theta .ne. 0.0d0)) then
        write(*,*) "ERROR: bad parameters"
     endif

     do a = 1, num_fields
        dval = eqn%operators(val, diff, grad, t)
        do i = 1, num_quad
           val_tmp = val(a,i) + dval(a,i)*dt
           val(a,i) = val_tmp
        enddo
     enddo
  end subroutine
end module homer
