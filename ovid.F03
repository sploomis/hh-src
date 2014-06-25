module ovid
  implicit none

  type sparse
     integer :: num_row, num_col, num_elem
     integer, dimension(:,:), allocatable :: indices
     real(kind(0.0d0)), dimension(:), allocatable :: values

  contains
     procedure :: init => sparse_init
     procedure :: add => sparse_add
     procedure :: addOne => sparse_addOne
     procedure :: rm => sparse_rm
     procedure :: verify => sparse_verify
     procedure :: elemAt => sparse_elemAt
  end type sparse

contains

  subroutine LU_decompose(n, mat, out_mat)
     integer, intent(in) :: n
     real(kind(0.0d0)), dimension(n,n), intent(in) :: mat
     real(kind(0.0d0)), dimension(n,n), intent(out) :: out_mat
     real(kind(0.0d0)) :: tmp
     real(kind(0.0d0)), dimension(n) :: l
     integer :: m, i, j

     out_mat = mat

     do m = 1, n-1
        do i = m+1, n
           l(i) = -out_mat(i,m)/out_mat(m,m)
           out_mat(i,m) = -l(i)
           do j = m+1, n
              tmp = out_mat(i,j)
              out_mat(i,j) = tmp + out_mat(m,j)*l(i)
           enddo
        enddo
     enddo
  end subroutine

  subroutine mat_inv(n, mat, inv)
     integer, intent(in) :: n
     real(kind(0.0d0)), dimension(n,n), intent(in) :: mat
     real(kind(0.0d0)), dimension(n,n), intent(out) :: inv
     real(kind(0.0d0)), dimension(n,n) :: LU
     real(kind(0.0d0)), dimension(n) :: y
     real(kind(0.0d0)) :: ysum, asum
     integer :: i,j,k

     call LU_decompose(n, mat, LU)

     do k = 1, n
        do i = 1, k-1
           y(i) = 0.0d0
        enddo
        y(k) = 1.0d0
        do i = k+1, n
           ysum = 0.0d0
           do j = k, i-1
              ysum = ysum - LU(i,j)*y(j)
           enddo
           y(i) = ysum
        enddo
        inv(n,k) = y(n)/LU(n,n)
        do i = n-1, 1, -1
           asum = y(i)
           do j = i+1, n
              asum = asum - LU(i,j)*inv(j,k)
           enddo
           inv(i,k) = asum/LU(i,i)
        enddo
     enddo
  end subroutine mat_inv

  recursive function det(n, mat) result(res)
     real(kind(0.0d0)) :: res
     integer, intent(in) :: n
     real(kind(0.0d0)), dimension(n,n) :: mat
     real(kind(0.0d0)), dimension(n-1,n-1) :: mat_nxt
     real(kind(0.0d0)) :: cf
     integer :: i

     if (n .ge. 1) then
        mat_nxt(1:(n-1),:) = mat(2:n,2:n)
        cf = det(n-1,mat_nxt)
        res = cf

        do i = 2, n-1
           mat_nxt(1:(i-1),:) = mat(1:(i-1),2:n)
           mat_nxt(i:(n-1),:) = mat((i+1):n,2:n)
           cf = (-1.0d0)**(i-1)*det(n-1,mat_nxt)
           res = res + cf
        enddo

        mat_nxt(1:(n-1),:) = mat(1:(n-1),2:n)
        cf = (-1.0d0)**(n-1)*det(n-1,mat_nxt)
        res = res + cf
        return
     else
        res = mat(1,1)
        return
     endif
  end function det

  subroutine sparse_init(this, num_row, num_col)
     class(sparse), intent(inout) :: this
     integer, intent(in) :: num_row, num_col

     this%num_row = num_row
     this%num_col = num_col
     this%num_elem = 0
  end subroutine sparse_init

  subroutine sparse_add(this, num_add, indices, values)
     class(sparse), intent(inout) :: this
     integer, intent(in) :: num_add
     integer, dimension(2, num_add), intent(in) :: indices
     real(kind(0.0d0)), dimension(num_add), intent(in) :: values
     integer, dimension(:,:), allocatable :: new_indices, tmp_indices
     integer, dimension(:), allocatable :: new_values, tmp_values

     if (this%num_elem .eq. 0) then
        this%num_elem = 1
        allocate(this%indices(2,1))
        allocate(this%values(1))
        this%indices(:,1) = indices(:,1)
        this%values(1) = values(1)
        if (num_add .ge. 1) then
           do i = 2, num_add
              this%addOne(indices(:,i),values(i))
           enddo
        endif
     else
        do i = 1, num_add
           this%addOne(indices(:,i),values(i))
        enddo
     endif
  end subroutine sparse_add

  subroutine sparse_addOne(this, ind, val)
     class(sparse), intent(inout) :: this
     integer, dimension(2), intent(in) :: ind
     real(kind(0.0d0)), intent(in) :: val
     integer, dimension(:,:), allocatable :: new_indices
     real(kind(0.0d0)), dimension(:), allocatable :: new_values
     integer :: info

     call this%verify(ind, val, info)
     if (info .eq. 0) then
        i = thi%elemAt(ind)

        if (i .ne. 0) then
           this%values(i) = val
        else
           allocate(new_indices(2,this%num_elem+1))
           allocate(values(this%num_elem+1))

           do j = 1, this%num_elem
              do k = 1, 2
                 new_indices(k,j) = this%indices(k,j)
              enddo
              new_values(j) = this%values(j)
           enddo

           deallocate(this%indices)
           deallocate(this%values)
           this%num_elem = this%num_elem + 1
           allocate(this%indices(2,this%num_elem))
           allocate(this%values(this%num_elem))

           do k = 1, 2
              new_indices(k,num_elem) = ind(k)
           enddo

           new_values(num_elem) = val
           this%indices = new_indices
           this%values = new_values
        endif
     elseif (info .eq. 1) then
        write(*,*) "Error: Indices are out of bounds"
     else
        write(*,*) "Error: Can't enter a zero value"
     endif
  end subroutine sparse_addOne
            

  function sparse_elemAt(this, ind) result(ans)
     class(sparse), intent(in) :: this
     integer, dimension(2), intent(in) :: ind
     integer :: ans
     integer :: i

     ans = 0
     do i = 1, this%num_elem
        if (this%elements(:,i) .eq. ind)
           ans = i
        endif
     enddo
  end function sparse_elemAt

  subroutine sparse_verify(this, ind, val, threshold, info)
     class(sparse), intent(in) :: this
     integer, dimension(2), intent(in) :: ind
     real(kind(0.0d0)), intent(in) :: val
     real(kind(0.0d0)), optional, intent(in) :: threshold
     real(kind(0.0d0)) :: th
     integer, intent(out) :: info
     
     if (present(threshold)) then
        th = threshold
     else
        th = 1.0d-5
     endif

     info = 0

     ! Check if there are any out-of-bounds indices
     if ((ind(1) .ge. this%num_row) .or. (ind(1) .le. 0)) then
        info = 1
     endif
     if ((ind(2) .ge. this%num_col) .or. (ind(2) .le. 0)) then
        info = 1
     endif

     ! Check if there is a zero value
     if (abs(val) .le. th) then
        info = 2
     endif
  end subroutine sparse_verify
end module ovid
