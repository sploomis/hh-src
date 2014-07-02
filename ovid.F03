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

  type sparse_vec
     integer :: size, num_elem
     integer, dimension(:) :: indices
     real(kind(0.0d0)), dimension(:), allocatable :: values

  contains
     procedure :: init => sparse_vec_init
     procedure :: add => sparse_vec_add
     procedure :: rm => sparse_vec_rm
     procedure :: check => sparse_vec_check
     procedure :: elemAt => sparse_vec_elemAt
  end type

  type sparse_rows
  end type

  type sparse_cols
  end type

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

  subroutine LU_decompose_sparse(n, mat, out_mat)
     integer, intent(in) :: n
     type(sparse), intent(in) :: mat
     type(sparse), intent(out) :: out_mat

     if ((mat%num_row .ne. n) .or. (mat%num_col .ne. n)) then
        write(*,*) "Error: Matrix bounds do not match given size"
     endif

     

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

  subroutine sparse_vec_init(this, sz)
     class(sparse_vec), intent(inout) :: this
     integer, intent(in) :: sz

     this%size = sz
     this%num_elem = 0
  end subroutine sparse_init

  subroutine sparse_vec_add(this, ind, val)
     class(sparse), intent(inout) :: this
     integer, intent(in) :: ind
     real(kind(0.0d0)), intent(in) :: val
     integer :: info

     this%check(ind,val,1.0d-10,info)

     if (info = 1) then
        write(*,*) "sparse_vec_add: Could not add: outside of vector bounds"
     elseif (info = 2) then
        this%rm(ind)
        
     endif
     
  end subroutine sparse_add

  function sparse_vec_elemAt(this, ind) result(ans)
     class(sparse_vec), intent(in) :: this
     integer, dimension(2), intent(in) :: ind
     integer :: ans
     integer :: i

     ans = 0
     do i = 1, this%num_elem
        if (this%indices(i) .eq. ind)
           ans = i
        endif
     enddo
  end function sparse_vec_elemAt

  subroutine sparse_vec_check(this, ind, val, threshold, info)
     class(sparse_vec), intent(in) :: this
     integer, intent(in) :: ind
     real(kind(0.0d0)), intent(in) :: val
     real(kind(0.0d0)), optional, intent(in) :: threshold
     real(kind(0.0d0)) :: th
     integer, intent(out) :: info
     
     if (present(threshold)) then
        th = threshold
     else
        th = 1.0d-10
     endif

     info = 0

     ! Check if there are any out-of-bounds indices
     if ((ind .ge. this%size) .or. (ind .le. 0)) then
        info = 1
     endif

     ! Check if there is a zero value
     if (abs(val) .le. th) then
        info = 2
     endif
  end subroutine sparse_verify
end module ovid
