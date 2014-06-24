module ovid
  implicit none

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

end module ovid
