program genMesh
  implicit none
  integer :: deg, num_elem1d, num_elem, num_pt1d, num_pt, i1, i2, m1, m2, p1, p2, i, m, p, num_quad
  double precision :: L
  double precision, dimension(:), allocatable :: quad, weights
  integer, dimension(:,:), allocatable :: structure
  double precision, dimension(:,:), allocatable :: coords
  character(len=50) :: mshname

  write(*,*) "Please input interpolation degree and sidelength number of elements."
  read(*,*) deg, num_elem1d

  num_pt1d = num_elem1d*deg+1
  num_pt = (num_elem1d*deg+1)**2
  num_elem = num_elem1d**2
  num_quad = (deg+1)**2

  write(*,*) "Please input the sidelength of the square domain."
  read(*,*) L

  allocate(quad(deg+1))
  allocate(weights(deg+1))
  allocate(structure(num_elem,deg+1))
  allocate(coords(num_pt,2))

  write(*,*) "Please insert the quadrature points."
  do i = 1, deg+1
     read(*,*) quad(i)
  enddo

  write(*,*) "Please insert the weights."
  do i = 1, deg+1
     read(*,*) weights(i)
  enddo

  do m1 = 1, num_elem1d
     do m2 = 1, num_elem1d
        do i1 = 1, deg+1
           do i2 = 1, deg+1
              m = (m2-1)*num_elem1d + m1
              i = (i2-1)*(deg+1) + i1
              p1 = (m1-1)*(deg+1) + i1 - m1 + 1
              p2 = (m2-1)*(deg+1) + i2 - m2 + 1
              p = (p2-1)*num_pt1d + p1
              structure(m,i) = p
           enddo
        enddo
     enddo
  enddo

  do m1 = 1, num_elem1d
     do m2 = 1, num_elem1d
        do i1 = 1, deg+1
           do i2 = 1, deg+1
              m = (m2-1)*num_elem1d + m1
              i = (i2-1)*(deg+1) + i1
              p = structure(m,i)
              coords(p,1) = (m1-1)*(1.0d0/num_elem1d) + 0.5d0*(L/num_elem1d)*(quad(i1)+1.0d0)
              coords(p,2) = (m2-1)*(1.0d0/num_elem1d) + 0.5d0*(L/num_elem1d)*(quad(i2)+1.0d0)
           enddo
        enddo
     enddo
  enddo

  write(*,*) "Please name the mesh file."
  read(*,*) mshname

  open(unit=111,file=trim(mshname))

  write(111,*) 2, deg, num_elem, num_pt
  write(111,*) (quad(i),i=1,deg+1)
  write(111,*) (weights(i),i=1,deg+1)
  do m = 1, num_elem
     write(111,*) (structure(m,i),i=1,num_quad)
  enddo
  do p = 1, num_pt
     write(111,*) coords(p,1), coords(p,2)
  enddo
  close(111)
  write(*,*) "Thanks!"
end program genMesh
