program genMesh
  implicit none
  integer :: deg, num_elem, num_pt, i, m, p
  double precision :: L
  double precision, dimension(:), allocatable :: quad, weights
  integer, dimension(:,:), allocatable :: structure
  double precision, dimension(:,:), allocatable :: coords
  character(len=50) :: mshname

  write(*,*) "Please input interpolation degree and number of elements."
  read(*,*) deg, num_elem

  num_pt = num_elem*deg+1

  write(*,*) "Please input the length of the periodic domain."
  read(*,*) L

  allocate(quad(deg+1))
  allocate(weights(deg+1))
  allocate(structure(num_elem,deg+1))
  allocate(coords(num_pt,1))

  write(*,*) "Please insert the quadrature points."
  do i = 1, deg+1
     read(*,*) quad(i)
  enddo

  write(*,*) "Please insert the weights."
  do i = 1, deg+1
     read(*,*) weights(i)
  enddo

  do m = 1, num_elem
     do i = 1, deg+1
        structure(m,i) = (m-1)*(deg+1) + i - m + 1
     enddo
  enddo

  do m = 1, num_elem
     do i = 1, deg
        p = structure(m,i)
        coords(p,1) = (m-0.5d0+0.5d0*quad(i))*(L/num_elem)
     enddo
  enddo
  coords(num_pt,1) = L

  write(*,*) "Please name the mesh file."
  read(*,*) mshname

  open(unit=111,file=trim(mshname))

  write(111,*) 1, deg, num_elem, num_pt
  write(111,*) (quad(i),i=1,deg+1)
  write(111,*) (weights(i),i=1,deg+1)
  do m = 1, num_elem
     write(111,*) (structure(m,i),i=1,deg+1)
  enddo
  do p = 1, num_pt
     write(111,*) coords(p,1)
  enddo
  close(111)
  write(*,*) "Thanks!"
end program genMesh
