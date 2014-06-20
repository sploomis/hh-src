!---------------------------------------------------------!
!--------------------- HESIOD v1.0 -----------------------!
!---------------------------------------------------------!
!  Defines the types that will be used in the HOMER       !
!  evolution method.
!
!  TYPES:
!        mesh
!        fields
!        fn_ptr
!        equation
!
module hesiod
  implicit none

  type mesh(num_elem, num_pt, num_lyr, deg, num_dim)
     integer, len :: num_elem, num_pt, num_lyr, deg, num_dim
     integer, dimension(num_elem,(deg+1)**3) :: structure
     real(kind(0.0d0)), dimension(deg+1) :: quad, weights
     real(kind(0.0d0)), dimension(num_pt,num_dim) :: coords

     real(kind(0.0d0)), dimension(deg+1,deg+1) :: derivs
     real(kind(0.0d0)), dimension(num_dim,(deg+1)**3,(deg+1)**3)
  end type mesh

  type fields(num_fields, num_pt)
     integer, len :: num_fields
     real(kind(0.0d0)), dimension(num_fields, num_pt) :: values
  end type

  type fn_ptr
     procedure, pointer :: fp
  end type

  type equation(num_split, num_fields)
    integer, len :: num_split, num_fields
    type(fn_ptr), dimension(num_fields,num_split) :: Ops
    character(len=10) :: split_type
    character(len=10) :: evolve_type
  end type

contains

  subroutine mesh_calcDPhi(this)
     implicit none
     class(mesh), intent(inout) :: this
     integer :: i, j, k, l
     real(kind(0.0d0)) :: psum, prod

     do i = 1, this%deg+1
        do j = 1, this%deg+1
           psum = 0
           do k = 1, this%deg+1
              if (i .ne. k) then
                 prod = 1.0d0/(this%quad(i) - this%quad(k))
                 do l = 1, this%deg+1
                    if ((l .ne. i) .and. (l .ne. k)) then
                       prod = prod*(this%quad(j) - this%quad(k))*1.0d0/(this%quad(i) - this%quad(k))
                    endif
                 enddo
              endif
              psum = psum + prod
           enddo
           this%derivs(j,i) = psum
        enddo
     enddo
  end subroutine mesh_calcDPhi

  subroutine mesh_calcDiff(this)
     implicit none
     class(mesh), intent(inout) :: this
     integer :: a, iv, jv

     do a = 1, this%num_dim
        do iv = 1, (d+1)**3
           do jv = 1, (d+1)
           enddo
        enddo
     enddo

        
  end subroutine mesh_calcDiff

  subroutine mesh_calcJac
  end subroutine mesh_calcJac

  subroutine m
end module hesiod
