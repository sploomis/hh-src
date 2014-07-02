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

!!!!!!!!!!!!!!!!!!!!!
! TYPES DEFINITIONS !
!!!!!!!!!!!!!!!!!!!!!

  ! type :: mesh
  !
  ! Contains data about the mesh on which the problem is defined

  type mesh
     ! BASIC PROPERTIES
     integer :: num_elem, num_pt, num_lyr, deg, dim, num_quad
     integer, dimension(:,:), allocatable :: structure 				!num_elem, num_quad
     real(kind(0.0d0)), dimension(:), allocatable :: quad, weights			!(deg+1)
     real(kind(0.0d0)), dimension(:,:), allocatable :: coords			!num_pt,dim

     ! DERIVED PROPERTIES
     integer :: num_bd
     integer, dimension(:), allocatable :: bd						!num_bd
     real(kind(0.0d0)), dimension(:), allocatable :: weightProds			!num_quad
     real(kind(0.0d0)), dimension(:,:), allocatable :: derivs			!deg+1,deg+1
     real(kind(0.0d0)), dimension(:,:,:), allocatable :: diffs			!dim, num_quad, num_quad
     real(kind(0.0d0)), dimension(:,:,:,:), allocatable :: jac_mat, jac_inv 		!num_elem, dim, dim, num_quad
     real(kind(0.0d0)), dimension(:,:,:,:), allocatable :: grad 			!num_elem, dim, num_quad, num_quad
     real(kind(0.0d0)), dimension(:,:), allocatable :: jac				!num_elem, num_quad

  contains
     ! PROCEDURES
     procedure :: setParam => mesh_setParam
     procedure :: setForm => mesh_setForm
     procedure :: getCoords => mesh_getCoords
     procedure :: getIndex => mesh_getIndex
     procedure :: calcDPhi => mesh_calcDPhi
     procedure :: calcDiff => mesh_calcDiff
     procedure :: calcJacMat => mesh_calcJacMat
     procedure :: calcJacInv => mesh_calcJacInv
     procedure :: calcJac => mesh_calcJac
     procedure :: calcGrad => mesh_calcGrad

  end type mesh

  ! type :: fields
  !
  ! Contains data about the fields to be evolved in the problem

  type fields
     ! BASIC PROPERTIES
     integer :: num_fields
     real(kind(0.0d0)), dimension(:,:) :: values 		!num_fields, num_pt
     type(mesh) :: base_mesh
     
     ! DERIVED PROPERTIES
     real(kind(0.0d0)), dimension(:,:,:) :: elem_values 	!num_elem, num_fields, num_quad
  contains
     ! PROCEDURES
     procedure :: setMesh => fields_setMesh
     procedure :: init => fields_init
  end type fields

  ! type :: fn_ptr
  !
  ! Contains a pointer to a function. Allows us to create functional arrays.

  type fn_ptr(dim, num_fields, num_quad)
     integer, len :: dim, num_fields, num_quad
     procedure(RHS), pointer :: fp

     abstract interface
        function RHS(values, diff, grad) result(ans)
           real(kind(0.0d0)), dimension(num_fields, num_quad), intent(in) :: u
           real(kind(0.0d0)), dimension(dim, num_quad, num_quad) :: diff, grad
           real(kind(0.0d0)), dimension(num_quad) :: ans
        end function RHS
     end interface
  end type fn_ptr

  ! type :: equation
  !
  ! Contains data about the equations governing the problem, as well as how we'd like to evolve them.

  type equation(dim, num_fields, num_quad)
    ! BASIC PROPERTIES
    integer, len :: dim, num_fields, num_quad
    type(fn_ptr(dim, num_fields, num_quad)), dimension(:,:), allocatable :: operators
    integer :: num_split
    double precision :: theta

  contains
    procedure :: addOps => equation_addOps
    procedure :: addInfo => equation_addInfo
  end type equation

  contains

!!!!!!!!!!!!!!!!!!!!
! MESH SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!

  ! subroutines mesh%getCoords and mesh%getIndex
  !
  ! Within an element, nodal points are parametrized by either the master index
  ! or the mesh%dim-tuplet of indices, each of which runs from 1 to (mesh%deg + 1).
  !
  ! getCoords and getIndex determine the relationship between these two indexing methods.

  subroutine mesh_setParam(this, dim, deg, num_elem, num_pt, num_lyr)
     class(mesh), intent(inout) :: this
     integer, intent(in) :: dim, deg, num_elem, num_pt, num_lyr
     integer :: num_quad
     num_quad = (deg+1)**dim

     this%dim = dim
     this%deg = deg
     this%num_elem = num_elem
     this%num_pt = num_pt
     this%num_lyr = num_lyr
     this%num_quad = num_quad

     allocate(this%structure(num_elem, num_quad))
     allocate(this%quad(deg+1))
     allocate(this%weights(deg+1))
     allocate(this%coords(num_pt,dim))
     allocate(this%weightProds(num_quad))
     allocate(this%derivs((deg+1),(deg+1)))
     allocate(this%diffs(dim,num_quad,num_quad))
     allocate(this%jac_mat(num_elem,dim,dim,num_quad))
     allocate(this%jac_inv(num_elem,dim,dim,num_quad))
     allocate(this%grad(num_elem,dim,num_quad,num_quad))
     allocate(this%jac(num_elem,num_quad))
  end subroutine mesh_param

  subroutine mesh_setForm(this, structure, coords, quad, weights)
     class(mesh), intent(inout) :: this
     integer, dimension(this%num_elem,this%num_quad), intent(in) :: structure
     real(kind(0.0d0)), dimension(this%num_pt,this%dim), intent(in) :: coords
     real(kind(0.0d0)), dimension(this%deg+1), intent(in) :: quad, weights

     this%structure = structure
     this%coords = coords
     this%quad = quad
     this%weights = weights
  end subroutine mesh_form

  subroutine mesh_getCoords(this, ind, indices)
     implicit none
     class(mesh), intent(in) :: this
     integer, intent(in) :: ind
     integer, dimension(this%dim), intent(out) :: indices
     integer :: i, r

     do i = this%dim, 1, -1
        indices(i) = ind/((this%deg+1)**(i-1)) + 1
        r = mod(ind, (this%deg+1)**(i-1))
        ind = r
     enddo
     indices(1) = indices(1) + 1
  end subroutine

  subroutine mesh_getIndex(this, indices, ind)
     implicit none
     class(mesh), intent(in) :: this
     integer, dimension(this%dim), intent(in) :: indices
     integer, intent(out) :: ind
     integer :: i

     ind = indices(1)
     do i = 1, (this%dim-1)
        ind = ind + (indices(i+1) - 1)*(this%deg + 1)**i
     enddo
  end subroutine

  ! subroutine mesh%calcDPhi
  !
  ! Calculates the derivatives of the 1D interpolation function at the 1D nodal points.

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

  ! subroutine mesh%calcDiff
  !
  ! Calculates the derivatives of the 3D interpolation functions at the 3D nodal points.
  !
  ! Derivatives calculated with respect to the reference cube.

  subroutine mesh_calcDiff(this)
     implicit none
     class(mesh), intent(inout) :: this
     integer, dimension(this%dim) :: iv, jv
     integer :: a, b, i, j
     integer :: n
     real(kind(0.0d0)) :: diff

     n = (this%deg+1)**this%dim

     do a = 1, this%dim
        do i = 1, n
           do j = 1, n
              call this%getCoords(i,iv)
              call this%getCoords(j,jv)
              diff = this%derivs(jv(a),iv(a))
              do b = 1, this%dim
                 if ((b .ne. a) .and. (jv(b) .ne. iv(b))) then
                    diff = 0.0d0
                 endif
              enddo
              this%diffs(a,i,j) = diff
            enddo
         enddo
      enddo
  end subroutine mesh_calcDiff

  ! subroutines mesh%calcJacMat, mesh%calcJacInv, mesh%calcJac
  !
  ! Calculates the Jacobian matrix detailing the map from the
  ! reference cube to the element on the mesh, as well as its
  ! inverse and determinant.

  subroutine mesh_calcJacMat(this)
     implicit none
     class(mesh), intent(inout) :: this
     integer :: m, a, b, i, j, n 
     real(kind(0.0d0)) :: jsum

     n = (this%deg+1)**this%dim

     do m = 1, this%num_elem
        do a = 1, this%dim
           do b = 1, this%dim
              do i = 1, n
                 jsum = 0.0d0
                 do j = 1, n
                   jsum = jsum + this%diffs(a,i,j)*this%coords(this%structure(m,j),b)
                 enddo
                 this%jac_mat(m,a,b,i) = jsum
              enddo
           enddo
        enddo
     enddo
  end subroutine mesh_calcJacMat

  subroutine mesh_calcJacInv(this)
     implicit none
     class(mesh), intent(inout) :: this

     do m = 1, this%num_elem
        do i = 1, this%num_quad
           call mat_inv(this%dim,this%jac_mat(m,:,:,i),this%jac_inv(m,:,:,i))
        enddo
     enddo
  end subroutine mesh_calcJacInv

  subroutine mesh_calcJac(this)
     implicit none
     class(mesh), intent(inout) :: this

     do m = 1, this%num_elem
        do i = 1, this%num_quad
           this%jac(m,i) = det(this%dim,this%jac_mat(m,:,:,i))
        enddo
     enddo
  end subroutine mesh_calcJac

  ! subroutine mesh_calcGrad
  ! 
  ! Calculates the gradient matrix, which can be used to get derivatives
  ! of functions with respect to the mesh coordinates.

  subroutine mesh_calcGrad(this)
     implicit none
     class(mesh), intent(inout) :: this
     real(kind(0.0d0)) :: dsum

     n = (this%deg+1)**this%dim
     
     do i = 1, n
        do j = 1, n
           do a = 1, dim
              dsum = 0.0d0
              do b = 1, dim
                 dsum = dsum + this%jac_inv(m,a,b,i)*this%diffs(b,i,j)
              enddo
              this%grad(m,a,i,j) = dsum
           enddo
        enddo
     enddo
  end subroutine mesh_calcGrad

  subroutine mesh_calcWeights(this)
     implicit none
     class(mesh), intent(inout) :: this
     real(kind(0.0d0)) :: prod
     integer, dimension(this%dim) :: indices
     integer :: i, a

     n = (this%deg+1)**this%dim

     do i = 1, n
        call this%getCoords(i, indices)
        prod = 1.0d0
        do a = 1, this%dim
           prod = prod*weights(indices(a))
        enddo
        this%weight_prods(i) = prod
     enddo
  end subroutine mesh_calcWeights

  subroutine mesh_calcProj(this)
     implicit none
     class(mesh), intent(inout) :: this
     real(kind(0.0d0)) :: psum
     real(kind(0.0d0)), dimension(this%num_pt) :: W
     integer :: p, m, n, i


     n = (this%deg+1)**this%dim

     do p = 1, this%num_pt
        W(p) = 0.0d0
     enddo

     do m = 1, this%num_elem
        do i = 1, n
           psum = W(this%structure(m,i))
           W(this%structure(m,i)) = psum + this%weight_prods(i)*this%jac(m,i)
        enddo
     enddo

     do m = 1, this%num_elem
        do i = 1, n
           projector(m,i) = this%weight_prods(i)*this%jac(m,i)/W(this%structure(m,i))
        enddo
     enddo
  end subroutine mesh_calcProj

!!!!!!!!!!!!!!!!!!!!!!
! FIELDS SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!

  subroutine fields_setMesh(this, base_mesh, num_fields)
     class(fields), intent(inout) :: this
     type(mesh), intent(in) :: base_mesh

     this%base_mesh = base_mesh

     allocate(this%values(num_fields,base_mesh%num_pt))
     allocate(this%elem_values(num_fields, base_mesh%num_elem, base_mesh%num_pt))
  end subroutine

  subroutine fields_init(this, values)
     class(fields), intent(inout) :: this
     real(kind(0.0d0)), dimension(this%num_fields,this%base_mesh%num_pt) :: values

     this%values = values
  end subroutine

  subroutine fields_getElemValues(this)
     implicit none
     class(fields), intent(inout) :: this
     integer :: p, i, f

     do m = 1, this%num_elem
        do f = 1, this%num_fields
           do i = 1, this%num_ept
              p = this%base_mesh%structure(m,i)
              this%elem_values(f, m, i) = this%values(f, p)
           enddo
        enddo
     enddo
  end subroutine fields_getElemValues

  subroutine fields_projectElemValues(this)
     implicit none
     class(fields), intent(inout) :: this
     integer :: n, p, f, m, i

     n = (this%deg+1)**this%dim     


     do p = 1, this%num_pt
        values(p) = 0.0d0
     enddo

     do f = 1, this%num_fields
        do m = 1, this%num_elem
           do i = 1, n
              psum = values(f,this%base_mesh%structure(m,i))
              values(f,this%base_mesh%structure(m,i)) = psum + projector(m,i)*elem_values(m,f,i)
           enddo
        enddo

        do m = 1, this%num_elem
           do i = 1, n
              p = this%base_mesh%structure(m,i)
              this%elem_values(m, f, i) = this%values(f, p)
           enddo
        enddo
     enddo
  end subroutine fields_projectElemValues

!!!!!!!!!!!!!!!!!!!!!!!!
! EQUATION SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equation_addOps(this, num_split, operators)
     class(equation(:,:,:)) :: this
     integer :: num_split
     real(kind(0.0d0)), dimension(this%num_fields, num_split) :: operators

     this%num_split = num_split
     allocate(this%operators(this%num_fields, num_split))
     this%operators = operators
  end subroutine equation_addOps

  subroutine equation_addTheta(this, theta)
     class(equation(:,:,:)) :: this
     real(kind(0.0d0)) :: theta

     this%theta = theta
  end subroutine equation_addInfo
end module hesiod

