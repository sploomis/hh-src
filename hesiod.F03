!--------------------*-------------------------------------!
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
  use ovid
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
     real(kind(0.0d0)), dimension(:,:,:), allocatable :: diff			!dim, num_quad, num_quad
     real(kind(0.0d0)), dimension(:,:,:,:), allocatable :: jac_mat, jac_inv 		!num_elem, dim, dim, num_quad
     real(kind(0.0d0)), dimension(:,:,:,:), allocatable :: grad 			!num_elem, dim, num_quad, num_quad
     real(kind(0.0d0)), dimension(:,:), allocatable :: jac, Jw				!num_elem, num_quad
     real(kind(0.0d0)), dimension(:,:), allocatable :: projector			!num_elem, num_quad

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
     procedure :: calcWeights => mesh_calcWeights
     procedure :: calcProj => mesh_calcProj
     procedure :: load => mesh_load

  end type mesh

  ! type :: fields
  !
  ! Contains data about the fields to be evolved in the problem

  type fields
     ! BASIC PROPERTIES
     integer :: num_fields
     real(kind(0.0d0)), dimension(:,:), allocatable :: values 		!num_fields, num_pt
     type(mesh) :: base
     
     ! DERIVED PROPERTIES
     real(kind(0.0d0)), dimension(:,:,:), allocatable :: elem_values 	!num_elem, num_fields, num_quad
  contains
     ! PROCEDURES
     procedure :: setMesh => fields_setMesh
     procedure :: init => fields_init
     procedure :: getElemValues => fields_getElemValues
     procedure :: projectElemValues => fields_projectElemValues
  end type fields

  ! type :: fn_ptr
  !
  ! Contains a pointer to a function. Allows us to create functional arrays.

  interface rhs_form
     function RHS(dim, num_fields, num_quad, diff,grad,Jw, u, t) result(ans)
        integer, intent(in) :: dim, num_fields, num_quad
        real(kind(0.0d0)), dimension(num_fields, num_quad) :: u
        real(kind(0.0d0)), dimension(dim, num_quad, num_quad) :: diff, grad
        real(kind(0.0d0)), dimension(num_quad) :: Jw
        real(kind(0.0d0)), dimension(num_quad) :: ans
        real(kind(0.0d0)) :: t
     end function RHS
  end interface rhs_form

  type fn_ptr
     procedure(RHS), pointer, nopass :: fp
  end type fn_ptr

  ! type :: equation
  !
  ! Contains data about the equations governing the problem, as well as how we'd like to evolve them.

  type equation
    ! BASIC PROPERTIES
    integer :: dim, num_fields, num_quad
    type(fn_ptr), dimension(:,:), allocatable :: operators
    integer :: num_split
    double precision :: theta

  contains
    procedure :: addOps => equation_addOps
    procedure :: addTheta => equation_addTheta

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
     allocate(this%diff(dim,num_quad,num_quad))
     allocate(this%jac_mat(num_elem,dim,dim,num_quad))
     allocate(this%jac_inv(num_elem,dim,dim,num_quad))
     allocate(this%grad(num_elem,dim,num_quad,num_quad))
     allocate(this%jac(num_elem,num_quad))
     allocate(this%jw(num_elem,num_quad))
     allocate(this%projector(num_elem,num_quad))
  end subroutine mesh_setParam

  subroutine mesh_setForm(this, structure, coords, quad, weights)
     class(mesh), intent(inout) :: this
     integer, dimension(this%num_elem,this%num_quad), intent(in) :: structure
     real(kind(0.0d0)), dimension(this%num_pt,this%dim), intent(in) :: coords
     real(kind(0.0d0)), dimension(this%deg+1), intent(in) :: quad, weights

     this%structure = structure
     this%coords = coords
     this%quad = quad
     this%weights = weights
  end subroutine mesh_setForm

  subroutine mesh_getCoords(this, ind, indices)
     implicit none
     class(mesh), intent(in) :: this
     integer, intent(in) :: ind
     integer, dimension(this%dim), intent(out) :: indices
     integer :: ind1, i, r

     ind1 = ind

     do i = this%dim, 1, -1
        indices(i) = ind1/((this%deg+1)**(i-1))
        r = mod(ind, (this%deg+1)**(i-1))
        ind1 = r
     enddo
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
           if (i .eq. j) then
              psum = 0.0d0
              do k = 1, this%deg+1
                 if (k .ne. i) then
                    psum = psum + 1.0d0/(this%quad(i) - this%quad(k))
                 endif
              enddo
              this%derivs(i,j) = psum
           else
              prod = 1.0d0/(this%quad(i) - this%quad(j))
              do k = 1, this%deg+1
                 if ((k .ne. j) .and. (k .ne. i)) then
                    prod = prod*(this%quad(j) - this%quad(k))/(this%quad(i) - this%quad(k))
                 endif
              enddo
              this%derivs(i,j) = prod
           endif
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
                 if (b .ne. a) then
                    if (jv(b) .ne. iv(b)) diff = 0.0d0
                 endif
              enddo
              this%diff(a,i,j) = diff
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
                   jsum = jsum + this%diff(a,i,j)*this%coords(this%structure(m,j),b)
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
     integer :: m, i

     do m = 1, this%num_elem
        do i = 1, this%num_quad
           call mat_inv(this%dim,this%jac_mat(m,:,:,i),this%jac_inv(m,:,:,i))
        enddo
     enddo
  end subroutine mesh_calcJacInv

  subroutine mesh_calcJac(this)
     implicit none
     class(mesh), intent(inout) :: this
     integer :: m, i

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
     integer :: m, n, i, j, a, b

     n = (this%deg+1)**this%dim
     do m = 1, this%num_elem
        do i = 1, n
           do j = 1, n
              do a = 1, this%dim
                 dsum = 0.0d0
                 do b = 1, this%dim
                    dsum = dsum + this%jac_inv(m,a,b,i)*this%diff(b,i,j)
                 enddo
                 this%grad(m,a,i,j) = dsum
              enddo
           enddo
        enddo
     enddo
  end subroutine mesh_calcGrad

  ! subroutines mesh%calcWeights and mesh%calcProj
  !
  ! Calculates the weights at the quadrature points and the resulting mass matrix.
  ! These values are then used in calculating the projection matrix.

  subroutine mesh_calcWeights(this)
     implicit none
     class(mesh), intent(inout) :: this
     real(kind(0.0d0)) :: prod
     integer, dimension(this%dim) :: indices
     integer :: n, i, a

     n = (this%deg+1)**this%dim

     do i = 1, n
        call this%getCoords(i, indices)
        prod = 1.0d0
        do a = 1, this%dim
           prod = prod*this%weights(indices(a))
        enddo
        this%weightProds(i) = prod
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
           W(this%structure(m,i)) = psum + this%weightProds(i)*this%jac(m,i)
        enddo
     enddo

     do m = 1, this%num_elem
        do i = 1, n
           this%projector(m,i) = this%weightProds(i)*this%jac(m,i)/W(this%structure(m,i))
           this%Jw(m,i) = this%weightProds(i)*this%jac(m,i)
        enddo
     enddo
  end subroutine mesh_calcProj

  subroutine mesh_load(this, mshfile)
     class(mesh), intent(inout) :: this
     character(len=*) :: mshfile
     integer :: dim, deg, n, m, nquad, en, em, i, a
     real(kind(0.0d0)), dimension(:), allocatable :: quad, weights
     integer, dimension(:,:), allocatable :: structure
     real(kind(0.0d0)), dimension(:,:), allocatable :: coords
        
     open(unit=111,file=trim(mshfile))
     read(111,*) dim, deg, em, en
     nquad = (deg+1)**dim
     
     allocate(quad(deg+1))
     allocate(weights(deg+1))
     allocate(structure(em,nquad))
     allocate(coords(en,dim))
     
     read(111,*) (quad(i),i=1,deg+1)
     read(111,*) (weights(i),i=1,deg+1)
     
     do m = 1, em
        read(111,*) (structure(m,i),i=1,nquad)
     enddo 

     do n = 1, en
        read(111,*) (coords(n,a),a=1,dim)
     enddo

     call this%setParam(dim, deg, em, en, 1)
     call this%setForm(structure, coords, quad, weights)
     call this%calcDPhi
     call this%calcDiff
     call this%calcJacMat
     call this%calcJacInv
     call this%calcJac
     call this%calcGrad
     call this%calcWeights
     call this%calcProj
     close(111)
  end subroutine


!!!!!!!!!!!!!!!!!!!!!!
! FIELDS SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!

  ! subroutines fields%setMesh and fields%init
  ! 
  ! Assigns the base mesh of the field and the initial values

  subroutine fields_setMesh(this, base, num_fields)
     class(fields), intent(inout) :: this
     type(mesh), intent(in) :: base
     integer, intent(in) :: num_fields

     this%base = base

     allocate(this%values(num_fields,base%num_pt))
     allocate(this%elem_values(num_fields, base%num_elem, base%num_pt))
  end subroutine

  subroutine fields_init(this, values)
     class(fields), intent(inout) :: this
     real(kind(0.0d0)), dimension(this%num_fields,this%base%num_pt) :: values

     this%values = values
  end subroutine

  ! subroutines fields%getElemValues and fields%projectElemValues
  !
  ! Converts between the nodal and elemental representation of the fields
  ! and enforces continuity on the elemental representation.

  subroutine fields_getElemValues(this)
     implicit none
     class(fields), intent(inout) :: this
     integer :: m, p, i, f

     do m = 1, this%base%num_elem
        do f = 1, this%num_fields
           do i = 1, this%base%num_quad
              p = this%base%structure(m,i)
              this%elem_values(f, m, i) = this%values(f, p)
           enddo
        enddo
     enddo
  end subroutine fields_getElemValues

  subroutine fields_projectElemValues(this)
     implicit none
     class(fields), intent(inout) :: this
     integer :: n, p, f, m, i
     real(kind(0.0d0)) :: psum

     n = (this%base%deg+1)**this%base%dim     

     do f = 1, this%num_fields
        do p = 1, this%base%num_pt
           this%values(f,p) = 0.0d0
        enddo
     enddo

     do f = 1, this%num_fields
        do m = 1, this%base%num_elem
           do i = 1, n
              psum = this%values(f,this%base%structure(m,i))
              this%values(f,this%base%structure(m,i)) = psum + this%base%projector(m,i)*this%elem_values(m,f,i)
           enddo
        enddo

        do m = 1, this%base%num_elem
           do i = 1, n
              p = this%base%structure(m,i)
              this%elem_values(m, f, i) = this%values(f, p)
           enddo
        enddo
     enddo
  end subroutine fields_projectElemValues

!!!!!!!!!!!!!!!!!!!!!!!!
! EQUATION SUBROUTINES !
!!!!!!!!!!!!!!!!!!!!!!!!

  ! subroutines equation%addOps and equation%addTheta
  !
  ! Lets the user define the equation, splitting order, and degree of implicity.

  subroutine equation_addOps(this, num_split, operators, dim, num_fields, num_quad)
     integer :: dim, num_fields, num_quad
     class(equation) :: this
     integer :: num_split
     type(fn_ptr), dimension(num_fields, num_split) :: operators

     this%dim = dim
     this%num_fields = num_fields
     this%num_quad = num_quad
     this%num_split = num_split
     allocate(this%operators(this%num_fields, num_split))
     this%operators = operators
  end subroutine equation_addOps

  subroutine equation_addTheta(this, theta)
     class(equation) :: this
     real(kind(0.0d0)) :: theta

     this%theta = theta
  end subroutine equation_addTheta
end module hesiod

