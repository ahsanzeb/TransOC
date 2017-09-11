
	module modmultiply
	implicit none

	contains

	subroutine multiply(row,col,nnz,mat,n1,n2,matf,n3,n4)
	!use kinds, only: EvecKind
	! use only when:
	!	1 - the sparse matrix has at most a single element in a row
	!			I think all Ht's have this property
	!	2 - the nonzero elements of the sparse matrix are all 1
	!			I will multiply t's later to get amplitudes, so they are 1 for Ht's
	! (3 - matrix mati is dense)
	!			Eigenvectors matrix Uf is dense
	implicit none
	integer, intent(in) :: nnz,n1,n2,n3,n4
	integer, dimension(nnz), intent(in) :: row ! Ht map
	integer, dimension(nnz), intent(in) :: col ! Ht map
	double precision, dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	double precision, dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i

	!write(*,*)"nnz,n1,n2,n3,n4 = ",nnz,n1,n2,n3,n4
	!write(*,*)minval(row), minval(col)
	!write(*,*) col

	matf = 0.d0
	do i=1,nnz
		matf(row(i),:) = mat(col(i),:)
	end do
	end subroutine multiply
!-----------------------------------	
	subroutine multiplyd(row,nnz,mat,n1,n2,matf,n3,n4)
	!use kinds, only: EvecKind
	! use only when:
	!	0 - the sparse matrix is diagonal
	!	1 - the sparse matrix has at most a single element in a row
	!			I think all Ht's have this property
	!	2 - the nonzero elements of the sparse matrix are all 1
	!			I will multiply t's later to get amplitudes, so they are 1 for Ht's
	! (3 - matrix mati is dense)
	!			Eigenvectors matrix Uf is dense
	! in: row contains row(=col) indexes for nonzero elements
	implicit none
	integer, intent(in) :: nnz,n1,n2,n3,n4
	integer, dimension(nnz), intent(in) :: row ! Ht map
	double precision, dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	double precision, dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i,j

	matf = 0.0d0
	do i=1,nnz ! diagonal but not full diagonal, selected rows/cols
		j = row(i)
		matf(j,:) = mat(j,:)
	end do

	end subroutine multiplyd
!-----------------------------------	



	subroutine multiplydc(col,nnz,mat,n1,n2,matf,n3,n4)
	!use kinds, only: EvecKind
	! use only when:			
	!	1 - the sparse matrix has exactly a single element in each row
	!			D,Phi annihilation is an example. 
	!			nnz = dim of initial Hilbert space
	!	2 - the nonzero elements of the sparse matrix are all 1
	!			I will multiply t's later to get amplitudes, so they are 1 for Ht's
	! (3 - matrix mat is dense)
	!			Eigenvectors matrix Uf is dense
	! in: col contains col indexes for rows=1,nnz; (nnz=n1)
	implicit none
	integer, intent(in) :: nnz,n1,n2,n3,n4
	integer, dimension(nnz), intent(in) :: col ! Ht map
	double precision, dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	double precision, dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i

	matf = 0.0d0
	do i=1,nnz ! nnz = dim of initial Hilbert space
		matf(i,:) = mat(col(i),:)
	end do

	!write(*,*) "mult: n1,n2, n3,n4", n1,n2, n3,n4
	!write(*,*) "mult: nnz, shape(col)",nnz, shape(col)
	!write(*,*) "mult: col",col

	!write(*,*) "mult: mat = "
	!write(*,*) mat
	!write(*,*) "mult: matf = "
	!write(*,*) matf

	end subroutine multiplydc
!-----------------------------------	

	end module modmultiply
