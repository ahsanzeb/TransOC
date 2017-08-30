
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
	real(kind=4), dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	real(kind=4), dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i

	matf = 0
	!write(*,*) "col : ", col
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
	implicit none
	integer, intent(in) :: nnz,n1,n2,n3,n4
	integer, dimension(nnz), intent(in) :: row ! Ht map
	real(kind=4), dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	real(kind=4), dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i,j

	matf = 0
	do i=1,nnz ! diagonal but not full diagonal, selected rows/cols
		j = row(i);
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
	! (3 - matrix mati is dense)
	!			Eigenvectors matrix Uf is dense
	implicit none
	integer, intent(in) :: nnz,n1,n2,n3,n4
	integer, dimension(nnz), intent(in) :: col ! Ht map
	real(kind=4), dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	real(kind=4), dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i

	matf = 0
	do i=1,nnz ! nnz = dim of initial Hilbert space
		matf(i,:) = mat(col(i),:)
	end do

	end subroutine multiplydc
!-----------------------------------	
