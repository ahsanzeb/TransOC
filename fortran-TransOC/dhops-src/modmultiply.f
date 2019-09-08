
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
	subroutine multiplydz(row,nnz,mat,n1,n2,matf,n3,n4)
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
	double complex, dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	double complex, dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i,j

	matf = 0.0d0
	do i=1,nnz ! diagonal but not full diagonal, selected rows/cols
		j = row(i)
		matf(j,:) = mat(j,:)
	end do

	end subroutine multiplydz
!-----------------------------------	


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
	subroutine multiplydk(row,nnz,mat,n1,n2,matf,n3,n4)
	! for cavity losses: includes the sqrt prefactors
	! use only when:
	!	0 - the sparse matrix is diagonal
	!	1 - the sparse matrix has at most a single element in a row
	!	2 - the nonzero elements of the sparse matrix are all *?
	! (3 - matrix mati is dense)
	!			Eigenvectors matrix Uf is dense
	! in: row contains row(=col) indexes for nonzero elements
	implicit none
	integer, intent(in) :: nnz,n1,n2,n3,n4
	double precision, dimension(nnz), intent(in) :: row ! Ht map
	double precision, dimension(n1,n2), intent(in) :: mat  ! Uf eigenvectors
	double precision, dimension(n3,n4), intent(out) :: matf! Ht.Uf output

	! local
	integer :: i

	matf = 0.0d0
	do i=1,nnz ! diagonal but not full diagonal, first nnz elements
		matf(i,:) = mat(i,:) * row(i)
	end do
	! since row/first index moves continuously,
	!efficient if this is inner loop, ??? below: 
	!do i=1,n2
	!	matf(1:nnz,i) = mat(1:nnz,i) * row(:) ! elem by elem multipl
	!enddo
	end subroutine multiplydk
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
!  	for rho_init to rho_final:
! 		multiply2, multiplyd2, multiplydc2
! 		two steps:
! 			1. ir-th row of rho_init to ic-th row of aux
!			2. ic-th col of aux to ir-th col of rho_final;
!		ir,ic are row and col indices of 
!		non-zero (=1) elements of Ht
!-----------------------------------	
	subroutine multiply2(row,col,nnz,mat,n1,n2,matf)
	! Ht's in the format of multiply, row and col indices
	implicit none
	integer, intent(in) :: nnz,n1,n2
	integer, dimension(nnz), intent(in) :: row ! Ht map, row indices
	integer, dimension(nnz), intent(in) :: col ! Ht map, col indices
	double precision, dimension(n1,n1), intent(in) :: mat ! rho_init
	double precision, dimension(n2,n2), intent(out) :: matf ! Ht^T.rho_init.Ht

	! local
	double precision, dimension(n2,n1):: aux! Ht.Uf
	integer :: i
	matf = 0.d0;
	aux = 0.0d0;
	! aux = HtT.rho; HtT = Transpose(Ht)
	! ir-th row to ic-th row
	do i=1,nnz
		aux(col(i),:) = mat(row(i),:)
	end do
	!	aux.Ht; 
	! ic-th col to ir-th col; (ic,ir map of Ht)
	do i=1,nnz
		matf(:,row(i)) = aux(:,col(i))
	end do
	
	end subroutine multiply2
!-----------------------------------	
	subroutine multiplyd2(row,nnz,mat,n1,n2,matf)
	! Ht format given in multiplyd()
	implicit none
	integer, intent(in) :: nnz,n1,n2
	integer, dimension(nnz), intent(in) :: row ! Ht map
	double precision, dimension(n1,n1), intent(in) :: mat  !
	double precision, dimension(n2,n2), intent(out) :: matf! 

	! local
	integer :: i,j
	double precision, dimension(n2,n1):: aux

	aux = 0.0d0
	do i=1,nnz
		j = row(i)
		aux(j,:) = mat(j,:)
	end do

	matf = 0.0d0
	do i=1,nnz 
		j = row(i)
		matf(:,j) = aux(:,j)
	end do
	
	end subroutine multiplyd2
!-----------------------------------	
	subroutine multiplydc2(col,nnz,mat,n1,n2,matf)
	! Ht format given in multiplydc()
	implicit none
	integer, intent(in) :: nnz,n1,n2
	integer, dimension(nnz), intent(in) :: col ! Ht map
	double precision, dimension(n1,n1), intent(in) :: mat 
	double precision, dimension(n2,n2), intent(out) :: matf
	! local
	integer :: i
	double precision, dimension(n2,n1):: aux

	aux = 0.0d0
	do i=1,nnz
		aux(col(i),:) = mat(i,:)
	end do
	
	matf = 0.0d0
	do i=1,nnz
		matf(:,i) = aux(:,col(i))
	end do

	end subroutine multiplydc2
!-----------------------------------	



	end module modmultiply
