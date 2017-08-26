
	program TransOC
	use modmain
	use basisstates, only: mkbasis
	use hamiltonian, only: makeHg
	use hoppings, only: dhops
	implicit none
	integer i,nnz,j,n1,n2,k
	real(kind=4), allocatable,dimension(:,:):: mat,matf
	integer(kind=4), allocatable,dimension(:):: row,col

	! set no of active sites and excitations
	na = 10; nx = 5;

	crosshops = .true.;
	detuning = .true.;
	call mkbasis(na,nx)
	write(*,*) " basis done....  "

	!do k=1,m1,1	
	!	ntot = basis(i)%pntr(j+2) - basis(i)%pntr(j+1); ! for config of this type
	!	if(allocated(sets))deallocate(sets)
	!	allocate(sets(ntot,k))
	!	call mksets(n,k,ntot,sets)
	





	if (1==1) then
		do j=1,5
			call makeHg(j,nnz)
			write(*,*) " itype, nnz= ",j, nnz
		end do
	end if


	if (1==1) then
		call dhops()
	end if

	
	write(*,*) " dhops done....  "

	n1 = Hg(1)%ntot
	n2 = Hg(2)%ntot
	! Ht(n1 x n2) . mat( n2 x n2 ) ==> matf(n1 x n2 )
	allocate(mat(n2,n2))
	allocate(matf(n1,n2))
	call random_number(mat)

	nnz = hop(1)%ht(1,1)%nnz;
	write(*,*) "Ht chal 1,3 nnz: ",nnz
	allocate(row(nnz))
	allocate(col(nnz))

	row = hop(1)%ht(1,1)%row
	col = hop(1)%ht(3,1)%col

	write(*,*) "nnz,n1,n2 = ",nnz,n1,n2
  !	hop 1, chan 3, is=1 => site 1
  ! Ht(n1 x n2) . mat( n2 x n2 ) ==> matf(n1 x n2 )
	write(*,*) "max: row, col= ",maxval(row),maxval(col)
	call multiply(row,col,nnz,mat,n2,n2,matf,n1,n2)



	nnz = hop(1)%ht(4,1)%nnz
	write(*,*) "Ht chan 2,4 nnz: ",nnz

	write(*,*) " multiply done....  "
	







	
	
	end program
