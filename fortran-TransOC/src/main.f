
	program TransOC
	use modmain
	use basisstates, only: mkbasis
	use hamiltonian, only: makeHg
	!use lists
	implicit none
	integer i,nnz,j


	! set no of active sites and excitations
	na = 10; nx = 5;


	crosshops = .true.;
	detuning = .true.;
	call mkbasis(na,nx)

	do j=1,5
		call makeHg(j,nnz)
		write(*,*) " itype, nnz= ",j, nnz
	end do

	
	end program
