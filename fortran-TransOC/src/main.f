
	program TransOC
	use modmain
	use basisstates, only: mkbasis
	use hamiltonian, only: makeHg,MakeHgMulti
	use hoppings, only: dhops
	implicit none
	integer i,nnz,j,n1,n2,k,ih,is,ib,itype1,itype2,nt
	real(kind=4), allocatable,dimension(:,:):: mat,matf
	integer(kind=4), allocatable,dimension(:):: row,col
	
	! to test make-hg-multi
	integer(kind=4), dimension(3):: itlist


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
	

	!wj=10; wc=1;
	! update mapb and mapt
	!call UpdateMapB(wj,wc,mapb,notusedb,nnub)
	!call UpdateMapT(wj,wc,mapt,notusedt,nnut)





	itlist = (/ 1,2,3 /);
	nt = 3;
	if (1==1) then
			call MakeHgMulti(itlist,nt)
			write(*,*) "main: done..."
			write(*,*) " >>>> MakeHgMulti(itlist,nt) for itlist=",itlist
	end if


	if (1==1) then
		do j=1,13
			call makeHg(j,nnz)
			write(*,*) " itype, nnz= ",j, nnz
		end do
	end if


	if (1==0) then
		call dhops()
		write(*,*) " dhops done....  "
	end if
	


	if (1==1) then
		call DPhiAn1()
		write(*,*) " DPhiAn1() done.... "
	end if

	if (1==0) then
		call DPhiAn2()	
		write(*,*) " DPhiAn2() done.... "
	end if


!============================================
! dphi annihilation: ih = 5; multiply test 
!============================================
	IF (1==1) then
	ih = 5; is=1;	
	itype1 = 1; itype2= 4 
	
	n1 = Hg(itype1)%ntot
	n2 = Hg(itype2)%ntot
	! Ht(n1 x n2) . mat( n2 x n2 ) ==> matf(n1 x n2 )
	allocate(mat(n2,n2))
	allocate(matf(n1,n2))
	call random_number(mat)

	nnz = hop(ih)%ht(1,is)%nnz;
	write(*,*) "Ht chal 1,2 nnz: ",nnz
	allocate(row(nnz))
	allocate(col(nnz))

	row = (/ (i, i=1,n1) /)
	col = hop(ih)%ht(1,is)%col

	write(*,*) "nnz,n1,n2 = ",nnz,n1,n2
	call multiply(row,col,nnz,mat,n2,n2,matf,n1,n2)
	write(*,*) " multiply done....  "
	endif
!============================================



	allocate(hop(25)%ht(1,1))
	call LossKappa()


	allocate(hop(26)%ht(1,1))
	call LossGamma(1)

	write(*,*) " kappa, gamma done .... "




	
	
	end program
