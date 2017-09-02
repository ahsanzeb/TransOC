
	program TransOC
	use modmain
	use init, only: initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: makeHg,MakeHgMulti
	use Hoppings, only: AllHops ! 
	use rates, only: CalRates
	use selection, only: ihSelect, icsSelect
	use modways, only: UpdateWays, UpdateOcc
	use readinput, only: input
	implicit none
	integer i,nnz,j,n1,n2,k,ih,is,ib,itype1,itype2,nt
	real(kind=4), allocatable,dimension(:,:):: mat,matf
	integer(kind=4), allocatable,dimension(:):: row,col
	integer :: wj,wc !,itype
	! to test make-hg-multi
	integer(kind=4), dimension(13):: itlist



	crosshops = .true.;
	detuning = .true.;

	!	initialise maps etc
	call initialise()

	! set no of active sites and excitations
	na = 10; nx = 5;

	call mkbasis(na,nx)
	write(*,*) " basis done....  "





	call AllHops()


	!do k=1,m1,1	
	!	ntot = basis(i)%pntr(j+2) - basis(i)%pntr(j+1); ! for config of this type
	!	if(allocated(sets))deallocate(sets)
	!	allocate(sets(ntot,k))
	!	call mksets(n,k,ntot,sets)
	


	wj=3; wc=4;
	write(*,*) "--------wj,wc=3,4----------"
	write(*,*) mapb%map
	write(*,*) mapt%map
	
	! update mapb and mapt
	call UpdateMapB(wj,wc)
	call UpdateMapT(wj,wc)
	call UpdateGroupTB()
	write(*,*) "----- updated mapb,mapt--------"
	write(*,*) mapb%map
	write(*,*) mapt%map
	write(*,*) "----- updated grouptb --------"
	write(*,*) "ntb = ",mapt%ntb
	do i=1,5
		write(*,*) mapt%grouptb(i,:)
	end do




	wj=5; wc=3;
	write(*,*) "------ wj,wc=5,4------------"
	write(*,*) mapb%map
	write(*,*) mapt%map
	
	! update mapb and mapt
	call UpdateMapB(wj,wc)
	call UpdateMapT(wj,wc)
	call UpdateGroupTB()
	write(*,*) "----- updated mapb,mapt--------"
	write(*,*) mapb%map
	write(*,*) mapt%map
	write(*,*) "----- updated grouptb --------"
	write(*,*) "ntb = ",mapt%ntb
	do i=1,5
		write(*,*) mapt%grouptb(i,:)
	end do


	! update basis
	itype = itypes(wj,wc);
	na = na + dna(itype);
	nx = na + dnx(itype);
	call mkbasis(na,nx)
	write(*,*) " basis done....  "

	
	
	! update hamiltonian
	if (1==1) then
		write(*,*) "main: ntb = ",mapt%ntb
		do ib=1,5
			nt = mapt%ntb(ib);
			if(nt > 0) then
				itlist(1:nt) = mapt%grouptb(ib,1:nt)
				call MakeHgMulti(itlist(1:nt),nt)
				write(*,*) "main: done..."
				write(*,*) " >>>> MakeHgMulti for itlist=",itlist(1:nt)
			end if
		end do
	end if

	stop



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

	!nnz = hop(ih)%ht(1,is)%nnz;
	write(*,*) "Ht chal 1,2 nnz: ",nnz
	allocate(row(nnz))
	allocate(col(nnz))

	row = (/ (i, i=1,n1) /)
	!col = hop(ih)%ht(1,is)%col

	write(*,*) "nnz,n1,n2 = ",nnz,n1,n2
	!call multiply(row,col,nnz,mat,n2,n2,matf,n1,n2)
	write(*,*) " multiply done....  "
	endif
!============================================



	!allocate(hop(25)%ht(1,1))
	!call LossKappa()


	!allocate(hop(26)%ht(1,1))
	!call LossGamma(1)

	write(*,*) " kappa, gamma done .... "




	
	
	end program
