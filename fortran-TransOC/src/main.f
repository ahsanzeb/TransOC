
	program TransOC
	use modmain
	use init, only: initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: mkHamilt,DegenSectors
	use Hoppings, only: AllHops ! 
	use rates, only: CalRates
	use selection, only: ihSelect, icsSelect
	use modways, only: UpdateWays, UpdateOcc
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt=0
	integer:: nev,ncv

	write(*,*) "transoc: in amplitudes: test HtUf = 1.0d0;"

	! start message
	write(*,*) "transoc: started.... "

	! readinput file
	call input()
	
	!	initialise maps etc
	! set no of active sites and excitations
	call initialise()

	! main loop over number of hops asked
	do iter=1,niter

		write(*,*) "======================================="
		write(*,*) "          iter = ",iter
		write(*,*) "          N, m = ",na,nx
		
		! make basis states ksub
		call mkbasis(na,nx)
		write(*,*) " basis done....  "

		! make hamiltonian
		call mkHamilt()


		ntot = Hg(1)%ntot
		write(*,*) "ntot = ",ntot ! 176?

		nev = 10;
		ncv = 25; ! make it bigger??? not bigger than ntot?!!!
		call diagonalise(1, ntot, nev, ncv)

		stop

		
		!-----------------------------------
		! diagonalise
		! choose Psi if first iteration
		!-----------------------------------
		! set dummy eig for testing... 
		do i=1,13
		if (Hg(i)%xst) then
		ntot = Hg(i)%ntot
		eig(i)%ntot=ntot
		eig(i)%n1=ntot
		eig(i)%n2=ntot
		if (allocated(eig(i)%evec)) deallocate(eig(i)%evec)
		if (allocated(eig(i)%eval)) deallocate(eig(i)%eval)
		allocate(eig(i)%evec(ntot,ntot))
		allocate(eig(i)%eval(ntot))
		!eig(i)%evec = 0.0d0;
		!eig(i)%eval = 0.0d0;
		call random_number(eig(i)%evec)
		eig(i)%eval = (/ (j*1.0d-5, j=1,ntot) /) ! ordered

		! degenerate sectors
		if(allocated(eig(i)%esec))deallocate(eig(i)%esec)
		if(allocated(eig(i)%ind))deallocate(eig(i)%ind)
		allocate(eig(i)%esec(ntot))
		allocate(eig(i)%ind(ntot+1))
		call DegenSectors(eig(i)%eval,ntot,
     .   eig(i)%nsec,eig(i)%esec,eig(i)%ind) ! make degenerate sectors

		else
				write(*,*)"main: itype, ntot = ",i,ntot	
		endif ! ntot > 0

		enddo

		
		if(allocated(psi)) deallocate(psi)
		if(iter == 1 .or. fixmap) then
			allocate(psi(1,eig(1)%ntot))
		else
			write(*,*) "main: itype =",itype
			allocate(psi(1,eig(mapt%map(1))%ntot))
		endif
		!psi(eig(1)%ntot) = 0.0d0
		!call random_number(psi)
		psi = 1.0d0
		Einit = 0.0d0
		!-----------------------------------

		write(*,*)"sys%occ = ",sys%occ



		!allocate qt (req maph ==> ia )


		if(iter==1) then
			do ia=1,14
				allocate(qt(ia)%cs(4,1)) ! alloc/deallocate at every iteration...
			enddo
			! allocate space for rates
			do ih=1,26 ! testing... alloc 26 all 
				allocate(rate(ih)%rcs(4,1))! alloc/deallocate at every iteration...
			enddo
		endif
		! All available hops: 
		!	transition amplitudes and amp^2 for degenerate sectors
		! and allocate space for rates???
		call AllHops()
		write(*,*) "main:   AllHops done... "	
		
		! rates from am2 and energetic penalties
		call CalRates()
		write(*,*) "main:   CalRates done... "	

		! select a hop based on rates
		ih = ihSelect() 
		write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		write(*,*) "main:   icsSelect: ic,is =  ",ic,is

		! perform the transition... update XXXXXXX

		! update na, nx
		itype = itypes(ih,ic);
		na = na + dna(itype);
		nx = nx + dnx(itype);

		if(ic==4) zt = zt+1


		! update occupations, sys%occ, Asites etc
		call UpdateOcc(ih,is)

		call UpdateWays()
		write(*,*) "main:   UpdateWays done... "	

		if (fixmap) then
			call DONTUSEmaps() ! dont reuse basis/hg
		else
		
			call UpdateReqType ! ReqType
			call UpdateMapT ! mapt%map, mapt%cal
			call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
			call UpdateMapB
			write(*,*) "----- updated mapb,mapt--------"
		endif

	enddo ! iter

	write(*,*)"transoc: niter hops done.... " 

	!	do postprocessing....
	write(*,*)"zt = ",zt
	

	! completion message...
	write(*,*)"transoc: everything done.... " 

	stop
	end program
