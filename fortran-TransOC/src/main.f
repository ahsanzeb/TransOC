
	program TransOC
	use modmain
	use init, only: initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: mkHamilt
	use Hoppings, only: AllHops ! 
	use rates, only: CalRates
	use selection, only: ihSelect, icsSelect
	use modways, only: UpdateWays, UpdateOcc
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt=0
	integer:: nev,ncv,it
	integer :: stath(26),statc(4)


	stath = 0; statc = 0; zt = 0;

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

		write(*,'(a)') ". . . . . . . . . . . . "
		write(*,'(a,i10,a,2i10)') " iter = ",iter,"  N, m = ",na,nx
		write(*,*)          
		
		! make basis states ksub
		call mkbasis(na,nx)
		!write(*,*) " basis done....  "

		! make hamiltonian
		call mkHamilt()
		!write(*,*) "main: Hamiltonians done....  "

		call diagonalise()
		!write(*,*) "main: diagonalisation done.... "

		!stop
		
		if(allocated(psi)) deallocate(psi)
		it = mapt%map(1);
				
		allocate(psi(1,eig(it)%n1))
		psi(1,:) = eig(it)%evec(:,1)
		Einit = eig(it)%eval(1)

		write(*,'(a,i5,x,i10,x,f10.5)')
     .      "main: it,ntot, Ei = ",it,eig(it)%n1,Einit

		!write(*,*)"main: psi="	,psi
		!write(*,*)"main: eval 1="	,eig(it)%eval
		
		write(*,*)"sys%occ = ",sys%occ


		!====== if PermSym then only a single site cases ====== 
		!	check if anything to do about it related to qt/rate allocation etc
		
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
		!write(*,*) "main:   AllHops done... "	
		
		! rates from am2 and energetic penalties
		call CalRates()
		!write(*,*) "main:   CalRates done... "	

		! select a hop based on rates
		ih = ihSelect() 
		!write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		!write(*,*) "main:   icsSelect: ic,is =  ",ic,is

		! perform the transition... update XXXXXXX

		! update na, nx
		itype = itypes(ih,ic);
		na = na + dna(itype);
		nx = nx + dnx(itype);

		if(ic==4) zt = zt+1


		stath(ih) = stath(ih) + 1;
		statc(ic) = statc(ic) + 1;

		write(*,'(a,10i3)')"stath = ",stath(1:8),stath(25:26)
		write(*,'(a,4i5)')"statc = ",statc
		!write(*,'(a,i5)')"zt = ",zt

		! update occupations, sys%occ, Asites etc
		call UpdateOcc(ih,is)

		call UpdateWays()
		!write(*,*) "main:   UpdateWays done... "	

		if (fixmap) then
			call DONTUSEmaps() ! dont reuse basis/hg
		else
		
			call UpdateReqType ! ReqType
			call UpdateMapT ! mapt%map, mapt%cal
			call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
			call UpdateMapB
			!write(*,*) "----- updated mapb,mapt--------"
		endif

	enddo ! iter

	write(*,*)"transoc: niter hops done.... " 

	!	do postprocessing....
	write(*,*)"zt = ",zt
	

	! completion message...
	write(*,*)"transoc: everything done.... " 

	end program
