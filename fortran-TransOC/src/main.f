
	program TransOC
	use modmain
	use init, only: initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: mkHamilt
	use Hoppings, only: AllHops ! 
	use rates, only: CalRates
	use selection, only: ihSelect,icsSelect,getpsi2
	use modways, only: UpdateWays, UpdateOcc
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt
	integer:: nev,ncv,it
	integer :: stath(26),statc(4),ntrap,stathc(26,4)
	integer :: totcharge,charge
	double precision:: tottime, dt
	integer:: trapiter,s

	stath = 0; statc = 0; zt = 0; ntrap = 0;
	stathc = 0;
	
	totcharge=0; tottime=0.0d0
	trapiter=0
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
10		call AllHops()
		!write(*,*) "main:   AllHops done... "	
		
		! rates from am2 and energetic penalties
		call CalRates()
		!write(*,*) "main:   CalRates done... "	


		if ( sum(rate(:)%r) < 1.0d-10) then
			!write(*,*)"main:sum(rate(:)%r) < 1.0d-10 "
			if (trapiter == iter) then
				s = s + 1;
				write(*,*) "main: trap for sector excited sector",s-1," too"
				write(*,*) "main: s, eig(it)%nsec=",s, eig(it)%nsec
				if (s > eig(it)%nsec) then
					write(*,*)"main: no more sectors! Aborting!"
					stop
				endif
			else
				write(*,*) "main: trap encountered! , iter = ",iter
				s = 2; ! second sec =  first excited sec
				ntrap = 	ntrap +1;
				trapiter = iter;
			endif		
			! psi2, Einit2, 
			!	in case the lowest state psi sees a trap,
			! we will use this state to get us out
			if(s <= eig(it)%nsec) then
				psi(1,:) = psi2(s,:);
				Einit = Einit2(s);
			else
				write(*,*)"main: s <= eig(it)%nsec ! Aborting..."
				stop
			endif
			goto 10 ! calculate amplitudes and rates again
							!	with this excited state
			! if needed, can we keep repeating this
			! for higher and higher degen sectors
			! until we are out of trap?
		elseif(s>1)then
				write(*,*) "main: got out of trap, s=",s	
				s = 0;	
		endif


		!write(*,*) "main: ===1==> ",rate(:)%r
	
		! select a hop based on rates
		ih = ihSelect() 

		!write(*,*) "main: ===2==> ",rate(:)%r

		!write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		!write(*,*) "main:   icsSelect: ic,is =  ",ic,is

		!write(*,*) "main: ===3==> "!,rate(:)%r

		! to get out of traps, keep excited states.
		! psi2, Einit2, 
		!	in case the lowest state psi sees a trap,
		! we will use this state to get us out
		call getpsi2(ih,ic,is);

		!write(*,*) "main: ---3--"
		!write(*,*) "main: =====> ",rate(:)%r
		!write(*,*) "main: =====> "
		! write grand total and total rates for all hop types
		call writeout()
		! perform the transition... update XXXXXXX
		!write(*,*) "main: ---4--"


		!----------------------------------------------
		!	wirte current file: delta charge, delta t
		!----------------------------------------------
		if(periodic .or. onlybulk) then
			charge = 0
			select case(ih)
				case(1,4,5,7)
					charge = +1
				case(2,3,6,8)
					charge = -1
			end select
			dt=1.0d0/sum(rate(:)%r)
			totcharge = totcharge + charge
			tottime = tottime + dt
			open(200,file='current.out',action='write',position='append')
			write(200,*) !'(i5,5x,2f15.10)')
     .    charge, dt, totcharge*1.0d0/tottime ! so far
			close(200)
			!else ! will see it later
		endif
		!-------------------------------

		! update na, nx
		itype = itypes(ih,ic);
		na = na + dna(itype);
		nx = nx + dnx(itype);

		if(ic==4) zt = zt+1


		stath(ih) = stath(ih) + 1;
		statc(ic) = statc(ic) + 1;
		stathc(ih,ic) = stathc(ih,ic) + 1

		write(*,'(a,10i10)')"stath = ",stath(1:8),stath(25:26)
		write(*,'(a,4i10)')"statc = ",statc
		write(*,'(a,i10)')"traps = ",ntrap
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


	! hops statistics
	do i=1,8
		write(*,'(a,i5,5x,4i10)')'ih = ',i,stathc(i,:)
	enddo



	!	do postprocessing....
	write(*,*)"zt = ",zt
	

	! completion message...
	write(*,*)"transoc: everything done.... " 

	end program
