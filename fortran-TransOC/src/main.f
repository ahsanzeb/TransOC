
	program TransOC
	use modmain
	use init, only: initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: mkHamilt
	use Hoppings, only: AllHops ! 
	use rates, only: CalRates
	use selection, only: ihSelect,icsSelect,getpsi2
	use modways, only: UpdateWays,UpdateOcc,UpdateDEQs
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt,icl
	integer:: nev,ncv,it
	integer :: stath(26),statc(4),ntrap,stathc(26,4)
	integer :: totcharge,charge
	double precision:: tottime, dt
	integer:: trapiter,s,iss,nc,ns,itraj,nex,nelec,idw,ielec
	logical :: found,alloc
	integer :: x(3)
	! time stamp & random seed 
	call timestamp

	!stath = 0; statc = 0; 
	alloc = .false.
	
	! readinput file
	call input()

	if(nog) write(*,*)" **** No Coupling! **** "
	!------------------------------------------------------

	do nex=mexmin,mexmax,dmex
	!nex = nx; ! copy for each trajectory

	! detuning
	do idw=1,ndw
	dw = dwmin + (idw-1)*ddw;
	ielec=0
	do nelec = nelmin,nelmax,dnelec
		if (nelec == 0 .and. onlydoped) cycle
		ielec=	ielec+1;

		!write(*,*) 'dw, nelec = ', dw, nelec
		stathc = 0;
		do itraj=1,ntraj
	!********************** trajectories *************
		zt = 0; ntrap = 0;
		totcharge=0; tottime=0.0d0
		trapiter=0

	!	initialise maps etc
	! set no of active sites and excitations
	nx = nex;
	call initialise(nelec)

	! put this in init?
	!if(nog ) then
	!	ipsi = 1;
	!	Einit = nx*dw;
	!endif

	x = 0;
	!========== iterations over number of hops
	! main loop over number of hops asked
	do iter=1,niter

		!write(*,'(a)') ". . . . . . . . . . . . "
		!write(*,'(a,i10,a,i10,a,2i10)') " itraj = ",itraj,
    ! .               " iter = ",iter,"  N, m = ",na,nx
		!write(*,*)          

		! make basis states ksub
		call mkbasis(na,nx)
		if(debug)write(*,*) " basis done....  "
		
		! make hamiltonian
		call mkHamilt()
		if(debug)write(*,*) "main: Hamiltonians done....  "

		call diagonalise()
		if(debug)write(*,*) "main: diagonalisation done.... "
		
		!uncoupled case: start with excited sites
		if(nog .and. iter==1) then
			call chooseipsi(ipsi,Einit);
			!ipsi = eig(mapt%map(1))%ntot
			!Einit = eig(mapt%map(1))%eval(ipsi)
		endif

		if (.not. nog) then
			if(allocated(psi)) deallocate(psi)
			it = mapt%map(1);
			allocate(psi(1,eig(it)%n1))
			psi(1,:) = eig(it)%evec(:,1)
			Einit = eig(it)%eval(1)
		endif

		!write(*,'(a,i5,x,i10,x,f10.5)')
    ! .      "main: it,ntot, Ei = ",it,eig(it)%n1,Einit

		!write(*,*)"main: psi="	,psi
		!write(*,*)"main: eval 1="	,eig(it)%eval
		
		!write(*,*)"sys%occ = ",sys%occ

		if(.not. alloc) then
			! at most nsites ways for any hop???
				do ia=1,14
					allocate(qt(ia)%cs(4,sys%nsites)) ! alloc/deallocate at every iteration...
				enddo
				! allocate space for rates
				do ih=1,26 ! testing... alloc 26 all 
					allocate(rate(ih)%rcs(4,sys%nsites))! alloc/deallocate at every iteration...
				enddo
				alloc = .true.
		endif

		do ih=1,26
			rate(ih)%rcs(:,:) = 0.0d0
		enddo

	!if(nog)write(*,*) "main: nog=T; ipsi =  ",ipsi


		! Before calculating the rates:
		! Net charge on the system
		Qnet = (nsites-nelec)-sum(sys%occ);
		! charging energies for contact hops
		call UpdateDEQs(Qnet, nelec)	

	if (leads)
     . x = x + (/nelec,sum(sys%occ),(nsites-nelec)-sum(sys%occ)/);


		! All available hops: 
		!	transition amplitudes and amp^2 for degenerate sectors
		! and allocate space for rates???
		s = 1;
10		call AllHops()
		if(debug)write(*,*) "main:   AllHops done... "	
		
		! rates from am2 and energetic penalties
		call CalRates()
		if(debug)write(*,*) "main:   CalRates done... "	

		if ( sum(rate(:)%r) < 1.0d-14 ) then
			write(*,*)"main: sum(rate(:)%r) = ",sum(rate(:)%r)
			
			! to avoid mem error if trap on first iteration
			!if(iter==1 .and. (.not. nog)) call getpsi2(1,1,1); 
			
			if (trapiter == iter) then
				s = s + 1;
				write(*,*) "main: trap for excited sector",s-1," too"
				write(*,*) "main: s, eig(it)%nsec=",s, eig(it)%nsec
				if (s > eig(it)%nsec) then
					write(*,*)"main: no more sectors! Aborting!"
					goto 99
				endif
			else
				write(*,*) "main: trap encountered! , iter = ",iter
				s = 2; ! second sec =  first excited sec
				ntrap = 	ntrap +1;
				trapiter = iter;
			endif		
			if (nog) then
				write(*,*) "main: trap encountered! , iter = ",iter
				goto 99
			endif
			! psi2, Einit2, 
			!	in case the lowest state psi sees a trap,
			! we will use this state to get us out
			if(s <= eig(it)%nsec) then
				if(iter==1) then
					Einit = eig(mapt%map(1))%eval(s);
					psi(1,:) = 0.0d0;
					!random superposition 
					do i=eig(it)%ind(s),eig(it)%ind(s+1)-1 ! second degenerate sector
						psi(1,:) = psi(1,:) + eig(it)%evec(:,i) * rand()
					enddo
					! normalise ?
					psi(1,:) = psi(1,:)/sum(psi(1,:)**2);			
				else
					psi(1,:) = psi2(s,:);
					Einit = Einit2(s);
				endif
			else
				write(*,*)"main: s <= eig(it)%nsec ! Aborting..."
				stop
			endif
			goto 10 ! calculate amplitudes and rates again
							!	with this excited state
			! if needed, can we keep repeating this
			! for higher and higher degen sectors
			! until we are out of trap?
		elseif(s>1 .and. (.not. nog))then
				write(*,*) "main: got out of trap, s=",s	
				s = 0;	
		endif


		!write(*,*) "main: ===1==> ",rate(:)%r
	
		! select a hop based on rates
		ih = ihSelect() 

		!write(*,*) "main: ===2==> ",rate(:)%r
		if(debug)write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		if(debug)write(*,*) "main:   icsSelect: ic,is =  ",ic,is

		! if nog, the final state index from selected ih,ic
		if (nog) then
			if(PermSym) then
				iss = 1;
			else
				iss = is
			endif
			
			itype = itypes(ih,ic); ! type for ih,ic selected
			it =mapt%map(itype); ! location of types
			ia = maph(ih,ic); ! location of amplitudes
			icl=mapc(ih,ic);
			!write(*,*)"main: amp",qt(ia)%cs(ic,iss)%amp2
			found=.false.
			do i=1,Hg(it)%ntot
				!if(allocated(qt(ia)%cs(ic,iss)%amp))then
				!else
				!	write(*,*)"main: amp not allc: ",itype,it,ia
				!endif
				if(qt(ia)%cs(icl,iss)%amp2(i) > 0.5d0) then ! amp = 0.0 or 1.0
					ipsi = i;
					Einit = eig(it)%eval(i)
					found=.true.
					exit
				endif
			enddo

			if(.not. found) then
				write(*,*)"main: ipsi not found!!!"
				write(*,*)	"main: sum(amp2)=", sum(qt(ia)%cs(icl,iss)%amp2)
				stop
			endif
			
		endif





		!write(*,*) "main: ===3==> "!,rate(:)%r

		! to get out of traps, keep excited states.
		! psi2, Einit2, 
		!	in case the lowest state psi sees a trap,
		! we will use this state to get us out
		if (.not. nog) call getpsi2(ih,ic,is);
		if(debug)write(*,*) "main: getpsi2 done...."
		!write(*,*) "main: ---3--"
		!write(*,*) "main: =====> ",rate(:)%r
		!write(*,*) "main: =====> "
		if(debug)write(*,*) "main: writeout done...."
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
		elseif(leads) then ! just monitor the left contact
			charge = 0
			select case(ih)
				case(15,16,21,23)
					charge = +1
				case(11,12,22,24)
					charge = -1
			end select	
		endif
		dt=1.0d0/sum(rate(:)%r)
		totcharge = totcharge + charge
		tottime = tottime + dt
		!-------------------------------

		! update na, nx
		itype = itypes(ih,ic);
		na = na + dna(itype);
		nx = nx + dnx(itype);

		if(ic==4) zt = zt+1


		!stath(ih) = stath(ih) + 1;
		!statc(ic) = statc(ic) + 1;
		stathc(ih,ic) = stathc(ih,ic) + 1

		!write(*,'(a,10i10)')"stath = ",stath(1:8),stath(25:26)
		!write(*,'(a,4i10)')"statc = ",statc
		!write(*,'(a,i10)')"traps = ",ntrap
		!write(*,'(a,i5)')"zt = ",zt


		! write grand total and total rates for all hop types
		! and zt, ntraps
		if(iter==niter) then
			call writeout(.true.,zt,ntrap)
		else
			call writeout(.false.,0,0)
		endif



		! update occupations, sys%occ, Asites etc
		call UpdateOcc(ih,is)
		if(debug)write(*,*) "main: UpdateOcc done...."

		call UpdateWays()
		if(debug)write(*,*) "main:   UpdateWays done... "	
		
		if (fixmap) then
			call DONTUSEmaps() ! dont reuse basis/hg
		else
		
			call UpdateReqType ! ReqType
			call UpdateMapT ! mapt%map, mapt%cal
			call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
			call UpdateMapB
			if(debug)write(*,*) "----- updated mapb,mapt--------"
		endif

	enddo ! iter
	!write(*,*)"transoc: niter hops done.... " 
	!===================================================== 

99		continue


	open(200,file='current.out',action='write',position='append')
	if(itraj == 1) then
		write(200,*) "# idw, ielec",idw,ielec
		write(200,*) 
		write(200,*) 
	endif				
	!write(200,*) !'(i5,5x,2f15.10)')
  !   . totcharge, tottime, net charge on the system
	if (leads) then
  		! Net charge on the system if contact are present
		! intrinsic #electrons - present #electrons
		write(200,*) totcharge, tottime, x/niter
	else
		write(200,*) totcharge, tottime
	endif
	close(200)
	!else ! will see it later





	enddo ! itraj
	!********************** trajectories *************

	! hops statistics
	open(10,file='stat.out',action='write',position='append')
	if(leads) then
		do i=1,24
			write(10,*)stathc(i,:)
		enddo
	else
		do i=1,8
			write(10,*)stathc(i,:)
		enddo
		write(10,*)stathc(25,:)
		write(10,*)stathc(26,:)
	endif
	close(10)

	
	enddo ! ielectron
	enddo ! idw
	enddo ! nex
	!	do postprocessing....
	!write(*,*)"zt = ",zt


	! completion message...
	write(*,*)"transoc: everything done.... " 

	end program
