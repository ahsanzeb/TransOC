
	program TransOC
	use modmain
	use init, only: initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: mkHamilt
	use Hoppings, only: AllHops ! 
	use rates, only: CalRates, setrates
	use selection, only: ihSelect,icsSelect,getpsi2
	use modways, only: UpdateWays,UpdateOcc,UpdateDEQs
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	use lindblad, only: initme, mesolve
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt,icl
	integer:: nev,ncv,it
	integer :: stath(26),statc(4),ntrap,stathc(26,4)
	integer :: totcharge,charge
	double precision:: tottime, dtau
	integer:: trapiter,s,iss,nc,ns,itraj,nex,nelec,idw,ielec
	logical :: found,alloc
	integer :: x(3)
	integer :: nms(2), wayss(26)
	double precision:: rout(10)
	
	! time stamp & random seed 
	!call timestamp

	!stath = 0; statc = 0; 
	alloc = .false.
	master = .false.
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

	write(*,*) "# dw = ",dw
	write(*,*)
	write(*,*)

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

		nms = 0; wayss = 0;
		rout = 0.0d0
		
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

		if(1==0) then
			write(*,'(a)') ". . . . . . . . . . . . "
			write(*,'(a,i10,a,i10,a,2i10)') " itraj = ",itraj,
     .               " iter = ",iter,"  N, m = ",na,nx
			write(*,*)          
		endif


		!write(*,*) iter, na, nx

		
		nms = nms + (/na,nx/);
		wayss = ways(:)%ns;
		
		!if(mod(iter-1,20)==0) write(*,*)"wayss=", wayss
		
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
					allocate(qt(ia)%cs(4,sys%nsites))
				enddo
				! allocate space for rates
				do ih=1,26 ! testing... alloc 26 all 
					if(PermSym)then
						allocate(rate(ih)%rcs(4,1))
						allocate(mrate(ih)%rcst(4,1,ntcoarse))
					else
						allocate(rate(ih)%rcs(4,sys%nsites))
						allocate(mrate(ih)%rcst(4,sys%nsites,ntcoarse))
					endif
					allocate(mrate(ih)%r(ntcoarse))
				enddo
				alloc = .true.
		endif

		do ih=1,26
			rate(ih)%rcs(:,:) = 0.0d0
			rate(ih)%r = 0.0d0
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
		call AllHops()
		if(debug)write(*,*) "main:   AllHops done... "	
		
		! rates from am2 and energetic penalties
		call CalRates()
		if(debug)write(*,*) "main:   CalRates done... "	

		if ( sum(rate(:)%r) < 1.0d-25 ) then
		!if ( mod(iter,3)==0) then
			!write(*,*)" ===== * ===== * ===== * ===== * ===== *"

			write(*,*) "main: trap encountered! , iter = ",iter
			
			if (nog) then
				goto 99
			endif

			write(*,*)"main: sum(rate(:)%r) = ",sum(rate(:)%r)
			write(*,*)"main: r = ",rate(1:8)%r

			ntrap = 	ntrap +1;
			!  Rt , delta ?!
			! set to 0.0d0
			do ih=1,26
				mrate(ih)%rcst = 0.0d0
				mrate(ih)%r = 0.0d0		
			enddo

			! what if evolution 1/R > tmax that rho is evolved up to?
			!	if lowest eigenstate is a trap, then higher states
			! will still be decaying to it, so the dtau=1/R needs
			! to consider this..... 
			write(*,*) "mesolve called... "
			call mesolve()

			write(*,*) "mesolve done... "

			! then go to AllHops() with master=.true.
			! and after than call findtau() instead of CalRates()
			! then things are usual.
			master = .true.
			call AllHops()
			write(*,*) "Allhops with master=T done... "

			master = .false.
			call setrates();
			write(*,*) "setrates done... "
			write(*,*)"main: mesolve: sum(rate(:)%r) = ",sum(rate(:)%r)

		endif


		!write(*,*) "main: ===1==> ",rate(:)%r
	
		! select a hop based on rates
		ih = ihSelect() 

		!if(na==10) write(*,*)"R = ",sum(rate(:)%r),rate(1:8)%r


		!write(*,*) "main: ===2==> ",rate(:)%r
		if(debug)write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		if(debug)write(*,*) "main:   icsSelect: ic,is =  ",ic,is


		rout(1:8) = rout(1:8) + rate(1:8)%r 
		rout(9:10) = rout(9:10) + rate(25:26)%r 
		
		! saves amp etc for later use in case it's a trap  
		call initme(ih,ic,is) 
		
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
		dtau=1.0d0/sum(rate(:)%r)
		totcharge = totcharge + charge
		tottime = tottime + dtau
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
	!------------------------------------------------------------ 
	open(200,file='R.out',action='write',position='append')
	if(itraj == 1) then
		write(200,*) "# idw, ielec",idw,ielec
		write(200,*) 
		write(200,*) 
	endif				
	write(200,*) sum(rout)/(niter*1.d0), rout/(niter*1.d0)
	close(200)

	
	open(200,file='ways.out',action='write',position='append')
	if(itraj == 1) then
		write(200,*) "# idw, ielec",idw,ielec
		write(200,*) 
		write(200,*) 
	endif				
	if (leads) then
		write(200,*) nms/(niter*1.d0), wayss/(niter*1.d0)
	else
		write(200,*) nms/(niter*1.d0), wayss(1:8)/(niter*1.d0)
	endif
	close(200)
	!------------------------------------------------------------
	open(200,file='current.out',action='write',position='append')
	if(itraj == 1) then
		write(200,*) "# nx, idw, ielec ",nex, idw, ielec
		write(200,*) 
		write(200,*) 
	endif				
	!write(200,*) !'(i5,5x,2f15.10)')
  !   . totcharge, tottime, net charge on the system
	if (leads) then
  		! Net charge on the system if contact are present
		! intrinsic #electrons - present #electrons
		write(200,*) totcharge, tottime, x/(niter*1.d0)
	else
		write(200,*) totcharge, tottime
	endif
	close(200)
	!else ! will see it later
	!------------------------------------------------------------
	enddo ! itraj
	!********************** trajectories *************

	! hops statistics
	open(10,file='stat.out',action='write',position='append')
	write(10,*) '# nx, idw, iel = ', nex, idw, ielec
	write(10,*) 
	if(leads) then
		do i=1,24
			write(10,*)stathc(i,:)
		enddo
	else
		do i=1,8
			write(10,*)stathc(i,:)
		enddo
	endif
	write(10,*)stathc(25,:)
	write(10,*)stathc(26,:)
	close(10)

	
	enddo ! ielectron
	enddo ! idw
	enddo ! nex
	!	do postprocessing....
	!write(*,*)"zt = ",zt


	! completion message...
	write(*,*)"transoc: everything done.... " 

	end program
