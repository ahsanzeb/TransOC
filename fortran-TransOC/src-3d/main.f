cc
! NOtes: Assume g, w0 same for all sites
! if not, DPhicreat/annihil, D/Phi hops, all would need seperate sets of eigenstates/values to compute the amplitudes making the comput cot very high.
! might do this in future, but at the moment, PermSym=F can only occur for reasons that do not require diff eig.



	program TransOC
	use modmain
	use init, only: init0, init1, init2, initialise
	use basisstates, only: mkbasis
	use hamiltonian, only: mkHamilt
	use Hoppings, only: AllHops0, AllHops1 ! 
	use rates, only: CalRates
	use selection, only: ihSelect,icsSelect,getpsi2
	use modways, only: UpdateWays,UpdateOcc,UpdateDEQs
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	use modq, only: SetEcoul, qtime
	use mpi
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt,icl
	integer:: nev,ncv,it,ier
	integer :: ntrap,stathc(42,4), stathc0(42,4)
	integer :: totcharge,charge
	double precision:: tottime, dtau
	integer:: trapiter,s,iss,nc,ns,itraj,nex,nelec,idw,ielec
	logical :: found,alloc
	integer :: x(3)
	integer :: nmways(44),nmways0(44) ! N,m,ways(1:42)
	double precision:: rout(10)
	! mpi
	integer:: ierr, ntp, num_procs, node, maxtraj
	double precision, allocatable, dimension (:):: Iav, Iall
	character*50 :: frmt
	integer :: itemp, ig, iEbl, iEbr, nphoton, idv

!======================================================================
	call MPI_INIT(ierr)
	!find out MY process ID, and how many processes were started.
	call MPI_COMM_RANK (MPI_COMM_WORLD, node, ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

	!write(*,*) "Node = ",node," num_procs, ntp = ",num_procs, ntp 
	! readinput file
	call input(node)
	alloc = .false.;

	! find ntp, number of traj per processor/node
	call setntp(ntraj,num_procs,ntp)
	maxtraj = ntp * num_procs; ! could be greater than ntraj
	allocate(Iav(ntp))

	call timestamp(node)
	if(node==0) then
		allocate(Iall(maxtraj))
		Iall = 0.0d0;
		write(frmt,'("(",I6.6,"G18.10)")') maxtraj+3 ! do we really need full out?
	endif

	! allocate space for basis, and calc maps, etc.
	call init0()  ! things that can be used during whole calculation
	!call printnode(node) ! XXX
	do ig=1,ng
		g = gmin + (ig-1)*dg;
		do idw=1,ndw ! detuning
			dw = dwmin + (idw-1)*ddw;
			call setdetuning(dw,detuning,node)
			! allocate space for Hamiltonian and eigstates, etc.
			call init1() ! because changed detuning ===> new Hg, eig.
			do nex=mexmin,mexmax,dmex ! excitations
			ielec=0;
			do nelec = nelmin,nelmax,dnelec
				if (nelec == 0 .and. onlydoped) cycle
				ielec=	ielec+1;
				do idv=1,dev%ntot ! devices: sets of barriers values to make global variables Ebl, Ebr 
				call init2(idv) ! sets Ebr, Ebl etc... 
				do itemp=1,ntemp
				beta = 1/(kb*(tmin + (itemp-1)* dtemp));

				do ier=1,ner ! Er: electric field * r_nns
					Er = Ermin + (ier-1)*dEr;
					if(node==0 .and. mod(ier,1)==0)then
						write(*,'("main: ier = ",i4," of ",i4)')ier,ner
						!write(*,*)"Er, beta, dw, g = ",Er, beta, dw, g
					endif
					!if (VRH) Efieldh = signEr*Er ! VRH= variable range hopping
					! runs ntp trajectories, output => current: Iav, h/c counts: stathc
					call trajectory(ntp,node,Iav,stathc,nmways)
					!write(*,*)"main: after traj"
					! gather all data
					call mpi_gather(Iav, ntp, mpi_double_precision, Iall, ntp, 
     .                mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
					call CheckError(ierr, node, 'gather')
					! sum up all stathc
					call mpi_reduce(stathc,stathc0,42*4, mpi_integer,mpi_sum,
     .                                      0, MPI_COMM_WORLD, ierr)
					call mpi_reduce(nmways,nmways0,44, mpi_integer,mpi_sum,
     .                                      0, MPI_COMM_WORLD, ierr)
					call CheckError(ierr, node, 'reduce')		

					if(node==0) then ! write output 
						call WriteIvsEr(maxtraj,Iall,stathc0)
						call WriteStat(stathc0)
						call writeways(nmways0)
					endif
				enddo ! Er
				enddo ! temp
				enddo ! idv
			enddo ! nelectron
		enddo ! nx
	enddo ! idw
	enddo ! g
	if(node==0) then
		write(*,*)"TransOC: everything done.... " 
		call timestamp(node)
	endif
	call MPI_FINALIZE(ierr)
!======================================================================
	contains
!======================================================================
	subroutine printnode(node)
	implicit none
	integer, intent(in) :: node
	if(node==0) then
		write(*,'("main: this is the master node ",I5)')node
	else
		write(*,'("main: this is the slave node ",I5)') node
	endif
	return
	end 	subroutine printnode
!=============================================	
	subroutine CheckError(ierr, node, oper)
	implicit none
	integer, intent(in):: ierr, node
	character(len=*):: oper
	if(ierr .ne. 0) then
		write(*,'("main: mpi operation ",a," at node ",I3," ier = ",I3)')
     .          oper, node, ierr
		stop
	endif
	
	! comment
	!write(*,'("main: mpi operation ",a," at node ",I3," ier = ",I3)')
  !   .          oper, node, ierr
	return
	end 	subroutine CheckError
!=============================================
	subroutine setntp(ntraj,num_procs,ntp)
	implicit none
	integer, intent(in) :: ntraj,num_procs
	integer, intent(out) :: ntp

	if( mod(ntraj,num_procs)==0 ) then
		ntp = ntraj/num_procs; ! integer div
		if(ntp==0) ntp = 1;
	else
		ntp = ntraj/num_procs + 1; ! integer div
		write(*,*) 'main: ntraj increased to ', ntp * num_procs
	endif
	write(*,*)"main: setntp: ntraj,num_procs,ntp=",ntraj,num_procs,ntp
	

	return
	end 	subroutine setntp
!=============================================
	subroutine setdetuning(dw,detuning,node)
	implicit none
	double precision, intent(inout) :: dw
	logical, intent(inout) :: detuning
	integer, intent(in) :: node
	if(abs(dw)<1.0d-6) then
		detuning = .false.;
		if(node==0) write(*,*) "# main(Warning): dw given, dw used = ",dw, 0.0
		dw = 0.0;
	else
		detuning = .true.;
	endif
	if(node==0) then
		write(*,*) "# dw = ",dw
		write(*,*)
		write(*,*)
	endif 
	return
	end subroutine setdetuning
!=============================================
	subroutine WriteIvsEr(maxtraj,Iall,stathc)
	implicit none
	integer, intent(in):: maxtraj
	double precision, dimension(maxtraj), intent(in) :: Iall
	integer, dimension(42,4), intent(in):: stathc
	open(300,file="I-vs-Er.out",position='append')
	if(ier==1) then
		write(300,'(a,5x,4I10)') "# ier, nx, idw, iel =  ", ier, 
     .         nex, idw, ielec
		write(300,'("     ")')
		write(300,'("     ")')
	endif
	! statistics in modmain: statistics = (/ average, variance^2 /)
	!write(300,frmt) Er, statistics(Iall,ntraj), Iall 
	write(300,'(3G18.10,3x,2I10)') Er, statistics(Iall,maxtraj),
     .                  stathc(25:26,1)
	close(300)
	return
	end subroutine WriteIvsEr
!===============================================
	subroutine WriteStat(stathc)
	implicit none
	integer, dimension(42,4), intent(in):: stathc
	! hops statistics
	open(10,file='stat.out',action='write',position='append')
	write(10,*) '# ier, nx, idw, iel = ', ier, nex, idw, ielec
	write(10,*) 
	if(leads) then
		do i=1,34
			write(10,*)i, stathc(i,:)
		enddo
	else
		do i=1,8
			write(10,*)i, stathc(i,:)
		enddo
		do i=27,34
			write(10,*)i, stathc(i,:)
		enddo
		write(10,*)25, stathc(25,:)
		write(10,*)26, stathc(26,:)
	endif

	if(impurity) then
		do i=35,42
			write(10,*)i, stathc(i,:)
		enddo
	endif

	close(10)
	return
	end subroutine WriteStat

!===============================================
	subroutine writeR()
	implicit none
	
	open(200,file='R.out',action='write',position='append')
	if(itraj == 1) then
		write(200,*) "# idw, ielec",idw,ielec
		write(200,*) 
		write(200,*) 
	endif				
	write(200,*) sum(rout)/(niter*1.d0), rout/(niter*1.d0)
	close(200)
	return
	end 	subroutine writeR
!===============================================
	subroutine writeways(nmways)
	implicit none
	integer,dimension(44), intent(in):: nmways
	double precision:: tot
	tot = niter*ntraj*1.d0;
	open(200,file='ways.out',action='write',position='append')

	if(ier==1) then
		write(200,'(a,5x,4I10)') "# ier, nx, idw, iel =  ", ier, 
     .         nex, idw, ielec
		write(200,'("     ")')
		write(200,'("     ")')
	endif

	if (leads) then
		if(impurity) then
			write(200,'(45G18.10)') Er, nmways/tot
		else
			write(200,'(37G18.10)') Er, nmways(1:36)/tot
		endif
	else ! periodic
		if(impurity)then
			write(200,'(29G18.10)') Er, nmways(1:10)/tot,
     .                        nmways(27:44)/tot
		else
			write(200,'(21G18.10)') Er, nmways(1:10)/tot,
     .                        nmways(27:36)/tot
		endif
	endif
	close(200)
	return
	end subroutine writeways
!======================================================================
! tranejctry runs ntp traj, output: Iav ==> their average currents
!======================================================================
	subroutine trajectory(ntp,node,Iav,stathc,nmways)
	implicit none
	integer, intent(in) :: ntp,node
	double precision, dimension(ntp), intent(out) :: Iav
	integer, dimension(42,4), intent(out) :: stathc
	integer, dimension(44), intent(out):: nmways
	stathc = 0;

	!write(*,'(42G10.5)')"main: Exb = ",Exb

	do itraj=1, ntp
		zt = 0; ntrap = 0;
		totcharge=0; tottime=0.0d0
		trapiter=0

		nmways=0;
		!rout = 0.0d0
		
	!	initialise maps etc
	! set no of active sites and excitations
	nx = nex;
	call initialise(nelec)

	!write(*,*) "main: BondLengths = ",BondLengths

	! put this in init?
	!if(nog ) then
	!	ipsi = 1;
	!	Einit = nx*dw;
	!endif

	!write(*,*)"coulomb: starting traj"
	!call qtime()

	x = 0;
	!========== iterations over number of hops
	! main loop over number of hops asked
	do iter=1,niter

		if(1==0 .and. mod(iter,1)==0 .and. node==0) then
			!write(*,'(a)') ". . . . . . . . . . . "
			write(*,'(a,i10,a,i10,a,i10,a,2i10)')"Node = ",node,
     . " itraj = ",itraj," iter = ",iter,"  N, m = ",na,nx
			!write(*,*)'Ebl, Ebr = ',Ebl, Ebr     
			write(*,'(a,3i3)')'N_d, N_phi, N_a = ', sys%n2,sys%n0,sys%n1
		endif


		!write(*,*) iter, na, nx
		!write(*,*) "in: Asites = ",Asites
		!write(*,'(3i1,x,3i1,x,3i1)')sys%occ

		nmways(1:2) = nmways(1:2) + (/na,nx/);
		nmways(3:44) = nmways(3:44) + ways(:)%ns;
		
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
			!write(*,*)'start: Nphoton = ', Nphoton(ipsi,1)

		endif

		if (.not. nog) then ! if same quantum state as last iter, why not avoid this alloc/copy etc.??
			!write(*,*)"main: not nog.... "	
			if(allocated(psi)) deallocate(psi)
			it = mapt%map(1);
			allocate(psi(1,eig(it)%n1))
			psi(1,:) = eig(it)%evec(:,1)
			Einit = eig(it)%eval(1)
			!write(*,*)"main: Einit block ends...."	
		endif

		if(debug)write(*,*) "main:   psi,Einit update done... "
		!write(*,'(a,i5,x,i10,x,f10.5)')
    ! .      "main: it,ntot, Ei = ",it,eig(it)%n1,Einit

		!write(*,*)"main: psi="	,psi
		!write(*,*)"main: eval 1="	,eig(it)%eval
		if(debug)then
			write(*,*)'--------------------------------------'
			i=1
			write(*,'(3i2)')sys%occ(i),sys%occ(i+3),sys%occ(i+6)
			i=2
			write(*,'(3i2)')sys%occ(i),sys%occ(i+3),sys%occ(i+6)
			i=3
			write(*,'(3i2)')sys%occ(i),sys%occ(i+3),sys%occ(i+6)
		endif
	
		if(.not. alloc) then
			! at most nsites ways for any hop???
				do ia=1,21
					allocate(qt(ia)%cs(4,sys%nsites))
					if (PermSym .and. ia > 9 ) exit
				enddo
				! allocate space for rates
				do ih=1,nproc ! testing... alloc 26 all 
					if(vrh .or. coulomb .or. (.not. PermSym))then
						allocate(rate(ih)%rcs(4,sys%nsites))
					else !if(PermSym) then 
						allocate(rate(ih)%rcs(4,1))
					endif
				enddo
				alloc = .true.
		endif

		do ih=1,nproc
			rate(ih)%rcs(:,:) = 0.0d0
			rate(ih)%r = 0.0d0
		enddo

	!if(nog)write(*,*) "main: nog=T; ipsi =  ",ipsi

		! Before calculating the rates:
		! Net charge on the system
		! nelec =  doping given in input; +ve for holes, -ve for electrons
		Qnet = (nsites-nelec)-sum(sys%occ);
		! charging energies for contact hops
		call UpdateDEQs(Qnet)	
		if(debug)write(*,*) "main:   UpdateDEQs done... "

		!write(*,*)"sys%occ = ", sys%occ
		!write(*,*)"main: ways = ",ways(:)%ns
		



	if (leads)
     . x = x + (/nelec,sum(sys%occ), Qnet /);


		! All available hops: 
		!	transition amplitudes and amp^2 for degenerate sectors
		! and allocate space for rates???
		s = 1;
		if(PermSym) then
			call AllHops0()
		else
			call AllHops1()
		endif
		if(debug)write(*,*) "main:   AllHops done... "

		!write(*,*)"main: not nog.... "	


		! Coulomb's interaction
		! sets Ecoul(ih=1:34)%dEq(is=1:ns(ih)) and Etotq for use in rates
		call SetEcoul(node,iter)
		if(debug)write(*,*) "main:   SetEcoul done... "	

		
		! rates from am2 and energetic penalties
		call CalRates()
		if(debug)write(*,*) "main:   CalRates done... "	

		!write(*,*) "main: ===1==> ",rate(:)%r
		!if(na==nsites) write(*,*)"main: N,m, R = ", na,nx,sum(rate(:)%r)

		!write(*,*)"main: ER, Rtot = ",Er, sum(rate(:)%r)

		! select a hop based on rates
		ih = ihSelect();

		if(ih == -1) then
			write(*,*)"main: R too small, so setting I=0 "
			Iav(itraj) = 0.0d0;
			exit ! exit iter loop; start next trajectory
		endif

		if(1==0 .and. na==nsites-1) then
			write(*,*)"main: Rtot = ",sum(rate(:)%r)
			write(*,*)"main: R_L/R = ",rate(1:8)%r
			write(*,*)"main: R_U/D = ",rate(27:34)%r
			write(*,*)"main: Rkappa, Rgamma= ",rate(25:26)%r
		endif
		!write(*,*) "main:   ih = ",ih
		!write(*,*) "main: ===2==> ",rate(:)%r
		if(debug)write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		if(debug)write(*,*) "main:   icsSelect: ic,is =  ",ic,is
		!write(*,*) "main:   ih,ic,is = ",ih,is

		!rout(1:8) = rout(1:8) + rate(1:8)%r 
		!rout(9:10) = rout(9:10) + rate(25:26)%r 
		
		! if nog, the final state index from selected ih,ic
		if (nog) then
			iss=is; ! forgot why iss was seperately defined months ago
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

			if(1==0)then
			write(*,*)'main: Nog: ih,ic =',ih,ic
			!write(*,*)'main: getNphoton arg ipsi,itype=',ipsi,itype
			call getNphoton(ipsi,itype,Nphoton)
			write(*,*)'main: Nphoton=',Nphoton
			endif
		endif

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
		dtau=1.0d0/sum(rate(1:nproc)%r); ! average used. 
		!can be picked from adistribution at random using log(eta)/R
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

	Iav(itraj) = totcharge/tottime

	!write(*,*)"main: dQ, dt, Iav = ", totcharge, tottime,Iav(itraj)
	enddo ! itraj



	!write(*,*)"main: end traj"	
	!call qtime()
	

	return
	end subroutine trajectory
!======================================================================

	end program
