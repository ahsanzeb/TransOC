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
	use rates, only: CalRates,dhopsrate
	use selection, only: ihSelect,icsSelect,getpsi2
	use modways, only: UpdateWays,UpdateOcc,UpdateDEQs
	use readinput, only: input
	use maps
	use diag, only: diagonalise
	use modq, only: SetEcoul
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

				!write(*,*)'main: dev%ntot=',dev%ntot

				do idv=1,dev%ntot ! devices: sets of barriers values to make global variables Ebl, Ebr 
				call init2(idv) ! sets Ebr, Ebl etc... 
				!write(*,*)'main: item loop  . . . . . . . ,ntemp=',ntemp

				do itemp=1,ntemp
				beta = 1/(kb*(tmin + (itemp-1)* dtemp));


				do ier=1,ner ! Er: electric field * r_nns
				!write(*,*)'main: ier loop  . . . . . . . ,ner=',ner

					Er = Ermin + (ier-1)*dEr;
					if(1==0 .and. node==0 .and. mod(ier,1)==0)then
						write(*,'("main: ier = ",i4," of ",i4)')ier,ner
						!write(*,*)"Er, beta, dw, g = ",Er, beta, dw, g
					endif
					!if (VRH) Efieldh = signEr*Er ! VRH= variable range hopping
					! runs ntp trajectories, output => current: Iav, h/c counts: stathc
					!write(*,*)'main: calling trajectory()  . . . . . . . '
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
		if(1==0 .and. node==0) 
     .  write(*,*) "# main(Warning): dw given, dw used = ",dw, 0.0
		dw = 0.0;
	else
		detuning = .true.;
	endif
	if(1==0 .and. node==0) then
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
	integer :: jsec,itermax

	stathc = 0;

	!write(*,'(42G10.5)')"main: Exb = ",Exb
	!write(*,*)'PermSym = ',PermSym
	
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
	x = 0;
	!========== iterations over number of hops
	! main loop over number of hops asked
	itermax = niter
	do iter=1,niter
	
		if(node==0 .and. na .ge. 10)write(*,'(a,i10)')" sector = ",iter

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


		if (.not. nog) then
			if(allocated(psi)) deallocate(psi)
			it = mapt%map(1);
			! maz: niter as dummy index to get amp FROM various eigenstates
			!write(*,*)'iter,  size(eig(it)%eval) = ',
     !.                          iter,  eig(it)%nsec
			if(iter > eig(it)%n2)then !.or. eig(it)%eval(iter) > 0.d0)then
				itermax = iter-1
				!write(*,*)'1: im,jsec = ',itermax,jsec
			 exit
			endif
			allocate(psi(1,eig(it)%n1))
			psi(1,:) = eig(it)%evec(:,iter)
			Einit = eig(it)%eval(iter);
			!if( Einit > 0.0d0) exit; ! only 0 and below states

			if(.not. allocated(drates)) then
				allocate(drates(eig(it)%nsec,eig(it)%nsec,2))
				drates = 0.0d0
			endif

			jsec = 0;
			do i=1,eig(it)%nsec
						if(iter .ge. eig(it)%ind(i) .and. 
     .         iter .lt. eig(it)%ind(i+1)) then
								jsec = i;
								exit
     				endif
			end do
			!write(*,*)'ind=',eig(it)%ind(1:eig(it)%nsec+1)
			!write(*,*)'iter, jsec = ',iter, jsec
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

		! All available hops: 
		!	transition amplitudes and amp^2 for degenerate sectors
		! and allocate space for rates???

		!write(*,*) "main:   AllHops starts... "
		s = 1;
		if(PermSym) then
			call AllHops0()
		else
			call AllHops1()
		endif

		! rates from am2 and energetic penalties
		!call CalRates()
		call dhopsrate(jsec)
		
	enddo ! iter

	!write(*,*)'2: im,jsec = ',itermax,jsec
	call writedrates(itermax,jsec)

	Iav(itraj) = 0.0d0

	enddo ! itraj



	!write(*,*)"main: end traj"	

	return
	end subroutine trajectory
!======================================================================
	subroutine writedrates(im,jsec)
	implicit none
	integer, intent(in) :: im,jsec
	integer :: it, nstates, jsecs
	logical :: fullout

		fullout = .false.

		! amplitudes^2 for j to j'/=j are the same for H-H and L-L channels
		! so we can sum them up, and write
		! nstates, Es, a2j+, a2j-, a2jj'max+-, a2jj'tot+-
		open(10,file='degen-sec-amp2.dat',
     .               action='write',position='append')
		if(fullout) then
		open(100,file='degen-sec-amp2-1.out',
     .               action='write',position='append')
		open(101,file='degen-sec-amp2-2.out',
     .               action='write',position='append')
		open(102,file='degen-sec-amp2-tot.out',
     .               action='write',position='append')
		endif

	jsecs = 1;
	it = mapt%map(1);
	! normalise with number of states in every dgenerate sector
	! to get average rate per state
	do i=1,jsec; !eig(it)%nsec
		! im: to fix: possibly incorrect average for the highest sector.
		! not all states of the highest sec may have been used in iter loop in trajectory().
		nstates = 0
		if(im .lt. eig(it)%ind(i+1)-1) then
			nstates = im - eig(it)%ind(i) + 1
			!write(*,*)'1 nstates = ',nstates
		else
			nstates = eig(it)%ind(i+1) - eig(it)%ind(i)
			!write(*,*)'2 nstates = ',nstates
		endif
		!drates(i,:,:) = 	drates(i,:,:)/nstates 
		!write(*,*)'im, isec, nstates=',im, i, nstates
		!write(*,*)'i1, i2=',eig(it)%ind(i), eig(it)%ind(i+1)
		if(fullout)then
		write(100,'(1000f15.8)') eig(it)%esec(i), drates(i,:,1)/nstates
		write(101,'(1000f15.8)') eig(it)%esec(i), drates(i,:,2)/nstates
		write(102,'(1000f15.8)') eig(it)%esec(i), 
     .                              sum(drates(i,:,:),dim=2)/nstates
		endif
		if (i==1) then
			write(10,'(i10,3x,1000f15.8)') nstates, eig(it)%esec(i), 
     . drates(i,i,1)/nstates, drates(i,i,2)/nstates, 0.0d0, 
     . sum(drates(i,1:i,1:2))/nstates
		elseif(i == eig(it)%nsec .or. nstates > 1) then
			jsecs = jsecs + 1
			write(10,'(i10,3x,1000f15.8)') nstates, eig(it)%esec(i), 
     . drates(i,i,1)/nstates, drates(i,i,2)/nstates,
     . 2*maxval(drates(i,1:i-1,1)/nstates), 
     . sum(drates(i,1:i,1:2))/nstates
		endif
	
	end do
	!write(10,*)
	!write(10,*)
	!write(*,*)'jsec, nsec = ', jsec, eig(it)%nsec

	if(fullout) then
	close(100)
	close(101)
	close(102)
	endif
	close(10)

	! for postprocessing by mathematica: could be in degen-sec-amp2.out
	open(11,file='degen-sec-amp2-header.dat',
     .               action='write',position='append')
	write(11,*) jsecs ! jsec
	close(11)

	
	deallocate(drates) ! so re-init=0 is done with new allocation
	return
	end subroutine  writedrates
!======================================================================

	
	end program
