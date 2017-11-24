cc
! NOtes: Assume g, w0 same for all sites
! if not, DPhicreat/annihil, D/Phi hops, all would need seperate sets of eigenstates/values to compute the amplitudes making the comput cot very high.
! might do this in future, but at the moment, PermSym=F can only occur for reasons that do not require diff eig.



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
	use modq, only: SetEcoul
	use mpi
	implicit none

	integer:: i,ic,is,ia,iter,ih,j,ntot,zt,icl
	integer:: nev,ncv,it,ier
	integer :: ntrap,stathc(34,4), stathc0(34,4)
	integer :: totcharge,charge
	double precision:: tottime, dtau
	integer:: trapiter,s,iss,nc,ns,itraj,nex,nelec,idw,ielec
	logical :: found,alloc
	integer :: x(3)
	integer :: nms(2), wayss(34)
	double precision:: rout(10)
	! mpi
	integer:: ierr, ntp, num_procs, node, maxtraj
	double precision, allocatable, dimension (:):: Iav, Iall
	character*50 :: frmt




	
	! time stamp & random seed 
	!call timestamp

	! readinput file
	call input()
	alloc = .false.;
!======================================================================
	call MPI_INIT(ierr)
	!find out MY process ID, and how many processes were started.
	call MPI_COMM_RANK (MPI_COMM_WORLD, node, ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
	! find ntp
	call setntp(ntraj,num_procs,ntp)
	maxtraj = ntp * num_procs; ! could be greater than ntraj
	allocate(Iav(ntp))

	!write(*,*) "Node = ",node," num_procs, ntp = ",num_procs, ntp 
	
	if(node==0) then
		call timestamp()
		allocate(Iall(maxtraj))
		write(frmt,'("(",I6.6,"G18.10)")') maxtraj+3 ! do we really need full out?
	endif
	
	!call printnode(node)

	! SetBondLengths() once for all loops or once for a trajectory?
	! currently, it's for a trajectory 

	do nex=mexmin,mexmax,dmex ! excitations
		do idw=1,ndw ! detuning
			dw = dwmin + (idw-1)*ddw;
			ielec=0;
			call setdetuning(dw,detuning,node)
			do nelec = nelmin,nelmax,dnelec
				if (nelec == 0 .and. onlydoped) cycle
				ielec=	ielec+1;
				do ier=1,ner ! Er: electric field * r_nns
					Er = Ermin + (ier-1)*dEr;
					if (VRH) Efieldh = signEr*Er ! VRH= variable range hopping
					! runs ntp trajectories, output => current: Iav, h/c counts: stathc
					call trajectory(ntp, Iav, stathc)
					! gather all data
					call mpi_gather(Iav, ntp, mpi_double_precision, Iall, ntp, 
     .                mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
					call CheckError(ierr, node, 'gather')
					! sum up all stathc
					call mpi_reduce(stathc,stathc0,34*4, mpi_integer,mpi_sum,
     .                                      0, MPI_COMM_WORLD, ierr)
					call CheckError(ierr, node, 'reduce')		
					if(node==0) then ! write output 
						call WriteIvsEr(maxtraj,Iall)
						call WriteStat(stathc0)
					endif
				enddo ! Er
			enddo ! ielectron
		enddo ! idw
	enddo ! nex
	if(node==0) then
		write(*,*)"TransOC: everything done.... " 
		call timestamp()
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
	subroutine WriteIvsEr(maxtraj,Iall)
	implicit none
	integer, intent(in):: maxtraj
	double precision, dimension(maxtraj), intent(in) :: Iall
	open(300,file="I-vs-Er.out",position='append')
	if(ier==1) then
		write(300,'(a,5x,4I10)') "# ier, nx, idw, iel =  ", ier, 
     .         nex, idw, ielec
		write(300,'("     ")')
		write(300,'("     ")')
	endif
	! statistics in modmain: statistics = (/ average, variance^2 /)
	!write(300,frmt) Er, statistics(Iall,ntraj), Iall 
	write(300,'(3G18.10)') Er, statistics(Iall,ntraj)
	close(300)
	return
	end subroutine WriteIvsEr
!===============================================
	subroutine WriteStat(stathc)
	implicit none
	integer, dimension(34,4), intent(in):: stathc
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
	subroutine writeways()
	implicit none
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
	return
	end subroutine writeways
!======================================================================
! tranejctry runs ntp traj, output: Iav ==> their average currents
!======================================================================
	subroutine trajectory(ntp,Iav, stathc)
	implicit none
	integer, intent(in) :: ntp
	double precision, dimension(ntp), intent(out) :: Iav
	integer, dimension(34,4), intent(out) :: stathc

	stathc = 0;
	
	do itraj=1, ntp
		zt = 0; ntrap = 0;
		totcharge=0; tottime=0.0d0
		trapiter=0

		nms = 0; wayss = 0;
		rout = 0.0d0
		
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

	x = 0;
	!========== iterations over number of hops
	! main loop over number of hops asked
	do iter=1,niter

		if(1==0) then
			!write(*,'(a)') ". . . . . . . . . . . . "
			write(*,'(a,i10,a,i10,a,2i10)') " itraj = ",itraj,
     .               " iter = ",iter,"  N, m = ",na,nx
			!write(*,*)        
		endif


		!write(*,*) iter, na, nx
		!write(*,*) "in: Asites = ",Asites
		!write(*,*) "in: occ = ",sys%occ

		
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
				do ia=1,19
					allocate(qt(ia)%cs(4,sys%nsites))
				enddo
				! allocate space for rates
				do ih=1,34 ! testing... alloc 26 all 
					if(PermSym .and. (.not. vrh))then
						allocate(rate(ih)%rcs(4,1))
					else
						allocate(rate(ih)%rcs(4,sys%nsites))
					endif
				enddo
				alloc = .true.
		endif

		do ih=1,34
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

		!write(*,*)"sum(sys%occ) = ", sum(sys%occ)

	if (leads)
     . x = x + (/nelec,sum(sys%occ), Qnet /);


		! All available hops: 
		!	transition amplitudes and amp^2 for degenerate sectors
		! and allocate space for rates???
		s = 1;
		call AllHops()
		if(debug)write(*,*) "main:   AllHops done... "	

		! Coulomb's interaction
		! sets Ecoul(ih=1:34)%dEq(is=1:ns(ih)) and Etotq for use in rates
		call SetEcoul()
		
		! rates from am2 and energetic penalties
		call CalRates()
		if(debug)write(*,*) "main:   CalRates done... "	

		!write(*,*) "main: ===1==> ",rate(:)%r
		!if(na==nsites) write(*,*)"main: N,m, R = ", na,nx,sum(rate(:)%r)

		! select a hop based on rates
		ih = ihSelect() 

		if(1==0 .and. na==nsites) then
			write(*,*)"main: Rtot = ",sum(rate(:)%r)
			write(*,*)"main: R_L/R = ",rate(1:8)%r
			write(*,*)"main: R_U/D = ",rate(27:34)%r
			write(*,*)"main: Rkappa, Rgamma= ",rate(25:26)%r
		endif
		
		!write(*,*) "main: ===2==> ",rate(:)%r
		if(debug)write(*,*) "main:   ihSelect: ih = ",ih
		call icsSelect(ih,ic,is)
		if(debug)write(*,*) "main:   icsSelect: ic,is =  ",ic,is

		rout(1:8) = rout(1:8) + rate(1:8)%r 
		rout(9:10) = rout(9:10) + rate(25:26)%r 
		
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
		dtau=1.0d0/sum(rate(1:nproc)%r)
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
	
	enddo ! itraj

	return
	end subroutine trajectory
!======================================================================

	end program
