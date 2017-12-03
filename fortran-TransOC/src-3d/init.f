
	module init
	use modmain
	use maps
	use modways
	implicit none
	double precision, allocatable, dimension(:,:) :: GaussArray
	!rnns values grid; a0-sigma0 to a0+sigma0
	integer:: nmax = 10000 ! large enough number; 

	
	contains


!-----------------------------------------
	subroutine init0()
	implicit none

	if (allocated(basis)) deallocate(basis)
	if (allocated(eig)) deallocate(eig)
	if (allocated(Hg)) deallocate(Hg)
	allocate(basis(NBasisSets))
	allocate(eig(NHilbertSpaces))
	allocate(Hg(NHilbertSpaces))
	! will need to compare n,m with 
	!req n,m in an iteration, so dont start with 0
	basis(:)%n = -1;
	Hg(:)%n = -1; Hg(:)%m = -1; Hg(:)%m1 = -1;
	!write(*,*) "ReqType = ",mapt%req
	! does not exist, set to false, for use in mapt/mapb
	Hg(:)%xst = .false.
	basis(:)%xst = .false.

	!-----------------------------------------
	! calculate maphc, the map from hopping index to amplitude index,
	!	and channel to channel location, needed for channel 8 only,
	! but does not cost much to make map for all 26 x 4.
	call calmaphc()
	!-----------------------------------------

	!-----------------------------------------
	! sets hopping parameters for various hops, channels
	!-----------------------------------------
	!{th, tl, tlh, thl, JhR, JlR, JhL, JlL} = tpar;
	!th = tpar(1); tl = tpar(2);
	!tlh = tpar(3); thl = tpar(4);
	!JhR = tpar(5); JlR = tpar(6);
	!JhL = tpar(7); JlL = tpar(8); 
	ts = reshape( (/
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . JlR,0.0d0,0.0d0,0.0d0, JhR,0.0d0,0.0d0,0.0d0,
     . JlL,0.0d0,0.0d0,0.0d0, JhL,0.0d0,0.0d0,0.0d0,
     . JhR,0.0d0,0.0d0,0.0d0, JlR,0.0d0,0.0d0,0.0d0,
     . JhL,0.0d0,0.0d0,0.0d0, JlL,0.0d0,0.0d0,0.0d0,
     . JlR,0.0d0,0.0d0,0.0d0, JhR,0.0d0,0.0d0,0.0d0,
     . JhR,0.0d0,0.0d0,0.0d0, JlR,0.0d0,0.0d0,0.0d0,
     . JlL,0.0d0,0.0d0,0.0d0, JhL,0.0d0,0.0d0,0.0d0,
     . JhL,0.0d0,0.0d0,0.0d0, JlL,0.0d0,0.0d0,0.0d0,
     . kappa,0.0d0,0.0d0,0.0d0, gamma,0.0d0,0.0d0,0.0d0,
     . th, tl, tlh, thl, ! 3D: up/down
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl, ! impurity; ic=1 applies only; let's use th at the moment
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl
     . /), (/ 42,4 /), order=(/2,1/) );
	!-----------------------------------------
	!Block injecton?? 
	!-----------------------------------------
	!ideally should also minimise computation of
	!corresponding transition amplitudes!
	!EBlock electron blocking layer; on left
	!HBlock hole blocking layer; on right
	! jump # 16,21 prob->0
	if(EBlock) then
		ts(16, 1) = 0.0d0;
		ts(21, 1) = 0.0d0;
	endif
	! jump # 10,18 prob->0
	if(HBlock) then
		ts(10, 1) = 0.0d0;
		ts(18, 1) = 0.0d0;
	endif
	!-----------------------------------------

	!-----------------------------------------
	! set dqc global vaiable 
	!-----------------------------------------
	call SetDQC()
	!-----------------------------------------

	!-----------------------------------------
	! set SigndEQ global vaiable in modways
	!-----------------------------------------
	call SetSigndEQ()

	return
	end subroutine init0

!------------------------------------------
	


!-----------------------------------------
!	initialise some global variables
!-----------------------------------------
	subroutine initialise(doping)
	implicit none
	integer, intent(in) :: doping ! +ve for holes, -ve for electrons
	! local
	integer :: i,ina, n0,n1,n2
	integer :: nelec ! number of electrons

	if(.not. impurity .and. iabs(doping) > nsites) then
		write(*,*)"Error(init): Nelectron <= 0 or >= 2*Nsites!"
		write(*,*)"Charge transport cannot occur!"
		stop	
	elseif(impurity .and. iabs(doping) > nsites-1) then
		write(*,*)"Error(init): Nelectron <= 0 or >= 2*Nsites!"
		stop	
	endif

	if(impurity) then		
		nelec = (nsites - 1) - doping; ! 1 site for impurity
	else
		nelec = nsites - doping;
	endif

	!------------------------------------------
	!	sys: nsites, occ; na, Asites
	!------------------------------------------
	call initOcc(nelec)
	!-----------------------------------------
	! initialise ways, mapb, mapt
	!-----------------------------------------
	call UpdateWays()
	call UpdateReqType ! ReqType
	call UpdateMapT ! mapt%map, mapt%cal
	call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
	call UpdateMapB
	!-----------------------------------------
	! set BondLengths global vaiable in modmain
	! considered positional disorder, and use vrh model for
	! prefactor of bare hopping rate vs dnns
	!-----------------------------------------
	!if (vrh) call SetBondLengths()
	call SetLattice(doping)

	return
	end subroutine initialise
!-------------------------------------------------------------
!	sets dqc array, energetic changes for various hops/channels
!-------------------------------------------------------------
	subroutine SetDQC() ! w0, Ebr, Ebl, Er 
  !local
	double precision, dimension(4):: dEbulk
	double precision, dimension(18):: dEcont
	double precision, dimension(8):: dEimp

	integer :: i,j

	! Exciton binding energy: Exb; Exciton energy = w0;
	! (LUMO/D energy = w0 + Exb; and Contact barriers are defind w.r.t this level)
  	! Important NOTE: 
  	! adjustig the 'energy of a D/Lumo filled' for
  	! the change in the energy reference of the eigenstates,
  	! we only need to consider the rest of the changes.
	dEbulk = (/ 0.0d0, 0.0d0, -w0, w0 /);
	dEcont = (/
     .   -Ebr, -Ebr + w0, -Ebl, -Ebl + w0,
     .   Ebr-Exb-w0, Ebr-Exb, Ebl-Exb-w0, Ebl-Exb,
     .   Ebr, -Ebr + w0+Exb, Ebr - w0, -Ebr+Exb,
     .   Ebl, -Ebl + w0+Exb, Ebl - w0, -Ebl+Exb,
     .   -w0, -w0 /);
	!if(impurity) then
	dEimp = (/ Eimp-w0-Exb, Eimp-Exb, -Eimp, w0-Eimp,
     .         w0+Exb-Eimp, Eimp, Exb-Eimp, Eimp-w0 /);
	!endif

	! kappa, gamma: -w0 to reflect cange in reference for spectrum
	! otherwise the transition would just take the energy difference between
	! kappa: what would be the energy of emitted photon? ~ w0-wR
	! similarly, exciton loss will release energy to the lattice;
	! it will not cost energy so no energetic penalty supression etc
		dqc(:,:) = 0.0d0
		do i=1,nproc
			select case(i)
				case(1:4,27:30)
					dqc(i,:) = dEbulk(:);
				case(5,6,31,32)
					dqc(i,:) = dEbulk(:) - Exb
				case(7,8,33,34)
					dqc(i,:) = dEbulk(:) + Exb
				case(9:26)
					dqc(i,:) = dEcont(i - 8)
				case(35:42)
					dqc(i,:) = dEimp(i - 34)
			end select
			!write(*,'(a,i5,a,4f15.10)')"ih = ",i,"  dqc = ",dqc(i,:)
		end do

	return
	end subroutine SetDQC
!-----------------------------------------
!------------------------------------------
!	sys: nsites, occ; na, Asites
!------------------------------------------



	subroutine initOcc(nelec)
	implicit none
	integer, intent(in) :: nelec ! number of electrons
	! local
	integer :: i,ina, n0,n1,n2,maxocc
	
	sys%nsites = nsites
	if (allocated(sys%occ)) deallocate(sys%occ)
	allocate(sys%occ(nsites))

	!write(*,*)" init: Nsites, nelec = ",nsites, nelec

	ina = 0;
	sys%occ(:) = 0;
	maxocc = 2;
	if (impurity) 	sys%occ(4) = impocc;
	
	!write(*,*) "init: mincarriers = ",mincarriers

	if(mincarriers) then
		if(nelec < nsites) then
			maxocc = 1;
			!write(*,*) "init: nelec < nsites ? ",nelec, nsites
		else ! nelec >= nsites
			if(impurity) then
				sys%occ(:) = 1;
				sys%occ(4) = impocc;
				ina = nsites-1 + impocc;
			else
				sys%occ(:) = 1;
				ina = nsites
			endif	
		endif
	endif

	!write(*,*) "init: ina = ",ina
	if(EqualDistr) then ! Equally distribute the carriers/dopants
		!	so that all three chains can get more or less equal carriers
		!	this option will be useful for OneDChains=T case
		i = 0;
		do while (ina < nelec)
			i = i + 1;
			if (	sys%occ(i) < maxocc) then
				if((.not. impurity) .or. (i .ne. 4)) then
					sys%occ(i) = sys%occ(i) + 1;
					ina = ina + 1;
				endif
			endif	
		enddo
	else ! random, with occ <= maxocc
		do while (ina < nelec)
			i = int(1+nsites*rand(0));
			!write(*,*)"i  = ",i
			if (	sys%occ(i) < maxocc ) then
				if((.not. impurity) .or. (i .ne. 4)) then
					sys%occ(i) = sys%occ(i) + 1;
					ina = ina + 1;
				endif
			endif	
		enddo
	endif

	n0 = 0;n1=0;n2=0;
	do i=1,nsites
			if(sys%occ(i)==0) then
				n0 = n0 +1
			elseif(sys%occ(i)==1) then
				if((.not. impurity) .or. (i .ne. 4)) then
					n1 = n1 + 1
				endif
			elseif(sys%occ(i)==2) then
				n2 = n2 +1
			else
				write(*,*) "Error(init): sys%occ(i)>2 ? "
				stop 
			endif
	enddo

	! set Asites
	if (allocated(Asites)) deallocate(Asites)
	allocate(Asites(n1))
	ina = 1;
	do i=1,nsites
		if (sys%occ(i)==1) then
			if(.not. impurity .or. i .ne. 4) then
				Asites(ina) = i;
				ina = ina + 1;
			endif
		endif
	enddo

	! set global var
	sys%n0=n0;
	sys%n1=n1;
	sys%n2=n2;
	if(impurity) then
		sys%nimp=1;
	else
		sys%nimp=0;
	endif

	na = n1;
	!write(*,*)"init: na = ",na

	return
	end subroutine initOcc




!---------------------------------
	subroutine initOcc0(nelec)
	implicit none
	integer, intent(in) :: nelec ! number of electrons
	! local
	integer :: i,ina, n0,n1,n2,maxocc

	sys%nsites = nsites
	if (allocated(sys%occ)) deallocate(sys%occ)
	allocate(sys%occ(nsites))

	!write(*,*)" init: Nsites, nelec = ",nsites, nelec

	ina = 0;
	sys%occ(:) = 0;
	maxocc = 2;
	!write(*,*) "init: mincarriers = ",mincarriers

	if(mincarriers) then
		if(nelec < nsites) then
			maxocc = 1;
			!write(*,*) "init: nelec < nsites ? ",nelec, nsites
		else ! nelec >= nsites
			sys%occ(:) = 1;
			ina = nsites;
			!write(*,*) "init: **** ina = ",ina
		endif
	endif

	!write(*,*) "init: ina = ",ina
	

	if(EqualDistr) then ! Equally distribute the carriers/dopants
		!	so that all three chains can get more or less equal carriers
		!	this option will be useful for OneDChains=T case
		i = 0;
		do while (ina < nelec)
			i = i + 1;
			if (	sys%occ(i) < maxocc ) then
				sys%occ(i) = sys%occ(i) + 1;
				ina = ina + 1;
			endif	
		enddo
	else ! random, with occ <= maxocc
		do while (ina < nelec)
			i = int(1+nsites*rand(0));
			!write(*,*)"i  = ",i
			if (	sys%occ(i) < maxocc ) then
				sys%occ(i) = sys%occ(i) + 1;
				ina = ina + 1;
			endif	
		enddo
	endif

	n0 = 0;n1=0;n2=0;
	do i=1,nsites
			if(sys%occ(i)==0) then
				n0 = n0 +1
			elseif(sys%occ(i)==1) then
				n1 = n1 + 1
			elseif(sys%occ(i)==2) then
				n2 = n2 +1
			else
				write(*,*) "Error(init): sys%occ(i)>2 ? "
				stop 
			endif
	enddo

	! set Asites
	if (allocated(Asites)) deallocate(Asites)
	allocate(Asites(n1))
	ina = 1;
	do i=1,nsites
		if (sys%occ(i)==1) then
			Asites(ina) = i;
			ina = ina + 1;
		endif
	enddo

	! set global var
	sys%n0=n0;
	sys%n1=n1;
	sys%n2=n2;

	na = n1;
	!write(*,*)"init: na = ",na

	return
	end subroutine initOcc0
!----------------------------------------
!	sign of charging energy for various contact processes
! i.e., electron/hole extraction/injection
	subroutine SetSigndEQ()
	implicit none
	integer:: ih

	do ih = 9,24
		select case(ih)
			case(9:12,18,20,22,24)
				! e extraction/hole injection
				SigndEQ(ih-8) = 1
			case(13:16,17,19,21,23)
				! e injection/hole extraction
				SigndEQ(ih-8) = -1
		end select
	enddo

	return
	end 	subroutine SetSigndEQ
!----------------------------------------
! bond lengths from a distribution
!	stored in the array BondLengths.
! first dim for site index, second for an index given based on
! which nns is the second site in the bond.
!	basis idea: for a site is1, the second site is2 can be in the range is1-3: is1+3
!	shift this by -is1+4 to move it to 1:7, four out of these 7 
!	correspond to the nns of is1.
! To read dij use ReadBondLengths() in rates module
	subroutine SetBondLengths() 
	implicit none
	integer:: i, is1, is2, is2s(4), ds1(4),ds2(4),ds3(4), id1,id2
	double precision:: x

	! sets GaussArray, accumulated probabilities vs rnns
	if(allocated(GaussArray)) deallocate(GaussArray)
	allocate(GaussArray(nmax+1,2))
	call SetGaussArray(nmax)
	
	ds1 = (/ -3,3,-2,-1/);
	ds2 = (/ -3,3, 1, 2/);
	ds3 = (/ -3,3, 1,-1/);
	
	if(allocated(BondLengths)) deallocate(BondLengths)
	allocate(BondLengths(nsites,7))
	
	BondLengths = -1.0d2
		do is1=1,nsites
				if(mod(is1,3)==0) then ! 3,6,9,....
					is2s = is1 + ds1
				elseif(mod(is1-1,3)==0) then ! 1,4,7,....
					is2s = is1 + ds2
				else !if(mod(is1-2,3)==0) then ! 2,5,8,....
					is2s = is1 + ds3
				endif
			do i=1,4
				is2 = is2s(i); 
				! for is at ends triangles, is2 can be -ve or > nsites,
				! but it does not matter here for the purpose
				! of setting an index of bonds. 
				! But, be careful when retrieving dij if {is1,is2} are in {first,last} triangles, any order.
				!	in such a case, set one of the end triangles to a shifted 
				!	ficticious triangles below/above physical lattice and use
				! the is2 index of site in this ficticious triangle to find id1, and read dij.
				id1 = is1 - is2 + 4
				id2 = -id1 + 8; !is2 - is1 + 4
				x = dnns() ! calculated nns distance
				BondLengths(is1,id1) = x;
				if(is2 > 0 .and. is2 <= nsites) then
					BondLengths(is2,id2) = x; ! maybe overwritten once, but that's ok
					! as long as a given bond has a unique bondlength when read in
					! reference to any of the two connecting sites is1/is2!
				endif
			enddo		
		enddo
	! To read bond lengths:
	!if(periodic) then
		! when retrieving:
		! convention: {is1,is2} are in {first,last} triangles:
		!	is2 <==== is2-nsites to get left virtual triangle
	!else! leads
		! when retrieving:
		! assume virtual triangle with -2,-1,0 sites as left contact,
		! and that with n+1,n+2,n+3 as the right contact
		! first gives the bonds with sites 1,2,3; second with n-2,n-1,n
	!endif

	return
	end subroutine SetBondLengths
!----------------------------------------
!	ReadBondLengths not used now, code moved to moved to dp function dij(ih,is);
!----------------------------------------
! uses GaussArray calculated by SetGaussArray
	double precision function dnns()
	implicit none
	logical :: found
	double precision:: eta,p1,p2
	integer:: ir
	
	eta = rand(0) * GaussArray(nmax+1,2);
	! location of eta on accumulated probabilities
	found = .false.;
	do ir=1,nmax
		p1=GaussArray(ir,2); p2=GaussArray(ir+1,2);
		if (eta .ge. p1 .and. eta .lt. p2 ) then
			dnns = GaussArray(ir,1);
			found = .true.
			exit
		endif
	end do	
	! error?
	if (.not. found) then
		write(*,*) "dnns: something wrong...."
		stop
	endif

	!open(123,file='dnns.out',action='write',position='append')
	!write(123,*) dnns
	!close(123)
	
	return
	end function dnns
!----------------------------------------
	! a gaussian distribution; accumulated 
	subroutine SetGaussArray(nmax)!,a0,sigma0)
	implicit none
	integer, intent(in)::nmax
	!double precision, intent(in) :: a0, sigma0
	!double precision, dimension(nmax+1,2):: GaussArray
	double precision:: amin, da, r, pref, psum, tsigma02
	integer :: i

	tsigma02 = 2*sigma0**2;
	!pref = 1/dsqrt(3.1415927*tsigma02);
	
	amin=a0-nsigma*sigma0;
	da = 2*nsigma*sigma0/(nmax-1)

	GaussArray = 0.0d0; psum = 0.0d0;
	do i=1,nmax
		r = amin + (i-1)*da;
		GaussArray(i,1) = r;
		psum = psum + dexp(-(r-a0)**2/tsigma02); !pref*dexp(-(r-a0)**2/tsigma02);
		GaussArray(i+1,2) = psum;
	enddo

	return
	end 	subroutine SetGaussArray
!---------------------------------------
	subroutine SetLattice(doping)
	implicit none
	integer, intent(in) :: doping
	integer :: i,is,it,sgn
	double precision:: x,y,z

	if(.not. allocated(sys%r)) allocate(sys%r(-2:nsites+3,3))
	! redular lattice: a chain of equilateral trianlges
	sys%r(-2,:) = (/ 0.0d0, 0.0d0, 0.0d0 /); 
	sys%r(-1,:) = (/ 0.0d0, a0, 0.0d0 /); 
	sys%r( 0,:) = (/ 0.0d0, 0.5d0*a0, 0.866d0*a0 /);
	is = 0;
	do it = 1, nsites/3 + 1;
		do i=1,3
			is = is + 1;
			sys%r(is,:) = sys%r(is-3,:) + (/ a0,0.d0,0.d0 /);
			!write(*,'(a,i3,2x,3f10.5)')"init: is, sys%r = ",is,sys%r(is,:)
		enddo
	enddo

	! positonal disorder?
	if(vrh)then
		! sets GaussArray, accumulated probabilities vs rnns
		if(.not. allocated(GaussArray)) allocate(GaussArray(nmax+1,2))
		call SetGaussArray(nmax)

		!add disorder
		do is=1,nsites
			! dnns(): random selection from a gaussian distribution
			! peaked around a0
			x = dnns() - a0; 
			y = dnns() - a0;
			z = dnns() - a0;
			sys%r(3+is,:) = sys%r(3+is,:) + (/x,y,z/)/1.73d0; !sqrt(3)?
			!write(*,*)"init: disor dr = ",dsqrt(x**2+y**2+z**2)
		enddo
	endif

	! donor/acceptor positons? only a single type of doping at the moment
	! q0 ====> charged system?
	if(coulomb) then
		!if(allocated(sys%q0)) deallocate(sys%q0)
		!allocate(sys%q0(nsites))
		! also allocate sys%q for later use
		if(allocated(sys%q)) deallocate(sys%q)
		allocate(sys%q(nsites))
		!sys%q0 = 0.0d0 ! D/Phi ===> q = -1/+1; q0 = 0 for all
		!sgn = isign(1,doping);
		!do is = 1, doping*sgn
		!	! choose sites at random for donor/acceptor
		!	i = int(rand(0)*nsites) + 1;
		!	sys%q0(i) = -sgn*1.0d0; ! accepter is -ve charged, donor is +ve charged
		!enddo
	endif
	return
	end 	subroutine SetLattice
!---------------------------------------
	
	end module init
