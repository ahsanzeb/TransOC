
	module init
	use modmain
	use maps
	use modways
	implicit none
	contains
!-----------------------------------------
!	initialise some global variables
!-----------------------------------------
	subroutine initialise(doping)
	implicit none
	integer, intent(in) :: doping ! +ve for holes, -ve for electrons
	! local
	integer :: i,ina, n0,n1,n2
	integer :: nelec ! number of electrons

	nelec = nsites - doping;

	if(nelec <= 0 .or. nelec >= 2*nsites) then
		write(*,*)"Error(init): Nelectron <= 0 or >= Nsites !"
		write(*,*)"Charge transport cannot occur!"
		stop
	endif


	!-----------------------------------------
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
	!------------------------------------------

	!------------------------------------------
	!	sys: nsites, occ; na, Asites
	!------------------------------------------
	call initOcc(nelec)
	!-----------------------------------------
	! initialise ways, mapb, mapt
	!-----------------------------------------
	call UpdateWays()
	call UpdateReqType ! ReqType

	!write(*,*) "ReqType = ",mapt%req

	! does not exist, set to false, for use in mapt/mapb
	Hg(:)%xst = .false.
	basis(:)%xst = .false.
	
	call UpdateMapT ! mapt%map, mapt%cal

	!write(*,*) "mapt = ",mapt%map


	call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
	call UpdateMapB

	!-----------------------------------------
	! calculate maphc, the map from hopping index to amplitude index,
	!	and channel to channel location, needed for channel 8 only,
	! but does not cost much to make map for all 26 x 4.
	call calmaphc()
	!-----------------------------------------


	!-----------------------------------------
	! set occupations
	






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
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl
     . /), (/ 34,4 /), order=(/2,1/) );
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
	end subroutine initialise
!-------------------------------------------------------------
!	sets dqc array, energetic changes for various hops/channels
!-------------------------------------------------------------
	subroutine SetDQC() ! w0, Ebr, Ebl, Er 
  !local
	double precision, dimension(4):: dEbulk
	double precision, dimension(18):: dEcont
	integer, dimension(34):: signEr
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
	! kappa, gamma: -w0 to reflect cange in reference for spectrum
	! otherwise the transition would just take the energy difference between
	! kappa: what would be the energy of emitted photon? ~ w0-wR
	! similarly, exciton loss will release energy to the lattice;
	! it will not cost energy so no energetic penalty supression etc
	signEr =(/
     .   1, -1, -1, 1, 1, -1, -1, 1,
     .   1, 1, -1, -1, -1, -1, 1, 1,
     .   -1, 1, -1, 1, 1, -1, 1, -1,
     .   0, 0,
     .   0, 0, 0, 0, 0, 0, 0, 0 /);

		dqc(:,:) = 0.0d0
		do i=1,34
			select case(i)
				case(1:4,27:30)
					dqc(i,:) = dEbulk(:);
				case(5,6,31,32)
					dqc(i,:) = dEbulk(:) - Exb
				case(7,8,33,34)
					dqc(i,:) = dEbulk(:) + Exb
				case(9:26)
					dqc(i,:) = dEcont(i - 8)
			end select
			dqc(i,:) = dqc(i,:) - signEr(i)*Er
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
	end subroutine initOcc
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
	end module init
