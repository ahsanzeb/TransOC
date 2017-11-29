
	module rates
	use modmain, only: qt,eig,itypes, 
     .  mapb,mapt,maph,mapc,nokappa, nogamma
	implicit none

	public:: CalRates
	private:: totr, ratehcs, PenaltyArray !, Penalty

	contains
!----------------------------------------------------
! calculates total rates for all processes
!----------------------------------------------------
	subroutine CalRates()
	use modmain, only: crosshops,nolosses,leads,impurity
	implicit none
	! local
	integer :: nc

	!write(*,*) "rates  0"
	! nc: only for ih=1-8
	! no crosshops? nc=1-2, otherwise nc=4
	nc=2;
	if (crosshops) nc=4; 
	! calculate rate(:)%r and rate(:)%rcs(:,:) for all hops
	call totr('dphihops',nc)
	!write(*,*) "rates  1: "
	
	call totr('annihilation',nc)
	!write(*,*) "rates  2"
	call totr('creation',nc)
	!write(*,*) "rates  3"
	if (.not. nolosses) call totr('losses',1)
	if (leads) call totr('contacts',1)
	if (impurity) call totr('impurity',1)

	return
	end subroutine CalRates

!******************************************************
	subroutine totr(process,nc)
	! calculates c,s resolved rate and total rate
	!	rate(ih)%rcs(ic,is) and rate(ih)%r
	use modmain, only: rate, ways, PermSym,na,nx,ts,
     .              ihdphops,ihcreat,ihannih,vrh,coulomb
	implicit none
	character(len=*), intent(in):: process
	integer, intent(in) :: nc
	! local
	integer:: ih, ic, ih1, ih2,is,iss
	integer, dimension(8):: ihs
	logical :: nvrh, psnvrh
	
	nvrh  = (.not. vrh) .and. (.not. coulomb); ! no variable range hopping?
	psnvrh = PermSym .and. nvrh;
	! if PermSym, amp for is=1 only;
	! if vrh, is>1 is given in ratehcs, fix it there if PermSym.
	
	!write(*,*)"totr: ways(:)%ns = ", ways(:)%ns

	select case(process)
	case('dphihops','creation')
	!-------------------------------
	! ih:1-4, 27-30 DHopsR/L, PhiHopsR/L, DHops Up/Down, PhiHops Up/Down
	!-------------------------------
	! ih=7,8; 33,34 (D,Phi), (Phi,D) created : R/L; Up/Down
	!--------------------------------------------------------

	! If PermSym; rates for 27=28; 29=30; 33=34; avoid repeat calc
		if(process == 'dphihops') then
			ih2=8; ihs = ihdphops;
		else ! 'creation'
			ih2=4; ihs(1:4) = ihcreat;
		endif
		do ih1=1,ih2
			ih = ihs(ih1)
			rate(ih)%rcs(:,:) = 0.0d0
			do is=1,ways(ih)%ns,1
				do ic=1,nc
						!write(*,*)"** ih,is,ic=...**",ih,is,ic
					! PermSym: use equiv ih's if rates already calculated;
					!	ih ordered in 'ihs'
					! equiv because no E.r term for in-plane hops ===> Up==Down
					if(PermSym) then 
						! check if we already have rates for an equivalent hop?
						if((ih==28 .and. ways(27)%ns > 0) .or.
     .         (ih==30 .and. ways(29)%ns > 0) .or.
     .         (ih==34 .and. ways(33)%ns > 0) ) then
							rate(ih)%rcs(ic,is) = rate(ih-1)%rcs(ic,is); ! 27 for 28, etc.; valid even if vrh
						elseif(nvrh)then
							call ratehcs(ih,ic,is)
							rate(ih)%rcs(ic,is)=
     .             ts(ih,ic)*rate(ih)%rcs(ic,is)*ways(ih)%ns
						else ! vrh: calculate rates for all sites is=1,ns
							call ratehcs(ih,ic,is)
							rate(ih)%rcs(ic,is) = ts(ih,ic)*rate(ih)%rcs(ic,is)
						endif
					else
						call ratehcs(ih,ic,is)
						rate(ih)%rcs(ic,is)=ts(ih,ic)*rate(ih)%rcs(ic,is)
					endif
				end do ! ic
				if(psnvrh) exit; ! only a single site/case for each hop type
			end do ! is
			rate(ih)%r = sum(rate(ih)%rcs); ! total rate for ih hop
			!write(*,*) "ih, rate(ih)%r = ",ih,rate(ih)%r
		end do ! ih

	case('annihilation')
	!-------------------------------------------------
	! ih:5-6 D,Phi annihilation R,L ==> D,Phi | Phi,D
	! ih:31,32 D,Phi annihilation up,dn ==> D,Phi | Phi,D
	!-------------------------------------------------
		do ih1=1,4
			ih = ihannih(ih1)
			rate(ih)%rcs(:,:) = 0.0d0
			rate(ih)%r = 0.0d0;
			iss=1; ! annihilation, amp available for is=1 only
			do is=1,ways(ih)%ns
				do ic=1,nc
					call ratehcs(ih,ic,iss)
					if(psnvrh) then 
						! total rate = (rates for one site) * ns
						rate(ih)%rcs(ic,is)=
     .           ts(ih,ic)*rate(ih)%rcs(ic,is)*ways(ih)%ns
     			else ! no PermSym or vrh
						rate(ih)%rcs(ic,is) = ts(ih,ic)*rate(ih)%rcs(ic,is)
     			endif
				enddo! ic
				if(psnvrh) exit; ! only a single site/case for each hop type
			enddo ! is
			rate(ih)%r = sum(rate(ih)%rcs(:,:)); ! total rate for ih hop
		end do
	case('losses')
	!-------------------------------------------------
	! ih:25,26 kappa, gamma losses
	!-------------------------------------------------
		! caivty losses
		ih=25; ic=1; is=1
		rate(ih)%rcs(:,:) = 0.0d0
		if ((.not. nokappa) .and. nx > 0) then
			call ratehcs(ih,ic,is)
			rate(ih)%r = ts(ih,ic)*rate(ih)%rcs(ic,is); ! total rate for ih hop
			!write(*,*) "rate:  kappa ==> ",rate(ih)%r
		else
			rate(ih)%r = 0.0d0
		endif
		!write(*,*) "rate:  kappa ==> ",rate(ih)%r

		! exciton losses
		ih=26; ic=1;
		rate(ih)%rcs(:,:) = 0.0d0
		if ((.not. nogamma) .and. nx*na > 0 ) then ! both na, nx >0
			do is=1,na
				call ratehcs(ih,ic,is)
				if(PermSym) then ! total rate = (rates for one site) * ns
					rate(ih)%rcs(ic,is)=
     .       ts(ih,ic)*rate(ih)%rcs(ic,is)*na
     			exit; ! only a single site/case for each hop type; no ic loop
				else
					rate(ih)%rcs(ic,is)=ts(ih,ic)*rate(ih)%rcs(ic,is)
				endif
			end do
			rate(ih)%r = sum(rate(ih)%rcs(:,:)); ! total rate for ih hop
			!write(*,*) "rate:  gamma ==> ",rate(ih)%r
		else
			rate(ih)%r = 0.0d0
		endif
		!write(*,*) "rate:  gamma ==> ",rate(ih)%r

	case('contacts')
	!-------------------------------------------------
	! ih:9-24; R/L contacts injection/extraction of e/h
	!-------------------------------------------------
		ic=1;
		do ih=9,24
			rate(ih)%rcs(:,:) = 0.0d0
			do is=1,ways(ih)%ns
				if(ih <= 16) then ! annihilation, only qt%cs(ic,is=1)%amp so set site index=1
					iss=1;
				else
					iss = is;
				endif
				if(ih==19 .or. ih==20 .or. ih==23 .or. ih==24) then
					if (nx > 0) call ratehcs(ih,ic,iss)
				else
					call ratehcs(ih,ic,iss)
				endif
				if(psnvrh) then ! total rate = (rates for one site) * ns
					rate(ih)%rcs(ic,is)=
     .       ts(ih,ic)*rate(ih)%rcs(ic,is)*ways(ih)%ns
					exit; ! only a single site/case for each hop type
				else
					rate(ih)%rcs(ic,is)=ts(ih,ic)*rate(ih)%rcs(ic,is)
				endif
			enddo ! is		
			rate(ih)%r = sum(rate(ih)%rcs); ! total rate for ih hop
		enddo ! ih

	case('impurity')
	!-------------------------------------------------
	! ih:34-42; impurity level: injection/extraction of e/h
	!-------------------------------------------------
		ic=1;
		do ih=34,42
			rate(ih)%rcs(:,:) = 0.0d0
			do is=1,ways(ih)%ns
				if(ih==41 .or. ih==42) then ! D,Phi creat, m-1 case
					if (nx > 0) call ratehcs(ih,ic,is)
				else
					call ratehcs(ih,ic,is)
				endif
				if(psnvrh) then ! total rate = (rates for one site) * ns
					rate(ih)%rcs(ic,is)=
     .       ts(ih,ic)*rate(ih)%rcs(ic,is)*ways(ih)%ns
					exit; ! only a single site/case for each hop type
				else
					rate(ih)%rcs(ic,is)=ts(ih,ic)*rate(ih)%rcs(ic,is)
				endif
			enddo ! is		
			rate(ih)%r = sum(rate(ih)%rcs); ! total rate for ih hop
		enddo ! ih

	end select
	!-------------------------------------------------
	return
	end subroutine totr
!******************************************************
	subroutine ratehcs(ih,ic,is)
	! rates for ih hop, ic channel, is site
	! sets global variable rate(ih)%rcs(ic,is)
	use modmain, only: nx,qt,eig,itypes,beta,
     .  mapt,maph,mapc,Einit,rate,ways,dqc,dEQs,
     .  simplepf,a0,dinvl,vrh,PermSym,Ecoul, Etotq, Er
	implicit none
	integer, intent(in) :: ih,ic,is
	! local
	double precision, allocatable :: de(:)
	integer :: nsec, ia,icl,itl,i,isa
	double precision:: x,y,rnns,tvr
	logical :: dblehop, emptyhop

	! I think this if block is not needed, the condition in it
	! is already filtered before calling this subroutine.
	! if no available hops, set rate = 0

	!write(*,*) 'heloowwww 1'

	if ((ih <= 8 .or. ih >= 27) .and. ways(ih)%ns == 0) then
		rate(ih)%rcs(ic,is) = 0.d0
		return
	endif

	!write(*,*) 'heloowwww 2'


	! conditions on nx for ih=1-4
	if(nx==0) then
						! Dhops Homo-Homo channel (ic=1) blocked, 
						! Lumo-homo (ic=3) not available because
						! all spin down: homo filled on all active sites
						! similarly, D,Phi creation only ic=4 possible
						! (conditions in ratehcs()
						! amplitudes for these ic's and ih's not calculated
		if(((ih==1 .or. ih==2) .and. (ic==1 .or. ic==3)) .or. 
     . ((ih==3 .or. ih==4) .and. (ic==2 .or. ic==3))) then
			rate(ih)%rcs(ic,is) = 0.0d0 
			return
		endif

		if(((ih==27 .or. ih==28) .and. (ic==1 .or. ic==3)) .or. 
     . ((ih==29 .or. ih==30) .and. (ic==2 .or. ic==3))) then
			rate(ih)%rcs(ic,is) = 0.0d0 
			return
		endif

	endif

	!write(*,*) 'heloowwww 3'

	! conditions on nx for ih=7,8; 33,34
	if (ih==7 .or. ih==8 .or. ih==33 .or. ih==34) then
		if ( (ic < 3 .and. nx .lt. 1) .or. 
     .   (ic == 3 .and. nx .lt. 2) ) then
			rate(ih)%rcs(ic,is) = 0.d0;
			!write(*,*)"rate: nx cond; ih,ic,is=",ih,ic,is
			return 
		endif
	endif

	isa = is; ! site index for amplitudes
	! annihilation? set isa=1 if PermSym in case vrh is used
	if (PermSym .and. is > 1) then
		! this happens when vrh=T
		isa = 1;
	else
		isa = is;
	endif

	!write(*,*) 'heloowwww d'

	! locations of data
	itl = mapt%map(itypes(ih,ic)) ! location of final hilber space
	ia = maph(ih,ic); ! location of amplitudes
	icl = mapc(ih,ic); ! localtion of ic; icl=ic except ih=8, ic 1,2 swapped

	!write(*,*)"rate -------------"
	!if(ih==26) then
	!	write(*,*)"amp2 alloc:",allocated(qt(ia)%cs(icl,is)%amp)
	!	write(*,*)"amp2 alloc:",allocated(qt(ia)%cs(icl,is)%amp2)
	!endif

	! no of degenerate sectors
	nsec = eig(itl)%nsec 
	allocate(de(nsec))

	!write(*,*) 'heloowwww e'

	! change in energy for all transitions
	de = eig(itl)%esec(1:nsec) + dqc(ih,ic); ! changes due barries, Exciton binding, reference, etc.
	de = de - Einit; ! total change in energy

	! Change in Coulomb's energy + Efield energy:
	!	Ecoul is calculated in main, routine SetEcoul() in modq/coulomb.f
	de = de + Ecoul(ih)%dEq(is);

	! charging energy contact hops:
	! A cheaper rough approx of true coulomb energy
	if(ih .ge. 9 .and. ih .le. 24) then
		de = de + dEQs(ih-8)
	endif

	!write(*,*)"rates: ih, de, dEq= ",ih,de(1),Ecoul(ih)%dEq(is)
	

	!write(*,*) "rate: ih, de = ", ih, de

	! rates; total for this hop/hchannel/site

	if(.not. allocated(qt(ia)%cs(icl,isa)%amp2)) then
	write(*,*) "alloc amp2=>F: ih,ic,is,icl,ia",ih,ic,is,icl,ia
	stop
	endif

	! variable range hopping: 
	if(VRH) then
		! find 	nns distance for this hop
		!rnns = dij(ih,is);
		!rnns = dijl(ih,is);
		rnns = ways(ih)%rij(is); ! PermSym not used so is is true index of hop
		tvr = dexp(-dinvl*(rnns-a0)/a0); ! bare amplitude prefactor
	else
		tvr = 1.0d0;
	endif

	!write(*,*) "ih, Er, dEnet = ", ih, Er, de(1:min(5,nsec))

	! Energy of leacking photon ?????? 
	! ih=25; photon leackage: w_photon = |dE| if dE < 0;
	! so set pf = 1 for all negative dE, and Boltzmann factor for +ve.
	if(simplepf .or. ih==25) then 
		! simple penalty function: min(1,Exp[-beta*dE])
		! ih=25 because the de does not contain energy of leaking photon
		! which should make de
		rate(ih)%rcs(ic,is) =
     .  tvr * sum(PenaltyArray(de,nsec) * qt(ia)%cs(icl,isa)%amp2)
	else ! use bath spectral density to evaluate the penalty function
		rate(ih)%rcs(ic,is) =
     .  tvr * sum(PenaltyArray2(de,nsec) * qt(ia)%cs(icl,isa)%amp2)
	endif

	!write(*,*) '=========================='
	!write(*,*) 'PenaltyArray(de,nsec) = ',PenaltyArray(de,nsec)
	!write(*,*) 'PenaltyArray2(de,nsec) = ',PenaltyArray2(de,nsec)

	
	! for debugging, remove later
	if (isnan(rate(ih)%rcs(ic,is)) .or. 1==0) then
		write(*,*)"rates: rcs =",rate(ih)%rcs(ic,is)
		x = 0.0d0
		do i=1,nsec
			y = de(i)*beta;
			if (y < 0.0d0) then
				x = x + qt(ia)%cs(icl,is)%amp2(i)
			elseif(y < 20.0d0)then
				x = x + dexp(-y) * qt(ia)%cs(icl,is)%amp2(i)
			endif
		enddo
		write(*,*)"******************************"
		write(*,*)"rates: x = ",x
		write(*,*)"rates: ih,ic,is = ",ih,ic,is
		write(*,*)"rates: Einit=",Einit
		write(*,*)"rates: eig(itl)%esec(1)=",eig(itl)%esec(:)
		write(*,*)"rates: dqc(ih,ic)=",dqc(ih,ic)
		write(*,*)"rates: de=",de(1)			
		write(*,*)"PenaltyArray(de,nsec)=",PenaltyArray(de,nsec)
		write(*,*)"qt(ia)%cs(icl,is)%amp2=",qt(ia)%cs(icl,is)%amp2
		stop
	endif

	!write(*,*) "rates:********************"
	!write(*,*) "ih,ic, dqc(ih,ic) = ",ih,ic, dqc(ih,ic)
	!write(*,'(i5,5x,f12.6)') ic, rate(ih)%rcs(ic,is) 

	!	the final degenerate sector will be slected on basis of amp2
	!	once ih, ic,is are selected on basis of rate(ih)%rcs, and rate(ih)%r

	return
	end subroutine ratehcs
!******************************************************
	function PenaltyArray(de,ne)
	use modmain, only: beta, wcut !, J0 = 1 
	implicit none
		integer,intent(in) :: ne
		double precision, dimension(ne),intent(in):: de
		double precision, dimension(ne):: tmp
		double precision, dimension(ne):: PenaltyArray
		! local
		integer:: i
		double precision:: x, fac=40.0d0

		! Ohmic spectral density with high energy cutoff wcut:
		! J(w) = J0 * w * dexp(-(w/wcut)**2);
		!Sw = Jw(w) * [ nw(w) + 1 ]; w >= 0 ===> de < 0 here
		PenaltyArray = 1.d0;
		tmp = de*beta;
		do i=1,ne
			if(tmp(i) > fac) then
				PenaltyArray(i)=0.0d0
			elseif (tmp(i) > 0.0d0) then
				PenaltyArray(i)= dexp(-tmp(i))
			endif
		end do

	return
	end function PenaltyArray
!******************************************************
	function PenaltyArray2(de,ne)
	use modmain, only: beta, wcut !, J0 = 1 
	implicit none
		integer,intent(in) :: ne
		double precision, dimension(ne),intent(in):: de
		double precision, dimension(ne):: tmp
		double precision, dimension(ne):: PenaltyArray2
		! local
		integer:: i
		double precision:: x, fac=80.0d0, smalle=1.0d-3

		! Ohmic spectral density with high energy cutoff wcut:
		! J(w) = J0 * w * dexp(-(w/wcut)**2);
		!Sw = Jw(w) * [ nw(w) + 1 ]; w >= 0 ===> de < 0 here
		PenaltyArray2 = 1.d0;
		tmp = de*beta;
		do i=1,ne
			if(abs(tmp(i)) > fac) then
				PenaltyArray2(i)=0.0d0
			else
				x = dexp(-(de(i)/wcut)**2);
				if (abs(tmp(i)) <= smalle) then ! avoid 0/0
					PenaltyArray2(i) = x / (1.d0 + 0.5*tmp(i)) 
				elseif(tmp(i) > smalle) then ! positive de
					PenaltyArray2(i) = x * tmp(i)/(dexp(tmp(i)) - 1.0d0)
				else ! negative de
					PenaltyArray2(i)=-x * tmp(i)/(1.0d0-dexp(tmp(i)))
				endif
			endif
		end do

	return
	end function PenaltyArray2
!******************************************************
!	finds the inter-site distances
	double precision function dij(ih,is)
	use modmain, only: BondLengths, periodic, leads, nsites,ways
	implicit none
	integer, intent(in):: ih,is
	integer :: l1, l2
	integer::is1,is2,id
	!logical :: log1, log2

	l1 = ways(ih)%active(is); ! main site associated with the hop
	if(ih >= 9 .and. ih <= 24) then ! contact hops
		l2 = -100
	else
		l2 = ways(ih)%sites(is); ! the second site involved
	endif

	!----------------------------------------
	! call ReadBondLengths(l1,l2) 
	!----------------------------------------
	! is1 is always in 1:nsites
	! is2 can be <0 or >nsites; diff interpretation for periodic vs leads
	if (leads) then
		! when retrieving bondlengths:
		! assume virtual triangle with -2,-1,0 sites as left contact,
		! and that with n+1,n+2,n+3 as the right contact
		! first gives the bonds with sites 1,2,3; second with n-2,n-1,n
		if (l2 == -100) then ! contact hops: only possible when leads=T and l1 is nns of a lead/contact
			if(l1 <= 3) then ! left contact
				is1 = l1; is2 = is1 - 3; ! left virtual trianlge
			elseif(l1 >= nsites-2) then ! right contact
				is1 = l1; is2 = is1 + 3 ! right virtual trianlge
			endif
		else ! bulk hops;
			is1 = l1; is2 = l2; 
		endif
	elseif (periodic) then
		! when retrieving bondlengths:
		! convention: if {is1,is2} are in {first,last} triangles, use left virtual triangles for bonds.

		!i.e., if is1,is2 are in last phys triangle of sites, and next virtual triangle,
		! take it the same as left vitual trianlge and first phys trianlge	
		!	if is2> nsites, then is2 = is2-nsites to get left virtual triangle
		
		! if the bonds are between first and the last triangles?
		!log1 = l1 <= 3 .and. l2 >= nsites-2;
		!log2 = l2 <= 3 .and. l1 >= nsites-2;
		! if log1/log2; left virtual trianlges for bonds
		if( l1 <= 3 .and. l2 >= nsites-2) then
			is1 = l1; is2 = l2 - nsites;
		elseif(l2 <= 3 .and. l1 >= nsites-2) then
			is1 = l2; is2 = l1 - nsites;
		else
			is1 = l1; is2 = l2; 
		endif
	else !  onlybulk=T ??
		write(*,*)"Error(dij): onlybulk=T? not implemented",
     . " (only this dij case?!) "
		stop
	endif
	
	id = is1 - is2 + 4; ! id in 1:7

	if(id<0 .or. id>7) Write(*,*)"l1,l2,is1,is2,id=",l1,l2,is1,is2,id
	
	dij = BondLengths(is1, id);
	if(dij < 0.0d0) then
		write(*,*)"rates: dij not set correctlly! fix it! "
		write(*,*)"rates: l1,l2,is1,is2 = ",l1,l2,is1,is2
		write(*,*)"rates: is, BondLengths(is1, :) below"
		do is1=1,nsites
			write(*,'(I5,3x,7G10.4)')is1, BondLengths(is1, :)	
		enddo
		stop
	endif
	!----------------------------------------

	return
	end function dij
!******************************************************


!	finds the inter-site distances
	double precision function dijl(ih,is)
	use modmain, only: sys, nsites, ways, a0
	implicit none
	integer, intent(in):: ih,is
	integer :: l1, l2
	integer::is1,is2,id

	l1 = ways(ih)%active(is); ! main site associated with the hop
	if(ih < 9 .or. ih > 24) then ! bulk hop 
		l2 = ways(ih)%sites(is); ! the second site involved
		dijl = dsqrt(sum((sys%r(l1,:) - sys%r(l2,:))**2))
	else !contact hops
		select case(ih)
			case(11,12,15,16, 21:24) ! left
				dijl = dabs(sys%r(l1,1)); ! distance to the contact/lead
			case(9,10,13,14, 17:20) ! right
				dijl = dabs(a0*(nsites/3 + 1) - sys%r(l1,1));
		end select
	endif
	! sys%r(-2:nsites+3), shifted index
	!dijl = dsqrt(sum((sys%r(l1,:) - sys%r(l2,:))**2))

	return
	end function dijl
!******************************************************
	end module rates

