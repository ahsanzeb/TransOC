
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
	use modmain, only: crosshops,nolosses,leads
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
	!write(*,*) "rates  1"
	call totr('annihilation',nc)
	!write(*,*) "rates  2"
	call totr('creation',nc)
	!write(*,*) "rates  3"
	if (.not. nolosses) call totr('losses',1)
	if (leads) call totr('contacts',1)

	return
	end subroutine CalRates

!******************************************************
	subroutine totr(process,nc)
	! calculates c,s resolved rate and total rate
	!	rate(ih)%rcs(ic,is) and rate(ih)%r
	use modmain, only: rate, ways, PermSym,na,nx,ts
	implicit none
	character(len=*), intent(in):: process
	integer, intent(in) :: nc
	! local
	integer:: ih, ic, ih1, ih2,is

	select case(process)
	case('dphihops','creation')
	!-------------------------------
	! ih:1-4 DHopsR/L, PhiHopsR/L
	!-------------------------------
	! ih=7,8 (D,Phi), (Phi,D) created
	! ih=8 same amplitudes,chan 1,2 swapped: done in ratehcs()
	!--------------------------------------------------------
		ih1=1;ih2=4; ! 'dphihops'
		if(process == 'creation') then
			ih1=7;ih2=8;
		endif
		do ih=ih1,ih2
			!write(*,*)"ih, ways(ih)%ns = ",ih, ways(ih)%ns
			!write(*,'(a,i5,5x,4f10.5)')"ih, ts(ih,:) =",ih, ts(ih,:)
			rate(ih)%rcs(:,:) = 0.0d0
			do is=1,ways(ih)%ns,1
				do ic=1,nc
					!write(*,*) 'heloowwww 2' !ways(ih)%active	
					call ratehcs(ih,ic,is)
					!write(*,*)"rates:ih,ic,rcs=",ih,ic,rate(ih)%rcs(ic,is)
					if(PermSym) then ! total rate = (rates for one site) * ns
						rate(ih)%rcs(ic,is)=
     .             ts(ih,ic)*rate(ih)%rcs(ic,is)*ways(ih)%ns
					else
						rate(ih)%rcs(ic,is)=ts(ih,ic)*rate(ih)%rcs(ic,is)
					endif
				end do
				if(PermSym) exit; ! only a single site/case for each hop type
			end do
			rate(ih)%r = sum(rate(ih)%rcs); ! total rate for ih hop
			!write(*,*)"ih,rate(ih)%rcs=",ih, rate(ih)%rcs
			!write(*,*)"ih,rate(ih)%r=",ih, rate(ih)%r
		end do
	case('annihilation')
	!-------------------------------------------------
	! ih:5-6 D,Phi annihilation R,L ==> D,Phi | Phi,D
	!-------------------------------------------------

		do ih=5,6
			rate(ih)%rcs(:,:) = 0.0d0
			if (ways(ih)%ns>0) then
				do ic=1,nc
					call ratehcs(ih,ic,1)
					! total rate = (rates for one site) * ns
					rate(ih)%rcs(ic,1)=
     .           ts(ih,ic)*rate(ih)%rcs(ic,1)*ways(ih)%ns
				end do
				rate(ih)%r = sum(rate(ih)%rcs(:,:)); ! total rate for ih hop
				!write(*,*)"rates: ann:ih,ic, rcs=",ih,ic,rate(ih)%rcs
			else
				rate(ih)%r = 0.0d0;
			endif
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
				else
					rate(ih)%rcs(ic,is)=ts(ih,ic)*rate(ih)%rcs(ic,is)
				endif

				if(PermSym) exit; ! only a single site/case for each hop type

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
		do ih=9,24
			rate(ih)%rcs(:,:) = 0.0d0
			if(ways(ih)%ns==1) then
				if(ih==19 .or. ih==20 .or. ih==23 .or. ih==24) then
					if (nx > 0) call ratehcs(ih,1,1)
				else
					call ratehcs(ih,1,1)
				endif
			endif
			rate(ih)%r = rate(ih)%rcs(1,1); ! total rate for ih hop
		end do

	end select
	!-------------------------------------------------
	return
	end subroutine totr
!******************************************************
	subroutine ratehcs(ih,ic,is)
	! rates for ih hop, ic channel, is site
	! sets global variable rate(ih)%rcs(ic,is)
	use modmain, only: nx,qt,eig,itypes,beta,
     .  mapt,maph,mapc,Einit,rate,ways,dqc,dEQs
	implicit none
	integer, intent(in) :: ih,ic,is
	! local
	double precision, allocatable :: de(:)
	integer :: nsec, ia,icl,itl,i
	double precision:: x,y
	logical :: dblehop, emptyhop

	!write(*,*) 'heloowwww A'


	! if no available hops, set rate = 0
	if (ih <= 8 .and. ways(ih)%ns == 0) then
		rate(ih)%rcs(ic,is) = 0.d0
		return
	endif

	!write(*,*) 'heloowwww b'

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
	endif

	!write(*,*) 'heloowwww c'

	! conditions on nx for ih=7,8
	if (ih==7 .or. ih==8) then
		if ( (ic < 3 .and. nx .lt. 1) .or. 
     .   (ic == 3 .and. nx .lt. 2) ) then
			rate(ih)%rcs(ic,is) = 0.d0;
			!write(*,*)"rate: nx cond; ih,ic,is=",ih,ic,is
			return 
		endif
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
	de = eig(itl)%esec(1:nsec) + dqc(ih,ic); ! changes due to efield, barries, etc.
	de = de - Einit; ! total change in energy

	! charging energy contact hops
	if(ih .ge. 9 .and. ih .le. 24) then
		de = de + dEQs(ih-8)
	endif

	! rates; total for this hop/hchannel/site

	if(.not. allocated(qt(ia)%cs(icl,is)%amp2)) then
	write(*,*) "alloc amp2=>F: ih,ic,is,icl,ia",ih,ic,is,icl,ia
	stop
	endif

	rate(ih)%rcs(ic,is) =
     .   sum(PenaltyArray(de,nsec) * qt(ia)%cs(icl,is)%amp2)

	! for debugging, remove later
	if (isnan(rate(ih)%rcs(ic,is))) then
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
	use modmain, only: beta
	implicit none
		integer,intent(in) :: ne
		double precision, dimension(ne),intent(in):: de
		double precision, dimension(ne):: tmp
		double precision, dimension(ne):: PenaltyArray
		! local
		integer:: i
		double precision:: x, fac=30.0d0
		
		PenaltyArray = 1.0d0;
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
	end module rates

