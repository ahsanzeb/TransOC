
	module rates
	use modmain, only: qt,eig,itypes,mapb,mapt,maph,mapc
	implicit none

	public:: CalRates
	private:: totr, ratehcs, PenaltyArray !, Penalty

	contains

!----------------------------------------------------
! calculates total rates for all processes
!----------------------------------------------------
	subroutine CalRates()
	use modmain, only: crosshops,nolosses,nocontacts
	implicit none
	! local
	integer :: nc
	
	! nc: only for ih=1-8
	! no crosshops? nc=1-2, otherwise nc=4
	nc=2;
	if (crosshops) nc=4; 

	! calculate rate(:)%r and rate(:)%rcs(:,:) for all hops
	call totr('dphihops',nc)
	call totr('annihilation',nc)
	call totr('creation',nc)
	if (.not. nolosses) call totr('losses',1)
	if (.not. onlybulk) call totr('contacts',1)

	return
	end subroutine CalRates

!******************************************************
	subroutine totr(process,nc)
	! calculates c,s resolved rate and total rate
	!	rate(ih)%rcs(ic,is) and rate(ih)%r
	use modmain, only: rate, ways, PermSym
	implicit none
	character(len=*), intent(in):: process
	integer, intent(in) :: nc
	! local
	integer:: ih, is, ic, ih1, ih2

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
			do is=1,ways(ih)%ns,1
				do ic=1,nc
					call ratehcs(ih,ic,is)
					if(PermSym) then ! total rate = (rates for one site) * ns
						rate(ih)%rcs(ic,is)=rate(ih)%rcs(ic,is)*ways(ih)%ns
					endif
				end do
				if(PermSym) exit; ! only a single site/case for each hop type
			end do
		rate(ih)%r = sum(rate(ih)%rcs(ic,is)); ! total rate for ih hop
		end do
	case('annihilation')
	!-------------------------------------------------
	! ih:5-6 D,Phi annihilation R,L ==> D,Phi | Phi,D
	!-------------------------------------------------
		do ih=5,6
			if (ways(ih)%ns>0) then
				do ic=1,nc
					call ratehcs(ih,ic,1)
					! total rate = (rates for one site) * ns
					rate(ih)%rcs(ic,1)=rate(ih)%rcs(ic,1)*ways(ih)%ns
				end do
			endif
			rate(ih)%r = sum(rate(ih)%rcs(ic,is)); ! total rate for ih hop
		end do
	case('losses')
	!-------------------------------------------------
	! ih:25,26 kappa, gamma losses
	!-------------------------------------------------
		! caivty losses
		ih=25; ic=1; is=1
		call ratehcs(ih,ic,is)
		rate(ih)%r = rate(ih)%rcs(ic,is); ! total rate for ih hop

		! exciton losses
		ih=25; ic=1;
		do is=1,ways(ih)%ns,1
			call ratehcs(ih,ic,is)
			if(PermSym) then ! total rate = (rates for one site) * ns
				rate(ih)%rcs(ic,is)=rate(ih)%rcs(ic,is)*ways(ih)%ns
			endif
			if(PermSym) exit; ! only a single site/case for each hop type
		end do
		rate(ih)%r = sum(rate(ih)%rcs(ic,is)); ! total rate for ih hop

	case('contacts')
	!-------------------------------------------------
	! ih:9-24; R/L contacts injection/extraction of e/h
	!-------------------------------------------------
		ic=1; is=1
		do ih=9,24
			call ratehcs(ih,ic,1)
			rate(ih)%r = rate(ih)%rcs(ic,is); ! total rate for ih hop
		end do

	end select
	!-------------------------------------------------
	return
	end subroutine totr
!******************************************************
	subroutine ratehcs(ih,ic,is)
	! rates for ih hop, ic channel, is site
	! sets global variable rate(ih)%rcs(ic,is)
	use modmain, only: nx,qt,eig,itypes,mapt,maph,mapc,Einit,rate
	implicit none
	integer, intent(in) :: ih,ic,is
	! local
	double precision, allocatable :: de(:)
	integer :: nsec


	! if no available hops, set rate = 0
	if (ways(ih)%ns == 0) then
		rate(ih)%rcs(ic,is) = 0.d0
		return
	endif

	! conditions on nx for ih=7,8
	if ((ih==7 .or. ih==8) then
		if ( (ic < 3 .and. nx .lt. 1) .or. 
     .   (ic == 3 .and. nx .lt. 2) ) then
			rate(ih)%rcs(ic,is) = 0.d0
			return 
		endif
	endif
	
	! locations of data
	itl = mapt%map(itypes(ih,ic)) ! location of final hilber space
	ia = maph(ih,ic); ! location of amplitudes
	icl = mapc(ih,ic); ! localtion of ic; icl=ic except ih=8, ic 1,2 swapped

	! no of degenerate sectors
	nsec = eig(itl)%nsec 
	allocate(de(nsec))
	
	! change in energy for all transitions
	de = eig(itl)%esec(:) + dqc(ih,ic); ! changes due to efield, barries, etc.
	de = de - Einit; ! total change in energy

	! rates; total for this hop/hchannel/site
	rate(ih)%rcs(ic,is) =
     .   sum(PenaltyArray(de,nsec) * qt(ia)%cs(icl,is)%amp2(i))

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
		double precision, dimension(ne):: PenaltyArray
		! local
		integer:: i
		penalty = 1.0d0;
		do i=1,ne
			if(de(i) > 0.0d0)PenaltyArray(i)=dexp(-de(i)*beta);
		end do
	return
	end function PenaltyArray
!******************************************************
	double precision function penalty(de)
	use modmain, only: beta
	implicit none
		penalty = 1.0d0;
		if(de > 0.0d0 ) penalty = dexp(-de*beta);
	return
	end function penalty
!******************************************************
	end module rates

