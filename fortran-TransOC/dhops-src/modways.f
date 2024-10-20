
	module modways
	use modmain, only: sys, ways, periodic, nsites,leads, impurity,
     .               nolosses, Asites
	implicit none
	integer, dimension(16) :: signdEQ ! can be moved to modmain, but no other module needs this so its better to put it here.

	public :: UpdateWays, UpdateOcc,UpdatedEQs
	private:: WhichBulkHop

	contains
!**********************************************
	subroutine UpdateWays()
	implicit none
	! local
	integer:: ntot, p,la,lo,is,is1, ndim,np
	integer, dimension(42) :: nps
	integer, dimension(24) ::
     .    plist = (/1,2,3,4,5,6,7,8,
     .       27,28,29,30,31,32,33,34,
     .       35,36,37,38,39,40,41,42 /);
	integer, allocatable, dimension(:,:) :: act,oth
	integer :: l,r, i,p2

	ways(:)%ns = 0;




	! set gloable variable ways	
	do i=1,1; !np;
		p = 1; !plist(i);
		ways(p)%ns = 1; !nps(p)
		if(allocated(ways(p)%active)) deallocate(ways(p)%active)
		if(allocated(ways(p)%sites)) deallocate(ways(p)%sites)
			allocate(ways(p)%active(1))
			allocate(ways(p)%sites(1))
			ways(p)%active = 2
			ways(p)%sites = 1
	end do


	return
	! maz, May 2019, dhops 
	




















	if(impurity) then ! include hop processes involving impurity: ih=35-42
		ndim = 42; np = 24; 
	else
		ndim = 34; np = 16;
	endif
	
	! initialise
	nps = 0	
	! copy current sys
	ntot = sys%nsites
	allocate(act(ndim,2*ntot)) ! every site linked with 4, so number of hops could be > nsites especially for two A,A => D,Phi etc. So, use a bigger number, say, 2*nsites to avoid such a situation.
	allocate(oth(ndim,2*ntot))
	act = -2
	oth = -2


	! left,right <===> up,down : left===up, right===down

	!-----------------------------------------
	! bulk processes
	! left-right
	do is=1,ntot ! loop over left/up site	
		! choose the right site
		if (is < ntot-2) then
			is1 = is+3;
		elseif (periodic) then
			is1 = is - ntot + 3;
		else
			cycle
		endif

		! a trap/dopant site (assumed to be site no 4) is involved?
		if(impurity .and. (is==4 .or. is1==4) ) then
			call WhichImpurityHop(is,is1,p,la,lo)
		else	 ! both sites regular molecules
			! find process for left-right sites
			call WhichBulkHop(is,is1,p,la,lo)
		endif
		
		!write(*,*) "=======>>> is, p, la,lo =",is, p, la,lo 
		if (p .gt. 0 ) then
				nps(p) = nps(p) + 1;
				act(p,nps(p)) = la
				oth(p,nps(p)) = lo
				if (p .eq. 7) then ! count both 7,8
					! la,lo = left,right sites for ih=7,8
					nps(8) = nps(8) + 1; 	!7: D,Phi; 8: Phi,D
					act(8,nps(8)) = lo 		! put d in active?
					oth(8,nps(8)) = la 
				elseif(p>34) then
					if(p == 39 .or. p==40) then ! count both  39,41; 40,42
						p2 = p+2;
					else !if(p==35) then
						p2 = p+1; ! 35,36;  37,38; 
					endif	
					nps(p2) = nps(p2) + 1;
					act(p2,nps(p2)) = la;
				endif		
		endif
		
	enddo
	!-----------------------------------------
	! bulk : up-down ! 
	! 27 Nov 2017: 	note for future changes in structure: 
	! L/R/U/D bulk hops can be treated on equal footing now since we use rij's
	do is=1,ntot
		! choose the down site
		if(mod(is,3)==0) then ! site=1,2,3 in first triangle; 4,5,6 in second, ...
			is1 = is-2; ! thrid in the triangles connect back to the first
		else
			is1 = is+1;
		endif

		! a trap/dopant site (assumed to be site no 4) is involved?
		if(impurity .and. (is==4 .or. is1==4) ) then
			call WhichImpurityHop(is,is1,p,la,lo)
		else	 ! both sites regular molecules
			! find process for up-down sites
			call WhichBulkHopUD(is,is1,p,la,lo)
		endif
		
		if (p .gt. 0 ) then
			nps(p) = nps(p) + 1;
			act(p,nps(p)) = la
			oth(p,nps(p)) = lo
			if (p .eq. 33) then ! count both 33,34
				! la,lo = up,down sites for ih=33,34
				nps(34) = nps(34) + 1; 	!33: D,Phi; 34: Phi,D
				act(34,nps(34)) = lo 		! put d in active?
				oth(34,nps(34)) = la 		
			elseif(p>34) then
				! no difference between L/R/U/D hops as now we use rij etc...
				if(p == 39 .or. p==40) then ! count both  39,41; 40,42
					p2 = p+2;
				else !if(p==35) then
					p2 = p+1; ! 35,36;  37,38; 
				endif	
				nps(p2) = nps(p2) + 1;
				act(p2,nps(p2)) = la;
			endif		

		endif
			
	end do
	!-----------------------------------------

	!write(*,*) 'nps(:) = ',nps(:)
	!write(*,*) 'sys%occ = ',sys%occ

	! set gloable variable ways	
	do i=1,np;
		p = plist(i);
		ways(p)%ns = nps(p)
		if(allocated(ways(p)%active)) deallocate(ways(p)%active)
		if(allocated(ways(p)%sites)) deallocate(ways(p)%sites)
		if (nps(p) .gt. 0) then
			!write(*,*) "nps(p) > 0 ?????   p, nps(p) = ",p, nps(p)		
			allocate(ways(p)%active(nps(p)))
			allocate(ways(p)%sites(nps(p)))
			ways(p)%active = act(p,1:nps(p))
			ways(p)%sites = oth(p,1:nps(p))
		endif
	end do


	!-----------------------------------
	! contact processes
	!-----------------------------------
	act = 0;
	
	if (leads) then
	ways(9:24)%ns = 0;

	do i=1,3
	l = sys%occ(i); ! site nns of left contact
	r = sys%occ(nsites-3+i); ! site nns of right contact

	select case(l)
		case(0)
			ways(15:16)%ns = ways(15:16)%ns + 1
			act(15:16,ways(15)%ns) = i; ! left contact
		case(1)
			ways(21:24)%ns = ways(21:24)%ns + 1
			act(21:24,ways(21)%ns) = i; ! left contact
		case(2)
			ways(11:12)%ns = ways(11:12)%ns + 1
			act(11:12,ways(11)%ns) = i; ! left contact
	end select

	select case(r)
		case(0)
			ways(13:14)%ns = ways(13:14)%ns + 1
			act(13:14,ways(13)%ns) = nsites-3+i;	! right contact
		case(1)
			ways(17:20)%ns = ways(17:20)%ns + 1
			act(17:20,ways(17)%ns) = nsites-3+i;	! right contact
		case(2)
			ways(9:10)%ns = ways(9:10)%ns + 1
			act(9:10,ways(9)%ns) = nsites-3+i;	! right contact
	end select

	enddo
	!-----------------------------------
	do p=9,24
		if(allocated(ways(p)%active)) deallocate(ways(p)%active)
		if (ways(p)%ns .gt. 0) then
			allocate(ways(p)%active(ways(p)%ns))
			ways(p)%active = act(p,1:ways(p)%ns)
		endif
	enddo
	!-----------------------------------
	endif !  leads

	deallocate(act,oth)

	if(.not. nolosses) then
		p=26;
		if(allocated(ways(p)%active)) deallocate(ways(p)%active)
		ways(p)%ns = sys%n1;
		if (ways(p)%ns .gt. 0) then
			allocate(ways(p)%active(ways(p)%ns))
			ways(p)%active = Asites
		endif
	endif

	return
	end subroutine UpdateWays
!**********************************************
	subroutine WhichImpurityHop(is,is1,p,la,lo)
	implicit none
	integer, intent(in) :: is,is1
	integer, intent(out) :: p,la ,lo
	! local
	integer:: l,r ! occupation on left, right sites
	integer:: iimp, irm ! ind_impurity, ind_regular-molecule
	! ih=35-42 similar to Left contact hop (11,12,15,16,21,22,23,24)
	! but ways conditional on impurity occupation

	! which site is impurity (no. 4)?
	if(is==4) then
		iimp = is; irm = is1;
	else
		iimp = is1; irm = is;
	endif

	l = sys%occ(iimp);
	r = sys%occ(irm);

	if (l==0) then
		! empty/0+ Phi
		if (r==0)then
			p=-1; la=-1;lo=-1; ! no way
			return
		elseif(r==1) then ! Phi creaation :: 22,24 like
			p=40; !40,42:  p=22,24, like Phi creation
		elseif(r==2) then
			p=35; ! 35,36: p=11,12 like D annihilation
		else
			write(*,*)"Error(modways): occ(irm) > 2 "
			stop
		endif
	elseif(l==1)then
		if (r==0)then
			p=37; !37,38: p=15,16 like Phi annihilation
		elseif(r==1) then 
			p=39; !39,41 p=21,23, like D creation
		elseif(r==2) then
			p=-1; la=-1;lo=-1; ! no way
			return
		else
			write(*,*)"Error(modways): occ(irm) > 2 "
			stop
		endif
	else
		write(*,*)"Error(modways): occ(iimp) > 1 "
		stop
	endif

	la = irm;
	lo = iimp; ! =4 at the moment. future?: for multiple impurities, generalise this. 
	return
	end subroutine WhichImpurityHop
!**********************************************






!-----------------------	-----------------------	
! 	WhichBulkHop(l,r):
!	 find which hop is possible
!	 for two nearest neighbour sites
!	 with occupancies l,r for left and right sites
!-----------------------	-----------------------	
! list of hops:
! 1, D hops Right
! 2, D hops Left
! 3, Phi hops right
! 4, Phi hops left
! 5, D, Phi annihilated
! 6, Phi,D annihilated
! 7, D, Phi created
! 8, Phi,D created
	subroutine WhichBulkHop(is,is1,p,la,lo)
	implicit none
	integer, intent(in) :: is,is1
	integer, intent(out) :: p,la,lo 
	! process ind, active site, other site (or phi if D,Phi annihilation)
	! local
	integer:: l,r ! occupation on left, right sites

	l = sys%occ(is);
	r = sys%occ(is1)
	! both D or both Phi
	if ((l==0 .and. r==0).or.(l==2 .and. r==2))then
		p=-1; la=-1;lo=-1; ! no way
	return
	endif

	! both active
	if ((l==1 .and. r==1))then
		p = 7; la=is;lo=is1 ! 7,8 both possible: phi,d | d,phi created
	return
	endif

	! left active
	if (l==1) then
		la=is;lo=is1;
		if (r==0) then ! phi
			p = 4 ! phi hops left
		else ! D
			p = 2 ! D hops left
		endif
	return
	endif

	! right active
	if (r==1) then
		la=is1;lo=is;
		if (l==0) then ! phi
			p = 3 ! phi hops right
		else ! D
			p = 1 ! D hops right
		endif
	return
	endif

	! left D
	if (l==2) then
		! right must be a phi, since all other options crossed above.
		la=is;lo=is1; ! la=D
		p = 5 ! d,phi annihilated
	elseif (l==0) then 	! left Phi
		! right must be a D, since all other options crossed above.
		la=is1;lo=is; ! la=D
		p = 6 ! phi,d annihilated
		!write(*,*) "=============== la,lo = ",la,lo		

	else
		p = -1; la=-1;lo=-1;
		write(*,*) "modways: something wrong .... "
		stop
	endif


	return
	end subroutine WhichBulkHop
!**********************************************


	subroutine WhichBulkHopUD(is,is1,p,la,lo)
	implicit none
	integer, intent(in) :: is,is1
	integer, intent(out) :: p,la,lo 
	! process ind, active site, other site (or phi if D,Phi annihilation)
	! local
	integer:: l,r ! occupation on up, down sites

	l = sys%occ(is);
	r = sys%occ(is1);
	! both D or both Phi
	if ((l==0 .and. r==0).or.(l==2 .and. r==2))then
		p=-1; la=-1;lo=-1; ! no way
	return
	endif

	! both active
	if ((l==1 .and. r==1))then
		p = 33; la=is;lo=is1 ! 33,34 both possible: phi,d | d,phi created
	return
	endif

	! left active
	if (l==1) then
		la=is;lo=is1;
		if (r==0) then ! phi
			p = 30 ! phi hops left
		else ! D
			p = 28 ! D hops left
		endif
	return
	endif

	! right active
	if (r==1) then
		la=is1;lo=is;
		if (l==0) then ! phi
			p = 29 ! phi hops right
		else ! D
			p = 27 ! D hops right
		endif
	return
	endif

	! left D
	if (l==2) then
		! right must be a phi, since all other options crossed above.
		la=is;lo=is1; ! la=D
		p = 31 ! d,phi annihilated
	elseif (l==0) then 	! left Phi
		! right must be a D, since all other options crossed above.
		la=is1;lo=is; ! la=D
		p = 32 ! phi,d annihilated
		!write(*,*) "=============== la,lo = ",la,lo		

	else
		p = -1; la=-1;lo=-1;
		write(*,*) "modways: something wrong .... "
		stop
	endif


	return
	end subroutine WhichBulkHopUD

!**********************************************

	subroutine UpdateOcc(ih,is)
	! only for bulk hoppings at the moment
	use modmain, only: sys,Asites,ways,nsites
	use lists, only: Drop4, Join, Substitute,MemberQ4
	implicit none
	integer, intent(in):: ih,is
	! local
	integer :: la,lo,n0,n1,n2,two=2
	integer, dimension(sys%n1+2) :: Asitesx,Asitesy
	integer:: i

	!write(*,*) "in: Asites = ",Asites
	!write(*,*) "in: occ = ",sys%occ
	
	if(ih <= 8 .or. (ih >= 27 .and. ih <= 34)) then
		la = ways(ih)%active(is)
		lo = ways(ih)%sites(is)
	elseif(ih ==25 .or. ih == 26) then
		! kappa/gamma, no change in occ
		return
	elseif(ih >= 9 .and. ih <= 24) then ! contacts
		la = ways(ih)%active(is) ! site involved in the process
	elseif(ih>=35 .and. ih <= 42) then ! impurity
		la = ways(ih)%active(is) ! regular site involved in the process
		lo = 4; ! impurity at site 4
	endif

	Asitesx = 0
	
	n0 = sys%n0; ! no of phi site
	n1 = sys%n1; ! no of active site
	n2 = sys%n2; ! no of d site

	!write(*,*)"la,lo=",la,lo
	!write(*,*)"n0,n1,n2=",n0,n1,n2
	!write(*,*) "sys%occ 1=",sys%occ
		
	select case(ih)	
		case(1,2,27,28)
			if(sys%occ(la)==1 .and. sys%occ(lo)==2)then
				sys%occ(la) = 2
				sys%occ(lo) = 1
			else
				write(*,*)"Error(UpdateOcc): case 1,2,27,28"
				stop
			endif
			call Substitute(Asites,n1,la,lo)
		case(3,4,29,30)
			if(sys%occ(la)==1 .and. sys%occ(lo)==0)then
				sys%occ(la) = 0
				sys%occ(lo) = 1
			else
				write(*,*)"Error(UpdateOcc): case 3,4,29,30"
				stop
			endif
			call Substitute(Asites,n1,la,lo)
		case(5,6,31,32)
			if(sys%occ(la)==2 .and. sys%occ(lo)==0)then
				sys%occ(la) = 1
				sys%occ(lo) = 1
			else
				write(*,*)"Error(UpdateOcc): case 5,6,31,32"
				stop
			endif
						
			Asitesx(1:n1) = Asites;
			! append first D then Phi; la=D for ih=5,6
			call Join(Asitesx(1:n1),n1,(/la,lo/),two,Asitesy(1:n1+2))
			deallocate(Asites); allocate(Asites(n1+2))
			Asites = Asitesy;

			sys%n0 = n0-1
			sys%n1 = n1+2
			sys%n2 = n2-1
		case(7,8,33,34)
			sys%occ(la) = 2 ! la=D,lo=phi for ih=7,8
			sys%occ(lo) = 0
			Asitesx(1:n1) = Asites;
			call Drop4(Asitesx(1:n1),n1,la,Asitesy(1:n1-1))
			call Drop4(Asitesy(1:n1-1),n1,lo,Asitesx(1:n1-2))
			deallocate(Asites); allocate(Asites(n1-2))
			Asites = Asitesx(1:n1-2);		
			sys%n0 = n0+1
			sys%n1 = n1-2
			sys%n2 = n2+1

		! contacts -----------------------------------
		case(9,10,13,14)! right contact annihilation, site=nsites
			sys%occ(la) = 1
			Asitesx(1:n1) = Asites;
			Asitesx(n1+1) = la;
			deallocate(Asites); allocate(Asites(n1+1))
			Asites = Asitesx(1:n1+1);
			sys%n1 = n1+1
			if(ih > 10) then
				sys%n0 = n0-1
			else
				sys%n2 = n2-1
			endif
		case(11,12,15,16) ! left contact annihilation, site=1
			sys%occ(la) = 1
			Asitesx(1:n1) = Asites;
			Asitesx(n1+1) = la;
			deallocate(Asites); allocate(Asites(n1+1))
			Asites = Asitesx(1:n1+1);
			sys%n1 = n1+1
			if(ih > 12) then
				sys%n0 = n0-1
			else
				sys%n2 = n2-1
			endif

		case(17:20) ! right contact creation, site=nsite
			call Drop4(Asites,n1,la,Asitesy(1:n1-1))
			deallocate(Asites); allocate(Asites(n1-1))
			Asites = Asitesy(1:n1-1);	
			sys%n1 = n1-1
			if(mod(ih,2)==0) then ! Phi created, ih=18,20
				sys%occ(la) = 0;
				sys%n0 = n0+1
			else ! D created, ih=17,19
				sys%occ(la) = 2;
				sys%n2 = n2+1
			endif

		case(21:24) ! left contact creation, site=1
			call Drop4(Asites,n1,la,Asitesy(1:n1-1))
			deallocate(Asites); allocate(Asites(n1-1))
			Asites = Asitesy(1:n1-1);	
			sys%n1 = n1-1
			if(mod(ih,2)==0) then ! Phi created, ih=22,24
				sys%occ(la) = 0;
				sys%n0 = n0+1
			else ! D created, ih=21,23
				sys%occ(la) = 2;
				sys%n2 = n2+1
			endif

		! impurity -----------------------------------
		case(35:38)! at impurity, annihilation of D or Phi
			sys%occ(la) = 1
			Asitesx(1:n1) = Asites;
			Asitesx(n1+1) = la;
			deallocate(Asites); allocate(Asites(n1+1))
			Asites = Asitesx(1:n1+1);
			sys%n1 = n1+1
			if(ih <= 36) then ! ih=35,36; D annihil
				sys%n2 = n2-1;
				sys%occ(lo) = 1; !0 -> 1; sys%occ(lo) + 1; ! impurity 1LS, max occ = 1;
			else ! ih=37,38; Phi annihil
				sys%n0 = n0-1
				sys%occ(lo) = 0; !1 -> 0; sys%occ(lo) - 1;
			endif
			
		case(39:42) ! at impurity, creation of D or Phi
			call Drop4(Asites,n1,la,Asitesy(1:n1-1))
			deallocate(Asites); allocate(Asites(n1-1))
			Asites = Asitesy(1:n1-1);	
			sys%n1 = n1-1
			if(mod(ih,2)==0) then ! Phi created, ih=40,42
				sys%occ(la) = 0;
				sys%n0 = n0+1;
				sys%occ(lo) = 1; ! 0 --> 1; impurity 1LS, max occ = 1;
			else ! D created, ih=39,41
				sys%occ(la) = 2;
				sys%n2 = n2+1;
				sys%occ(lo) = 0; ! 1 --> 0
			endif

	end select		

	if(sys%n0+sys%n1+sys%n2+sys%nimp .ne. sys%nsites)then
		write(*,*)"modways: nsites != n0+n1+n2+nimp "
		stop
	endif


	!write(*,*) "out: Asites = ",Asites
	!write(*,*) "out: occ = ",sys%occ
	!write(*,*) "-----------------------"	

	do i=1,nsites
		if(sys%occ(i) == 1) then
			if(.not. impurity .or. i .ne. 4) then
			if(.not. MemberQ4(Asites,sys%n1,i)) then
				write(*,*)"Error(UpdateOcc): Asites not updated correctly!"
				write(*,*)"Asites=",Asites
				write(*,*)"occ=",sys%occ
				write(*,*)" ih,is = ",ih,is
			endif
			endif
		endif
	enddo




	return
	end subroutine UpdateOcc
!**********************************************
! signdEQ = 	sign of charging energy for various contact processes
! i.e., electron/hole extraction/injection
!----------------------------------------
	subroutine UpdatedEQs(Qtot)
	use modmain, only: Eq, dEQs
	! in: Qtot, Eq;
	! out: dEQs
	implicit none
	integer, intent(in)::Qtot
	integer:: i, Qtotf
	double precision :: Eqi
	! Net charge on the system
	!	Qnet = (nsites-nelec)-sum(sys%occ);
	! charging energy
	Eqi = ChargingE(Qtot);
	! change in charging energy for various contact processes
	!	i.e., electron/hole injection/extraction
	! similar charging increase energy, opposite decreases
	do i=1,16
		Qtotf = Qtot + signdEQ(i);
		dEQs(i) = ChargingE(Qtotf) - Eqi;
	enddo
	!write(*,*)'Eqi, Eq = ',Eqi, Eq
	!write(*,*)"dEQs = ",dEQs
	!write(*,*) "Qtotf = ", Qtot + signdEQ
	!write(*,*) "Qtot = ", Qtot
	!write(*,*) "Sign = ", signdEQ
	
	return
	end 	subroutine UpdatedEQs
!-----------------------------------
	double precision function ChargingE(q)
	use modmain, only: Eq
	implicit none
	integer, intent(in) :: q
	integer :: absq
	absq = abs(q);
	ChargingE = Eq * (absq - 1) * absq / 2.0d0
	return
	end function ChargingE
!-----------------------------------

	end module modways
