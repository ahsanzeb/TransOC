
	module modq
	use modmain, only: sys, nsites, Kq, ways, periodic, leads,
     .    nproc, a0, Ecoul, Etotq, signEr, Er, vrh, coulomb
	implicit none
	double precision, allocatable, dimension(:) :: Vq ! Coulomb potential

	!public :: SetEcoul
	!private::
	contains
!******************************************************
!	calculates potential Vq and total coulomb energy Eqtot
! todo?: do we need to have Etotq at all? if not, remove it.
	subroutine SetEcoul() !Er)
	implicit none
	!double precision, intent(in) ::	Er
	integer :: ih,is, l1,l2,i,j
	double precision :: rij, rim(3), x(2)

	if(.not. coulomb) then ! applied Electric field changes & return
		do ih=1,nproc
			if(allocated(Ecoul(ih)%dEq)) deallocate(Ecoul(ih)%dEq)
			allocate(Ecoul(ih)%dEq(ways(ih)%ns)) ! allocate works even if 0?
			if(allocated(ways(ih)%rij)) deallocate(ways(ih)%rij)
			allocate(ways(ih)%rij(ways(ih)%ns))
			do is=1,ways(ih)%ns
				if(vrh) then ! positional disorder
					! find the indeces of sites l1,l2 for 
					!	electron jumping from l1 to l2;
					call ElectronJumpSites(ih,is,l1,l2)
					!	dE for all hops & their sites due to applied Electric field
					x = ErChange(l1,l2);
					ways(ih)%rij(is) = x(1);
					Ecoul(ih)%dEq(is) = x(2);
				else ! regular lattice/triangles, no positional disorder
					ways(ih)%rij(is) = a0;
					Ecoul(ih)%dEq(is) = -signEr(ih)*Er
				endif
				!write(*,*)"ih,is, Ecoul%dE = ", Ecoul(ih)%dEq(is)
			enddo
		enddo
		return
	endif

	!call qtime()
	
!	If coulomb's interaction is true, proceed as follows.

	!--------------------------------------------
	! set sys%q
	!--------------------------------------------
	do i=1,nsites
		sys%q(i) = sys%q0(i) -(sys%occ(i)-1)
	enddo
	!--------------------------------------------
	! set Vq and Etotq
	!--------------------------------------------	
	if(.not. allocated(Vq)) 	allocate(Vq(-2:nsites+3))
	Vq = 0.0d0;
	Etotq = 0.0d0;
	! ==> energy for real charges in the potential of image charges
	! is the same as image charges in the potential of real charges
	! ==> Vq at image sites will be used when calculating the 
	! change in energy due to change in image charges
	! without iterating over all real charges
	!--------------------------------------------
	! real sites
	do i=1,nsites;
		do j=1,nsites;
			if(i .ne. j) then
				rij = dsqrt(sum((sys%r(i,:) - sys%r(j,:))**2));
				Vq(i) = Vq(i) + sys%q(j)/rij
			endif
		enddo
		Vq(i) = Vq(i) * Kq;
		Etotq = 	Etotq + Vq(i)*sys%q(i);
	enddo
	! double counting correction
	Etotq = 0.5*Etotq;
	!write(*,*) "real sites: Vq = ",Vq
	!--------------------------------------------
	if(leads) then
		! image sites at left
		do i=-2,0; 
			rim = rimage(i+3);
			do j=1,nsites; ! only real, second mirror images ignored?
				rij = dsqrt(sum((rim - sys%r(j,:))**2))
				Vq(i) = Vq(i) + sys%q(j)/rij
			enddo
			Vq(i) = Vq(i) * Kq;
			Etotq = 	Etotq - Vq(i)*sys%q(i+3)		
		enddo
		!--------------------------------------------
		! image sites at right
		do i=nsites+1,nsites+3
			rim = rimage(i-3);
			do j=1,nsites; ! only real, second mirror images ignored?
				rij = dsqrt(sum((rim - sys%r(j,:))**2))
				Vq(i) = Vq(i) + sys%q(j)/rij
			enddo
			Vq(i) = Vq(i) * Kq;
			Etotq = 	Etotq - Vq(i)*sys%q(i-3)
		enddo
	endif

	!write(*,*) "Vq = ",Vq
	!--------------------------------------------
	! finally, set dEq
	!--------------------------------------------
	!	calculate change in coulomb energy for all hops & their sites
	do ih=1,nproc
		if(allocated(Ecoul(ih)%dEq)) deallocate(Ecoul(ih)%dEq)
		allocate(Ecoul(ih)%dEq(ways(ih)%ns)) ! allocate works even if 0?
		if(allocated(ways(ih)%rij)) deallocate(ways(ih)%rij)
		allocate(ways(ih)%rij(ways(ih)%ns))
		!Ecoul(ih)%dEq(:) = 0.0d0
		do is=1,ways(ih)%ns
			! find the indeces of sites l1,l2 for 
			!	electron jumping from l1 to l2;
			call ElectronJumpSites(ih,is,l1,l2)
			! calculate the coulomb energy changes
			Ecoul(ih)%dEq(is) = CoulombEnergyChange(l1,l2)
			!write(*,*)"ih,is, Ecoul%dE = ",ih,is, Ecoul(ih)%dEq(is)
			! Applied Electric field, Er term
			x = ErChange(l1,l2);
			ways(ih)%rij(is) = x(1);
			Ecoul(ih)%dEq(is) = Ecoul(ih)%dEq(is) + x(2);
			!write(*,*)"    Ecoul%dE +Er = ", Ecoul(ih)%dEq(is)
		enddo
	enddo
	!--------------------------------------------
	!write(*,*)"coulomb takes:"
	!call qtime()
	
	return
	end subroutine SetEcoul
!******************************************************
	function ErChange(l1,l2)
	implicit none
	!double precision, intent(in) :: Er
	integer, intent(in):: l1,l2
	double precision, dimension(2) :: ErChange
	double precision :: drx, r1(3),r2(3)
	integer :: z, z1, z2

	! zones l1,l2 lie in.
	z1 = zone(l1);	 z2 = zone(l2);	
	z = z1 + z2;

	if(periodic) then
		r1 = sys%r(l1,:); r2 = sys%r(l2,1);
		if (z==2) then ! same or opp ends? 
			! x-comp of displacement of hopping electron	
			if(l1<=3 .and. l2>= nsites-2)then ! opposite end layers, apply periodicity
				r2(1) = sys%r(l2,1) - a0*nsites/3;
				drx = r2(1) - sys%r(l1,1);
			elseif(l2<=3 .and. l1>= nsites-2)then ! opposite end
				r2(1) = sys%r(l2,1) + a0*nsites/3;
				drx = r2(1) - sys%r(l1,1);
			else ! same end
				drx = sys%r(l2,1) - sys%r(l1,1);
			endif
		else ! deep bulk, no periodic next cell involved
			drx = sys%r(l2,1) - sys%r(l1,1);
		endif
		ErChange(1) = dsqrt(sum((r1-r2)**2)); ! nns distance, not just x-comp
		ErChange(2) = - Er * drx/a0; ! Er is scaled with a0
		return
	endif

	if(z > 3)then ! a bulk hop, no contact involved
		! x-comp of displacement of hopping electron
		drx = sys%r(l2,1) - sys%r(l1,1);
		ErChange(1) = dsqrt(sum((sys%r(l1,:)-sys%r(l2,:))**2));
	elseif(z==3) then ! a contact hop
		if(z1==1) then ! l1 ==> contact, l2 ==> interface site
			! left/right contact?
			if(l1<1)then! left contact
				drx = sys%r(l2,1); ! +x direction;
			else ! right contact
				drx = -((nsites/3 + 1)*a0-sys%r(l2,1)); ! -x direction
			endif
		else ! z2==1; ! l2 ==> contact, l1 ==> interface site
			! left/right contact?
			if(l2<1)then! left contact
				drx = sys%r(l1,1); ! -x direction;
			else ! right contact
				drx = (nsites/3 + 1)*a0-sys%r(l1,1); ! +x direction
			endif
		endif
		ErChange(1) = dabs(drx);! dist to the contact = x-comp
	else
		write(*,*)"Error(ErChange): z1, z2 = ", z1, z2
		stop
	endif 
	
	ErChange(2) = - Er * drx/a0; ! Er is scaled with a0
	
	return
	end function ErChange
!-----------------------------------------------
! find the initial and final indeces of sites l1,l2
! for 'electron' jumps
	subroutine ElectronJumpSites(ih,is,l1,l2)
	implicit none
	integer, intent(in):: ih,is
	integer, intent(out):: l1,l2

	! Convention:
	! right contact identifier: site index > nsites 
	! left contact identifier: site index < 0 

	! See Processes Table for process indexes
	select case(ih)
	case(1,2,7,8,26,27,33,34)
		! la to lo electron jump
		l1 = ways(ih)%active(is); l2 = ways(ih)%sites(is);
	case(3:6,28:32)
		! lo to la electron jump
		l2 = ways(ih)%active(is); l1 = ways(ih)%sites(is);
	case(9:10,18,20)	! e extraction/hole injection at right contact
		! la to contact electron jump
		l1 = ways(ih)%active(is); l2 = nsites+1; 
	case(13:14,17,19)	! e injection/h extraction at right contact
		! contact to la electron jump
		l2 = ways(ih)%active(is); l1 = nsites+1; 
	case(11:12,22,24)! la to left contact electron jump
		l1 = ways(ih)%active(is); l2 = -1; 
	case(15:16,21,23)! left contact to la electron jump
		l2 = ways(ih)%active(is); l1 = -1; 
	end select

	return
	end 	subroutine ElectronJumpSites
!-----------------------------------------------
! find image charge position
	function rimage(l1)
	implicit none
	integer, intent(in)::l1
	double precision, dimension(3) :: rimage
	rimage = sys%r(l1,:);
		if(l1 <= 3 ) then
			rimage(1) = -sys%r(l1,1);! reflection along x axis
		else
			! 2 x mirror position - site's pos
			rimage(1) = 2*a0*(nsites/3 + 1) - sys%r(l1,1);
		endif
	return
	end function rimage
!====================================
! find image charge index
	double precision function indimage(l)
	implicit none
	integer, intent(in) :: l ! index of an interface site
		if(l <= 3 ) then
			indimage = l - 3;
		else
			indimage = l + 3;
		endif
	return
	end function indimage
!=====================================================
	! both sites away from the interface in deep bulk
	double precision function EqChange1(l1,l2)!,q1,q2,q1p,q2p)
	implicit none
	integer, intent(in) :: l1,l2
	!double precision, intent(in):: q1,q2,q1p,q2p
	double precision:: q1,q2,q1p,q2p
	double precision:: v1p,v2p, r12

	q1 = sys%q(l1); q2 = sys%q(l2);
	q1p = q1+1.0d0; q2p = q2-1.0d0;
	r12 = dsqrt(sum((sys%r(l1,:) - sys%r(l2,:))**2))
	
	v1p = Vq(l1) - Kq/r12 ! (-q2 + q2p) = -1
	v2p = Vq(l2) + Kq/r12 ! (-q1 + q1p) = +1
	EqChange1 = -(Vq(l1)*q1 + Vq(l2)*q2) + (v1p*q1p + v2p*q2p)
	return
	end function EqChange1
!======================================
	! one site at interface, other away from it
	! consider the image charge of the interface site 
	double precision function EqChange2(l1,l2)
	implicit none
	integer, intent(in) :: l1,l2
	!double precision, intent(in):: q1,q2,q1p,q2p
	double precision:: q1,q2,q1p,q2p, q30,q3,q3p
	double precision:: v1p,v2p,v3p, rim(3)
	double precision:: r12, r13, r23
	integer:: l3 
	q1 = sys%q(l1); q2 = sys%q(l2);
	q1p = q1+1.0d0; q2p = q2-1.0d0; ! e jumps from l1 to l2
	!q10 = 1.0d0; !q1p - q1; 
	!q20 = -1.0d0; !q2p - q2;
	
	r12 = dsqrt(sum((sys%r(l1,:) - sys%r(l2,:))**2))
	! identify the interfact site to find its image charge
	! l0 = interface site index
	if( l1<= 3 .or. l1 >= nsites) then ! l1 is at interface
		q3 = -q1; q3p = -q1p;
		q30 = -1.0d0; !-q10; !q1 - q1p;
		rim = rimage(l1);
		l3 = indimage(l1);
	else ! l2 is at interfact
		q3 = -q2; q3p = -q2p;
		q30 = 1.0d0; !-q20; !q2 - q2p;
		rim = rimage(l2);
		l3 = indimage(l2);
	endif
	r13 = dsqrt(sum((sys%r(l1,:) - rim)**2));
	r23 = dsqrt(sum((sys%r(l2,:) - rim)**2));
	!v1p = Vq(l1) + Kq*( q20/r12 + q30/r13 )
	!v2p = Vq(l2) + Kq*( q10/r12 + q30/r23 )
	v1p = Vq(l1) + Kq*( -1.0d0/r12 + q30/r13 )
	v2p = Vq(l2) + Kq*(  1.0d0/r12 + q30/r23 )
	! change in image charge field affects all charges in system
	! corresponding energy change == as if img charge was real
	v3p = Vq(l3) + Kq*( 1.0d0/r13 - 1.0d0/r23 )

	EqChange2 = -(Vq(l1)*q1 + Vq(l2)*q2 + Vq(l3)*q3) +
     .         (v1p*q1p + v2p*q2p - v3p*q1p)

	return
	end function EqChange2
!=====================================================
	! both site at interface
	! consider the image charge of both sites 
	double precision function EqChange3(l1,l2)
	implicit none
	integer, intent(in) :: l1,l2
	!double precision, intent(in):: q1,q2,q1p,q2p
	double precision:: q1,q2,q1p,q2p, q30,q40 !,q3,q3p
	double precision:: v1p,v2p, v3p, v4p, rim3(3), rim4(3)
	double precision:: r12, r13, r23, r14, r24
	integer :: l3, l4
	
	q1 = sys%q(l1); q2 = sys%q(l2);
	q1p = q1+1.0d0; q2p = q2-1.0d0; ! e jumps from l1 to l2
	r12 = dsqrt(sum((sys%r(l1,:) - sys%r(l2,:))**2));
	! q3=-q1; q4=-q2;
	! q3p=-q1p; q4p=-q2p
	
	q30 = q1 - q1p; rim3 = rimage(l1);
	q40 = q2 - q2p; rim4 = rimage(l2);

	r13 = dsqrt(sum((sys%r(l1,:) - rim3)**2));
	r23 = dsqrt(sum((sys%r(l2,:) - rim3)**2));
	r14 = dsqrt(sum((sys%r(l1,:) - rim4)**2));
	r24 = dsqrt(sum((sys%r(l2,:) - rim4)**2));
	
	v1p = Vq(l1) + Kq*(-1/r12 + q30/r13 + q40/r14 );
	v2p = Vq(l2) + Kq*( 1/r12 + q30/r23 + q40/r24 );

	l3 = indimage(l1); l4 = indimage(l2);
	v3p = Vq(l3) + Kq*( 1/r13 - 1/r23);
	v4p = Vq(l4) + Kq*( 1/r14 - 1/r24);
	
	EqChange3 = -(Vq(l1)*q1 + Vq(l2)*q2 - Vq(l3)*q1 - Vq(l4)*q2) +
     .         (v1p*q1p + v2p*q2p - v3p*q1p - v4p*q2p);
	
	return
	end function EqChange3
!=====================================================
	! contact hops: a site at interface
	! consider the image as the other site
	double precision function EqChange4(s1,s2)
	implicit none
	integer, intent(in) :: s1,s2
	double precision:: q1,q2,q1p,q2p
	double precision:: v1p,v2p, rim(3)
	double precision:: r12
	integer :: l1, l2 ! 1 real, 2 imag

	! local index: real 1, imag 2
	if(s1 < 0 .or. s1 > nsites) then ! s2 real site
		! s2 real site, e jumps from contact to the system
		! use l2 for image site of s2
		l1 = s2;
		l2 = indimage(l1);	! index of imag of l1, need for Vq
		q1 = sys%q(s2); q1p=q1-1; ! e jump to s2
	else  ! s1 real site
		l1 = s1;
		l2 = indimage(l1);	
		q1 = sys%q(s1); q1p=q1+1; ! e jump from s1
	endif
	! image charges
	q2 = -q1; q2p = -q1p;
	! pos of image of l1
	rim = rimage(l1);
	! distance between im and real
	r12 = dsqrt(sum((sys%r(l1,:) - rim)**2));
	
	v1p = Vq(l1) + Kq*(-q2+q2p)/r12;
	v2p = Vq(l2) + Kq*(-q1+q1p)/r12; !trt real, see EqChange2/3
	
	EqChange4 = -(Vq(l1)*q1 + Vq(l2)*q2) + (v1p*q1p + v2p*q2p);
	
	return
	end function EqChange4
!=====================================================
! calls appropriate function to calculate energy changes
	double precision function CoulombEnergyChange(l1,l2)
	implicit none
	integer, intent(in) :: l1,l2
	integer :: z, z1,z2

	! peridoc? warning in readinout.
	! interaction btw nns cells? do it later.... use some cutoff distance....
	! use ErChange() code with this, together sepereate routine when periodic.
	
	! zones l1,l2 lie in.
	!z1 = zone(l1);	 	z2 = zone(l2);	
	z = zone(l1) + zone(l2);
	! which function to call?
	if (z==6) then ! both bulk: 3+3
		CoulombEnergyChange = EqChange1(l1,l2)
	elseif(z==5) then ! 2+3 or 3+2
		CoulombEnergyChange = EqChange2(l1,l2)
	elseif(z==4) then ! 2+2
		CoulombEnergyChange = EqChange3(l1,l2)
	elseif(z==3) then ! 1+2 or 2+1
		CoulombEnergyChange = EqChange4(l1,l2)
	else
		write(*,*)'Error(CoulombEC): hop sites not consistent?!'
		write(*,*)' z1, z2 = ',zone(l1), zone(l2)
		stop
	endif

	!write(*,*)' z, dEq ',z,CoulombEnergyChange

	return
	end function CoulombEnergyChange
!============================================================
	integer function zone(l)
	implicit none
	integer, intent(in) :: l

	if(l < 1 .or. l > nsites) then ! imag site
		zone = 1
	elseif(l <= 3 .or. l >= nsites-2) then ! interface site
		zone = 2
	else ! 'bulk' site
		zone = 3
	endif
	return
	end function zone
!============================================================
	subroutine qtime()
	implicit none
	integer, save	:: dt(8), dti(8)=0
	character*10 :: bb(3)
	logical, save :: start = .true.;

	call date_and_time(bb(1), bb(2), bb(3), dt)

	write(6,
     . '("time: ",i2,":",i2.2,":",i2.2,"(h:m:s)")')
     .   dt(5:7)-dti(5:7)
	write(6,*) dt(5:8)-dti(5:8)
	dti = dt;
	return
	end subroutine qtime
!==========================================
	end module modq
