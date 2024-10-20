
	module modq
	use modmain, only: sys, nsites, Kq, ways, periodic, leads,
     .    nproc, a0, Ecoul, Etotq, signEr, Er, vrh, coulomb,
     .    impocc, impurity,vqout
	implicit none
	double precision, allocatable, dimension(:) :: Vq ! Coulomb potential

	!public :: SetEcoul
	!private::
	contains
!******************************************************
!	calculates potential Vq and total coulomb energy Eqtot
! todo?: do we need to have Etotq at all? if not, remove it.
	subroutine SetEcoul(node,iter) !Er)
	implicit none
	integer, intent(in):: node, iter
	!double precision, intent(in) ::	Er
	integer :: ih,is, l1,l2,i,j
	double precision :: rij, rim(3), x(2)
	character*20:: fname
	integer :: iu

	if(.not. coulomb) then ! applied Electric field changes & return
		do ih=1,nproc
			if(ih==25 .or. ih==26) cycle; ! kappa, gamma excluded. or set=0?
			!	ih=25,26 ways not defined, so better to exclude them compared to defining and setting to 0
			if(allocated(Ecoul(ih)%dEq)) deallocate(Ecoul(ih)%dEq)
			allocate(Ecoul(ih)%dEq(ways(ih)%ns)) ! allocate works even if 0?
			if(allocated(ways(ih)%rij)) deallocate(ways(ih)%rij)
			allocate(ways(ih)%rij(ways(ih)%ns))
			do is=1,ways(ih)%ns
				if((vrh .or. ih>34)) then ! positional disorder
					! find the indeces of sites l1,l2 for 
					!	electron jumping from l1 to l2;
					call ElectronJumpSites(ih,is,l1,l2)
					!	dE for all hops & their sites due to applied Electric field
					!x = ErChange(l1,l2);
					ways(ih)%rij(is) = sys%dij(l1,l2);!x(1);
					!write(*,*) "coulom: l1,l2, rij. dE = ",l1,l2, x
					Ecoul(ih)%dEq(is) = - Er * sys%xij(l1,l2)!x(2);
				else 	! regular lattice/triangles, no positional disorder
							!or impurity level hop that has no sense of L/R/U/D that
							! could allow using signEr(ih)*Er style
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
		sys%q(i) = -(sys%occ(i)-1)
	enddo
	if(impurity) then
		i = 4;
		sys%q(i) = -(sys%occ(i)-impocc); ! impocc is neutral state occup of imp level
	endif

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
		!Vq(i) = Vq(i) * Kq;
		!Etotq = 	Etotq + Vq(i)*sys%q(i);
	enddo
	! double counting correction
	!Etotq = 0.5*Etotq;
	!write(*,*) "real sites: Vq = ",Vq
	!--------------------------------------------
	if(leads) then ! if not needed, routine called only if device geometry 
									! but might need to test onlybulk geometry

	! Note: potential of real charge at image site is used to calculate dEq
	! as an equivalent of the change in energy of all real charges due to change in an image charge
	! It saves calc potential of image charge at locations of all real charges and then their energies.
									
		! image sites at left
		do i=-2,0; 
			rim = rimage(i+3);
			do j=1,nsites; ! only real, second mirror images ignored?
				rij = dsqrt(sum((rim - sys%r(j,:))**2))
				!write(*,*)'coulomb; im: rij',rij
				Vq(i) = Vq(i) + sys%q(j)/rij; ! potential of real charge at image site;
				Vq(j) = Vq(j) + sys%q(i+3)/rij; ! potential of image charge at real site
			enddo
			!Vq(i) = Vq(i) * Kq;
			!Etotq = 	Etotq - Vq(i)*sys%q(i+3)		
		enddo
		!--------------------------------------------
		! image sites at right
		do i=nsites+1,nsites+3
			rim = rimage(i-3);
			do j=1,nsites; ! only real, second mirror images ignored?
				rij = dsqrt(sum((rim - sys%r(j,:))**2))
				Vq(i) = Vq(i) + sys%q(j)/rij; ! potential of real charge at image site;
				Vq(j) = Vq(j) + sys%q(i-3)/rij; ! potential of image charge at real site
			enddo
			!Vq(i) = Vq(i) * Kq;
			!Etotq = 	Etotq - Vq(i)*sys%q(i-3)
		enddo
	endif

	! scale by Kq
	Vq = Kq * Vq;
	
	!write(*,*)"-------------------"
	!write(*,'(16f8.3)') 0.,0.,0.,sys%q,0.,0.,0.
	!write(*,'(16f8.3)') Vq


	!---------------------------------------
	if(vqout)then
	write(fname,'("Vq_",I4.4,".out")') node
	iu = 150+node;
	open(iu,file=trim(fname), action='write',position='append')
	if(iter==1)then
		write(iu,'("       ")')
		write(iu,'("       ")')
	endif
	do i=1,nsites
		write(iu,'(i10,3x,2G18.10)') (i-1)/3, Vq(i), sys%q(i)
	enddo
	write(iu,'("       ")')
	close(iu)
	endif
	!---------------------------------------

	!write(*,*) "Vq = ",Vq
	!--------------------------------------------
	! finally, set dEq
	!--------------------------------------------
	!	calculate change in coulomb energy for all hops & their sites
	do ih=1,nproc
		if(ih==25 .or. ih==26) cycle; ! kappa, gamma excluded. or set=0?
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
			!write(*,*)"ih,l1,l2,dEq = ",ih,l1,l2,Ecoul(ih)%dEq(is)
			if(vrh .or. ih > 34)then
				! Applied Electric field, Er term
				!x = ErChange(l1,l2);
				ways(ih)%rij(is) = sys%dij(l1,l2)!x(1);
				Ecoul(ih)%dEq(is) = Ecoul(ih)%dEq(is) - Er*sys%xij(l1,l2) !x(2);
			else! regular lattice/triangles, no positional disorder
				ways(ih)%rij(is) = a0;
				Ecoul(ih)%dEq(is) = Ecoul(ih)%dEq(is)-signEr(ih)*Er
			endif
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
	!write(*,*)"ErChange: nsites, l1,l2, z1,z2 = ",nsites, l1,l2, z1,z2 
	!-----------------------------------------
	! periodic
	if(periodic) then ! zone 1 not possible, only zone 2,3
		r1 = sys%r(l1,:); r2 = sys%r(l2,:);
		if (z==4) then ! 4=2+2, same or opp ends? 
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
		else ! z=5,6 ! deep bulk, no periodic next cell involved
			drx = sys%r(l2,1) - sys%r(l1,1);
		endif
		ErChange(1) = dsqrt(sum((r1-r2)**2)); ! nns distance, not just x-comp
		ErChange(2) = - Er * drx/a0; ! Er is scaled with a0
		!write(*,*)"ErChange: rij = ",ErChange(1)
		return
	endif
	!-------------------------------------
	!leads: device
	if(z > 3)then ! a bulk hop, no contact involved
		! x-comp of displacement of hopping electron
		drx = sys%r(l2,1) - sys%r(l1,1);
		ErChange(1) = dsqrt(sum((sys%r(l1,:)-sys%r(l2,:))**2));
		!write(*,*)"coulomb: l1,l2,rij = ",l1,l2,ErChange(1)
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
				drx = -sys%r(l1,1); ! -x direction;
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
	!write(*,*)" drx = ",drx
	return
	end function ErChange
!******************************************************

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

	! See Processes Table for process indexes, and UpdateWays() for conventions
	! D,Phi|Phi,D annihilation/creation:
	! creat: ih=7,8,33,34: D always in active (la=ways(ih)%active)
	! annih: ih=5,6,31,32: D always in active (la=ways(ih)%active)
	select case(ih) 
	case(1,2,7,8,27,28,33,34) ! electron jumps from active to other
		! lo to la electron jump
		l1 = ways(ih)%sites(is); l2 = ways(ih)%active(is); 
	case(3:6,29:32) ! electron jumps from other to active
		! la to lo electron jump
		l1 = ways(ih)%active(is); l2 = ways(ih)%sites(is); 
	case(9:10,18,20)	! e extraction/hole injection at right contact
		! la to contact electron jump
		l1 = ways(ih)%active(is); l2 = nsites+1; 
	case(13:14,17,19)	! e injection/h extraction at right contact
		! contact to la electron jump
		l2 = ways(ih)%active(is); l1 = nsites+1; 
	case(11:12,22,24)! la to left contact electron jump
		l1 = ways(ih)%active(is); l2 = 0; 
	case(15:16,21,23)! left contact to la electron jump
		l2 = ways(ih)%active(is); l1 = 0; 
	case(35:36,40,42) ! e jumps from molecule to impurity level: like 11,12,22,24
		! impurity is assumed to be on site 4
		l2 = 4; l1 = ways(ih)%active(is);
	case(37:39,41) ! e jumps to molecule from impurity level: like 15,16,21,23
		l1 = 4; l2 = ways(ih)%active(is);
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
		elseif(l1 >= nsites-2)then
			! 2 x mirror position - site's pos
			rimage(1) = 2*a0*(nsites/3 + 1) - sys%r(l1,1);
		else
			write(*,'("Error(rimage): site ",I3.3" not at interface!")')l1
			stop
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

	!write(*,*)'EqChange1: ',l1,l2,r12
	
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

	!write(*,*)'EqChange2: l1,l2 = ',l1,l2
	
	r12 = dsqrt(sum((sys%r(l1,:) - sys%r(l2,:))**2))
	! identify the interfact site to find its image charge
	! l0 = interface site index
	if( l1<= 3 .or. l1 >= nsites-2) then ! l1 is at interface
		q3 = -q1; q3p = -q1p;
		q30 = -1.0d0; !-q10; !q1 - q1p;
		rim = rimage(l1);
		l3 = indimage(l1);
	else !if( l2<= 3 .or. l2 >= nsites-2)then! l2 is at interfact
		q3 = -q2; q3p = -q2p;
		q30 = 1.0d0; !-q20; !q2 - q2p;
		rim = rimage(l2);
		l3 = indimage(l2);
	!else
	!	write(*,'("Error(EqChange2): None of l1,l2 at interface")')
	!	stop
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
     .         (v1p*q1p + v2p*q2p + v3p*q3p)

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

	!write(*,'(a,2i3.3,6f10.4)')'EqChange3:',l1,l2,r13,r23,r14,r24

	
	v1p = Vq(l1) + Kq*(-1/r12 + q30/r13 + q40/r14 ); !(-q2+q2p)=-1
	v2p = Vq(l2) + Kq*( 1/r12 + q30/r23 + q40/r24 ); !(-q1+q1p)=+1

	l3 = indimage(l1); l4 = indimage(l2);
	v3p = Vq(l3) + Kq*( 1/r13 - 1/r23);! (-q1+q1p)=+1; (-q2+q2p)=-1
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
		l1 = s2; ! real site index
		l2 = indimage(l1);	! index of imag of l1, need for Vq
		q1 = sys%q(l1); q1p=q1-1; ! e jump to s2(=l1)
	else  ! s1 real site
		l1 = s1;
		l2 = indimage(l1);	
		q1 = sys%q(l1); q1p=q1+1; ! e jump from s1(=l1)
	endif
	! image charges
	q2 = -q1; q2p = -q1p;
	! pos of image of l1
	rim = rimage(l1);
	! distance between im and real
	r12 = dsqrt(sum((sys%r(l1,:) - rim)**2));
	
	v1p = Vq(l1) + Kq*(-q2+q2p)/r12;
	v2p = Vq(l2) + Kq*(-q1+q1p)/r12; !treat image as real, see EqChange2/3, same energy changes
	
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
	end module modq
