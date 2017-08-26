

!------------------------------------------
!	Map for D Phi Annihilation for channel 1,2
!------------------------------------------
	subroutine AnMap1(ib,n,k,mapa,la)
	!	chan 1,2
	! D&Phi sites appended to active list.
	! out: index for D,Phi ==> up,dn
	!	add 1 to get index for ==> dn,up OR
	!	???use the same index but switch which
	!	???of D/Phi site is going to be added at N+1 as up.
	use modmain, only: basis
	use basisstates, only: LexicoIndex
	implicit none
	integer(kind=1), intent(in) :: ib,n,k
	integer(kind=4), intent(in) :: la
	integer(kind=4), dimension(la), intent(out):: mapa
	! local
	integer(kind=1), dimension(k+1) :: set2
	integer(kind=1):: n1,n2,k1
	integer:: i

	n1 = n+1;	n2 = n+2;
	k1 = k+1;
	do i=1,la
		set2(1:k) = basis(ib)%sec(k)%sets(i,:)
		set2(k1) = n1; ! an up at n+1, dn at n+2
		mapa(i) = LexicoIndex(set2,n2,k1)
	end do
	return
	end subroutine AnMap1
!------------------------------------------
!	Map for D Phi Annihilation for channel 1,2,3,4
!------------------------------------------
	subroutine AnMap2(ib,n,k,map,la)
	!	chan 1,2,3,4
	! D&Phi sites appended to active list.
	! out: index for D,Phi ==> up,dn
	!	add 1 to get index for ==> dn,up OR
	use modmain, only: basis
	use basisstates, only: LexicoIndex
	!use lists, only: MemberQ
	implicit none
	integer(kind=1), intent(in) :: ib,n,k
	integer(kind=4), intent(in) :: la
	integer(kind=4), dimension(3,la), intent(out):: map
	! local
	integer(kind=1), dimension(k+2) :: set2
	integer(kind=1):: n1,n2,k1,k2
	integer:: i

	n1 = n+1;	n2 = n+2;
	k1 = k+1; k2 = k+2;
	do i=1,la
		! for channel 1,2
		set2(1:k) = basis(ib)%sec(k)%sets(i,:)
		set2(k1) = n1; ! an up at n+1, dn at n+2
		map(1,i) = LexicoIndex(set2(1:k1),n2,k1)
		! for channel 3: dn,dn at n+1,n+2
		map(2,i) = LexicoIndex(set2(1:k),n2,k)
		! for channel 4: up,up at n+1,n+2
		set2(k2) = n2;
		map(3,i) = LexicoIndex(set2,n2,k2)
	end do
	return
	end subroutine AnMap2
!-----------------------------------




!-----------------------------------
!	D Phi Annihilation for channel 1,2
!------------------------------------------
	subroutine DPhiAn1()
	use modmain, only: basis,hop,na,nx
	implicit none
	!	local
	integer(kind=1):: ih=5, is=1, ib1=3, ib2=5 ! see dnalist5 in modmain
	integer(kind=1)::k,n,m,m1,n1,n2,m2,i
	integer :: ntot, ind, ntot1
	integer, allocatable, dimension(:) :: map
	integer, allocatable, dimension(:):: pntr1,pntr2

	!nc=4; 
	!hop(ih)%nc = nc
	!hop(ih)%ns = ns
	if (allocated(hop(ih)%ht))deallocate(hop(ih)%ht)
	allocate(hop(ih)%ht(4,1))
	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations
	if (n==0) then
!		write(*,*) "DPhiAn1: n = 0"
	! allocate transition matrix: diagonal format
		do i=1,2
			if (allocated(hop(ih)%ht(i,is)%col))
     .				deallocate(hop(ih)%ht(i,is)%col)
			allocate(hop(ih)%ht(i,is)%col(1)) 
		end do
		hop(ih)%ht(1,is)%col(1) = 1 + 1 ! n+1 = 1 ==> ind=1 in k=1 sec
		hop(ih)%ht(2,is)%col(1) = 1 + 2 ! n+2 = 2 ==> ind=2 in k=1 sec
														! 1 added for k=0 sector
		return
	endif

	m1 = min(m,n); ! max no of up spins possible
	n1 = n+1; n2 = n+2;
	! dm = 1 (chan 1,2), 0 (chan 3), 2 (chan 4)
	m2 = min(m+1,n2); ! chan 1,2

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2))
	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! only the relevant part
	pntr2(:) = basis(ib2)%pntr(1:m2+2) 
	ntot1 = pntr1(m1+2); ! total number of basis states
	!ntot2 = pntr2(m2+2);

!	write(*,*) "DPhiAn1: pntr done .... "



	! allocate transition matrix: diagonal format
		do i=1,2
			if (allocated(hop(ih)%ht(i,is)%col))
     .				deallocate(hop(ih)%ht(i,is)%col)
			allocate(hop(ih)%ht(i,is)%col(ntot1)) 
		end do

	!	calc the matrix

	! k==0
	hop(ih)%ht(1,is)%col(1) = n1;
	hop(ih)%ht(2,is)%col(1) = n2;

!	write(*,*) "DPhiAn1: hop k=0 done .... "

	ind = 2;
	do k=1,m1,1
		ntot = pntr1(k+2)-pntr1(k+1)
		allocate(map(ntot));
		call AnMap1(ib1,n,k,map,ntot)
		map = pntr2(k+2) + map
		hop(ih)%ht(1,is)%col(ind:ind+ntot-1) = map
		hop(ih)%ht(2,is)%col(ind:ind+ntot-1) = map + 1
		ind = ind+ntot
		deallocate(map)
	end do

!	write(*,*) "DPhiAn1: hop k>0 done .... "

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiAn1
!------------------------------------------





!------------------------------------------
!	D Phi Annihilation for channel 1,2,3,4
!------------------------------------------
	subroutine DPhiAn2()
	! chan 1,2,3,4 
	use modmain, only: basis,hop,na,nx
	implicit none
	!	local
	integer(kind=1):: ih=5, is=1, ib1=3, ib2=5 ! see dnalist5 in modmain
	!	itype? n,m, etc....?
	integer(kind=1)::k,n,m,m1,n1,n2,m2,m3,m4,i
	integer :: ntot, ind, ind2, ntot1 !, ntot2, ntot3, ntot4
	integer, allocatable, dimension(:,:) :: map
	integer, allocatable, dimension(:):: pntr1,pntr2,map2



	!nc=4; 
	!hop(ih)%nc = nc
	!hop(ih)%ns = ns
	if (allocated(hop(ih)%ht))deallocate(hop(ih)%ht)
	allocate(hop(ih)%ht(4,1))
	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations
	if (n==0) then
	! allocate transition matrix: diagonal format
		do i=1,4
			if (allocated(hop(ih)%ht(i,is)%col))
     .				deallocate(hop(ih)%ht(i,is)%col)
			allocate(hop(ih)%ht(i,is)%col(1)) 
		end do
		hop(ih)%ht(1,is)%col(1) = 1 + 1 ! n+1 = 1 ==> ind=1 in k=1 sec
		hop(ih)%ht(2,is)%col(1) = 1 + 2 ! n+2 = 2 ==> ind=2 in k=1 sec
																! 1 added for k=0 sector
		hop(ih)%ht(3,is)%col(1) = 1 			! no up ==> ind=1 in k=0 sec
		hop(ih)%ht(4,is)%col(1) = 1 + 2 ! 1,2 ==> ind=1 in k=2 sec (both up)
		return
	endif

	m1 = min(m,n); ! max no of up spins possible
	n1 = n+1; n2 = n+2;
	m2 = min(m+1,n2); ! chan 1,2
	m3 = min(m,n2); ! chan 3
	m4 = min(m+2,n2); ! chan 4


	allocate(pntr1(m1+2))
	allocate(pntr2(m4+2))

	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! only the relevant part
	pntr2(:) = basis(ib2)%pntr(1:m4+2) 
	ntot1 = pntr1(m1+2); ! total number of basis states

	! allocate transition matrix: diagonal format
		do i=1,4
			if (allocated(hop(ih)%ht(i,is)%col))
     .				deallocate(hop(ih)%ht(i,is)%col)
			allocate(hop(ih)%ht(i,is)%col(ntot1)) 
		end do

	!	calc the matrix

	! k==0
	hop(ih)%ht(1,is)%col(1) = 1	+ n1; ! 1 for k=0 sec in n+2 case
	hop(ih)%ht(2,is)%col(1) = 1 + n2; ! 
	hop(ih)%ht(3,is)%col(1) = 1; ! k=0 sec in n+2 case
	hop(ih)%ht(4,is)%col(1) = pntr2(3);! k=2 sec in n+2 case
	! pntr2(3) = pntr2(2) + ind of last basis of k=2	

	! k>0 	
	ind = 2;
	do k=1,m1,1
		ntot = pntr1(k+2)-pntr1(k+1)
		allocate(map(3,ntot));
		call AnMap2(ib1,n,k,map,ntot)
		allocate(map2(ntot));
		map2(:) = pntr2(k+2) + map(1,:)
		ind2 = ind+ntot-1;
		hop(ih)%ht(1,is)%col(ind:ind2) = map2
		hop(ih)%ht(2,is)%col(ind:ind2) = map2 + 1
		hop(ih)%ht(3,is)%col(ind:ind2) = pntr2(k+1) + map(2,:)
		hop(ih)%ht(4,is)%col(ind:ind2) = pntr2(k+3) + map(3,:)
		ind = ind2+1
		deallocate(map,map2)
	end do
	
	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiAn2
!------------------------------------------

