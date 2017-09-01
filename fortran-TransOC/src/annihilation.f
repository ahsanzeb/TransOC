
	module Annihilation
	implicit none

	public::DPhiAn1,DPhiAn2
	private::AnMap1,AnMap2

	contains
!------------------------------------------
!	Map for D Phi Annihilation for channel 1,2
!------------------------------------------
	subroutine AnMap1(ibl,n,k,mapa,la)
	!	chan 1,2
	! D&Phi sites appended to active list.
	! out: index for D,Phi ==> up,dn
	!	add 1 to get index for ==> dn,up OR
	!	???use the same index but switch which
	!	???of D/Phi site is going to be added at N+1 as up.
	use modmain, only: basis
	use basisstates, only: LexicoIndex
	implicit none
	integer(kind=1), intent(in) :: ibl,n,k
	integer(kind=4), intent(in) :: la
	integer(kind=4), dimension(la), intent(out):: mapa
	! local
	integer(kind=1), dimension(k+1) :: set2
	integer(kind=1):: n1,n2,k1
	integer:: i

	n1 = n+1;	n2 = n+2;
	k1 = k+1;
	do i=1,la
		set2(1:k) = basis(ibl)%sec(k)%sets(i,:)
		set2(k1) = n1; ! an up at n+1, dn at n+2
		mapa(i) = LexicoIndex(set2,n2,k1)
	end do
	return
	end subroutine AnMap1
!------------------------------------------
!	Map for D Phi Annihilation for channel 1,2,3,4
!------------------------------------------
	subroutine AnMap2(ibl,n,k,map,la)
	!	chan 1,2,3,4
	! D&Phi sites appended to active list.
	! out: index for D,Phi ==> up,dn
	!	add 1 to get index for ==> dn,up OR
	use modmain, only: basis
	use basisstates, only: LexicoIndex
	!use lists, only: MemberQ
	implicit none
	integer(kind=1), intent(in) :: ibl,n,k
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
		set2(1:k) = basis(ibl)%sec(k)%sets(i,:)
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
	use modmain, only: basis,na,nx,mapb
	implicit none
	!	local
	integer(kind=1):: ih=5, is=1, ib1i=3, ib2i=5 ! see dnalist5 in modmain
	integer(kind=1):: ib1,ib2
	integer(kind=1)::k,n,m,m1,m2,i
	integer :: ntot, ind, ntot1,n1,n2,n3
	integer, allocatable, dimension(:) :: map,col1,col2
	integer, allocatable, dimension(:):: pntr1,pntr2


	ib1= mapb%map(ib1i);
	ib2= mapb%map(ib2i);


	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations
	if (n==0) then
!		write(*,*) "DPhiAn1: n = 0"
	! allocate transition matrix: diagonal format
		allocate(col1(1)) 
		!allocate(col2(1)) 
		col1(1) = 1 + 1 ! n+1 = 1 ==> ind=1 in k=1 sec
		!col2(1) = 1 + 2 ! n+2 = 2 ==> ind=2 in k=1 sec
														! 1 added for k=0 sector
	else

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
	write(*,*) "=======>>>>>> pntr2(m2+2)=",pntr2(m2+2)
!	write(*,*) "DPhiAn1: pntr done .... "

	! allocate transition matrix: diagonal format
	allocate(col1(ntot1)) 
	allocate(col2(ntot1)) 

	!	calc the matrix

	! k==0
	col2(1) = n1;
	!col2(1) = n2;

!	write(*,*) "DPhiAn1: hop k=0 done .... "

	ind = 2;
	do k=1,m1,1
		ntot = pntr1(k+2)-pntr1(k+1)
		allocate(map(ntot));
		call AnMap1(ib1,n,k,map,ntot)
		map = pntr2(k+2) + map ! shift k by 1; final k+1 up
		col1(ind:ind+ntot-1) = map
		!col2(ind:ind+ntot-1) = map + 1
		ind = ind+ntot
		deallocate(map)
	end do

	endif

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	call CalAmp(ih,1,is,col1,n3,n3,"multiplydc") ! ic=1
	call CalAmp(ih,2,is,col1+1,n3,n3,"multiplydc") ! ic=2
	!---------------------------------

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiAn1
!------------------------------------------





!------------------------------------------
!	D Phi Annihilation for channel 1,2,3,4
!------------------------------------------
	subroutine DPhiAn2()
	! chan 1,2,3,4 
	use modmain, only: basis,na,nx,mapb
	implicit none
	!	local
	integer(kind=1):: ih=5, is=1, ib1i=3, ib2i=5 ! see dnalist5 in modmain
	!	itype? n,m, etc....?
	integer(kind=1):: ib1,ib2
	integer(kind=1)::k,n,m,m1,m2,m3,m4,i,ic
	integer :: ntot, ind, ind2, ntot1,n1,n2,n3 !, ntot2, ntot3, ntot4
	integer, allocatable, dimension(:,:) :: map,col
	integer, allocatable, dimension(:):: pntr1,pntr2,map2

	ib1= mapb%map(ib1i);
	ib2=mapb%map(ib2i);

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations
	!--------------- n=0 --------------------	
	if (n==0) then
	! allocate transition matrix: diagonal format
		!hop(ih)%ht(1,is)%nnz = 1;
		allocate(col(1,4)) 
		col(1,1) = 1 + 1 ! n+1 = 1 ==> ind=1 in k=1 sec
		col(1,2) = 1 + 2 ! n+2 = 2 ==> ind=2 in k=1 sec,
		! 1 added for k=0 sector
		col(1,3) = 1 			! no up ==> ind=1 in k=0 sec
		col(1,4) = 1 + 2 ! 1,2 ==> ind=1 in k=2 sec (both up)
	else
	!----------------- n>0 -----------------	
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
	!hop(ih)%ht(1,is)%nnz = ntot1;
	allocate(col(ntot1,4)) 
	!	calc the matrix
	! k==0
	col(1,1) = 1	+ n1; ! 1 for k=0 sec in n+2 case
	col(1,2) = 1 + n2; ! 
	col(1,3) = 1; ! k=0 sec in n+2 case
	col(1,4) = pntr2(3);! k=2 sec in n+2 case
	! pntr2(3) = pntr2(2) + ind of last basis of k=2	

	! k>0 	
	ind = 2;
	do k=1,m1,1
		ntot = pntr1(k+2)-pntr1(k+1)
		allocate(map(3,ntot));
		call AnMap2(ib1,n,k,map,ntot)
		allocate(map2(ntot));
		map2(:) = pntr2(k+2) + map(1,:) ! final k+1
		ind2 = ind+ntot-1;
		col(ind:ind2,1) = map2
		col(ind:ind2,2) = map2 + 1
		col(ind:ind2,3) = pntr2(k+1) + map(2,:) ! final k
		col(ind:ind2,4) = pntr2(k+3) + map(3,:) ! final k+2
		ind = ind2+1
		deallocate(map,map2)
	end do

	endif ! n==0
	!--------------------------------------	


	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	do ic=1,4
		call CalAmp(ih,ic,is,col(:,ic),n3,n3,"multiplydc") ! ic=1
	end do
	!call CalAmp(ih,1,is,col(:,1),n3,n3,"multiplydc") ! ic=1
	!call CalAmp(ih,2,is,col(:,2),n3,n3,"multiplydc") ! ic=2
	!call CalAmp(ih,3,is,col(:,3),n3,n3,"multiplydc") ! ic=2
	!call CalAmp(ih,4,is,col(:,4),n3,n3,"multiplydc") ! ic=2
	!---------------------------------

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiAn2
!------------------------------------------


	end module Annihilation
