
	module Dhops
	use amplitudes
	use lists, only: MemberQ,Drop,SortedInsert,GetPosition,MemberQ4
	
	implicit none
	integer::one=1,two=2,three=3,four=4

	public :: dhops1, dhops2
	private :: dhopsmap1, dhopsmap2

	contains
!-----------------------------------
	subroutine dhopsmap1(ib,k,l,mapa,la,mapb,lb)
	use modmain, only: basis,isk
	implicit none
	integer, intent(in) :: ib,k
	integer, intent(in) :: l
	integer, intent(in) :: la,lb
	integer, dimension(la), intent(out):: mapa
	integer, dimension(lb), intent(out):: mapb
	! local
	integer(kind=isk), dimension(k) :: set
	integer:: ntot,i,inda,indb

	ntot = la+lb;
	inda = 0; indb=0;
	do i=1,ntot
		set=basis(ib)%sec(k)%sets(i,:)
		if (MemberQ(set,k,l)) then
			inda = inda +1
			mapa(inda) = i
		else
			indb = indb +1
			mapb(indb) = i
		endif
	end do

	return
	end subroutine dhopsmap1
!-----------------------------------
	subroutine dhops1(ih,is)
	! chan 1,2 only
	use modmain,only: basis,na,nx,dna,dnx,ibs,mapb,ways,Asites,isk
	implicit none

	integer, intent(in) :: ih
	integer, intent(in) :: is
	!	local
	integer :: itype=1
	integer :: l
	integer::ib,ibi,k,n,m,m1,ia,ic
	integer :: ntot, lat,lbt,inda,indb,la,lb
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:) :: mapa, mapbb
	integer, allocatable, dimension(:):: pntr,row1,row2
	double precision, allocatable:: HtUf(:,:)
	integer:: n1,n2,n3
	integer:: itl
	logical :: dhop, phop


	dhop = MemberQ4((/1,2,27,28/),4,ih)
	phop = MemberQ4((/3,4,29,30/),4,ih)

	l = ways(ih)%active(is); ! localtion of active site in lattice
	l = GetPosition(Asites,na,l); ! localtion of active site in Asites
	!------------------------------------------	
	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible
	!------------------------------------------
	!	m1==0 case: no excitations: matrices are [[]] and [[1]];
	!------------------------------------------
	if (m1==0) then
		!allocate(row2(1))
		! m1=0: empty chan 1 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		!row2(1) = 1;
		!n3=1;
		! swap channel 1,2 for PhiHops
		if(dhop) then
			call CalAmp(ih,2,is,(/1/),1,1,'multiplyd') ! ic=2
		else	if(phop) then
			call CalAmp(ih,1,is,(/1/),1,1,'multiplyd') ! ic=1
		else
			write(*,*)"Error(dhops1): stopping!";
			stop
		endif
		! set rates for ic=1 to 0 ????
	else
	!------------------------------------------
	! m1 > 0 case
	!------------------------------------------
	! pointers for start index
	allocate(pntr(m1+2));
	! ib: itype ===> which of 5 N case?
	ibi = ibs(itype);
	ib = mapb%map(ibi);	
	pntr(:) = basis(ib)%pntr(1:m1+2) ! only the relevant part
	! dimensions of maps for diff k
	allocate(las(m1))
	allocate(lbs(m1))
	las(:) = 0;	lbs(:)=0;
	do k=1,m1
		ntot = pntr(k+2)-pntr(k+1)
		las(k) = k*ntot/n
		lbs(k) = (n-k)*ntot/n
	end do
	! dimensions of full transition matrices
	lat = sum(las)
	lbt = sum(lbs) + 1; ! +1 for k=0 term
	allocate(row1(lat)) 
	allocate(row2(lbt))

	!	calc the matrix
	inda = 1; indb=1;

	!k==0
	row2(indb) = 1;
	indb = indb+1

	do k=1,m1,1
		la = las(k); lb = lbs(k); 
		allocate(mapa(la));
		allocate(mapbb(lb))
		!	calc maps
		call dhopsmap1(ib,k,l,mapa,la,mapbb,lb)
		!	assign values to transition matrices
		row1(inda:inda+la-1) = pntr(k+1) + mapa
		row2(indb:indb+lb-1) = pntr(k+1) + mapbb
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapbb)
	end do

	

	n3=pntr(m1+2) ! dim of initial hilbert space
	deallocate(pntr,las,lbs)

	!write(*,*) "dhops: row1 = ",row1
	!write(*,*) "dhops: row2 = ",row2
	! calculate transition amplitudes
	! swap channel 1,2 for PhiHops
	if(dhop) then
		call CalAmp(ih,1,is,row1,lat,n3,'multiplyd') ! ic=1
		call CalAmp(ih,2,is,row2,lbt,n3,'multiplyd') ! ic=2
	else	if(phop) then
		call CalAmp(ih,2,is,row1,lat,n3,'multiplyd') ! ic=2
		call CalAmp(ih,1,is,row2,lbt,n3,'multiplyd') ! ic=1
	else
			write(*,*)"Error(dhops1): stopping!";
			stop
	endif


	endif
	
	return
	end subroutine dhops1
!------------------------------------------










!-----------------------------------
	subroutine dhopsmap2(ib,n,k,l,mapa,la,mapb,lb)
	use modmain, only: basis,isk
	use basisstates, only: LexicoIndex
	
	implicit none
	integer, intent(in) ::l
	integer, intent(in) :: ib,n,k
	integer, intent(in) :: la,lb
	integer, dimension(la,2), intent(out):: mapa
	integer, dimension(lb,2), intent(out):: mapb

	! local
	integer(kind=isk), dimension(k) :: set
	integer(kind=isk), dimension(k-1) :: seta
	integer(kind=isk), dimension(k+1) :: setb
	integer:: ntot,i,inda,indb,k1
	
	ntot = la+lb;
	inda = 0; indb=0;
	do i=1,ntot
		set=basis(ib)%sec(k)%sets(i,:)
		!write(*,*) "LexicoIndex = "	,LexicoIndex(set,n,k)	
		if (MemberQ(set,k,l)) then
			inda = inda +1
			mapa(inda,1) = i
			call Drop(set,k,l,seta)
			k1 = k-1
			mapa(inda,2) = LexicoIndex(seta,n,k1)	
		else
			indb = indb +1
			mapb(indb,1) = i
			! set is already sorted so use SortedInsert for append+sort, more efficient
			call SortedInsert(set,k,l,setb)
			k1 = k+1
			mapb(indb,2) = LexicoIndex(setb,n,k1)		
		endif
	end do

	return
	end subroutine dhopsmap2
!-----------------------------------





!-----------------------------------
	subroutine dhops2(ih,is)
	!	channel 1,2 have diagonal Ht, save row only
	!	channel 3,4 have the same row as 1,2, save col only
	use modmain,only: basis,na,nx,dna,dnx,ibs,mapb,ways,Asites,isk
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih
	!	local
	integer :: itype=1
	integer:: l
	integer::ib,ibi,k,n,m,m1,m3,ic
	integer :: ntot, lat,lbt,inda,indb,la,lb,ia,n1,n2,n3
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:,:) :: mapa, mapbb
	integer, allocatable, dimension(:):: pntr
	integer, allocatable, dimension(:):: row1,row2,col3,col4
	double precision, allocatable:: HtUf(:,:)
	integer :: itl
	logical :: dhop, phop

	dhop = MemberQ4((/1,2,27,28/),4,ih)
	phop = MemberQ4((/3,4,29,30/),4,ih)


	l = ways(ih)%active(is); ! localtion of active site in lattice
	l = GetPosition(Asites,na,l); ! localtion of active site in Asites

	!-------------------------------------------------------
	! calculate Transition matrix
	!-------------------------------------------------------
	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! chan 1,2, max no of up spins possible
	!m2 = min(m-1,n); ! chan 3, max no of up spins possible
	m3 = min(m+1,n); ! chan 4, max no of up spins possible
	!------------------------------------------
	!	m1==0 case: no excitations: matrices are [[]] and [[1]];
	!------------------------------------------
	if (m1==0) then
		allocate(row2(1))
		allocate(col4(1))
		! m1=0: empty chan 1,3 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		row2(1) = 1;
		col4(1) = 1 + is; ! added 1 for k=0 sec; LexicoIndex((/ is /),n,1)
		n3=1;
		lbt=1;
		! calculate transition amplitudes

		! swap channel 1,2 for PhiHops
		if(dhop) then
			call CalAmp(ih,2,is,row2,lbt,n3,'multiplyd') ! ic=2
		else	if(phop) then
			call CalAmp(ih,1,is,row2,lbt,n3,'multiplyd') ! ic=1
		else
			write(*,*)"Error(dhops2): stopping!";
			stop
		endif
		call CalAmp0(ih,4,is,row2,lbt,n3,col4) ! ic=4
	
		! set rates for ic=1 to 0 ????
	else
	!------------------------------------------
	! m1 > 0 case
	!------------------------------------------
	! ib: itype ===> which of 5 N case?
	ibi = ibs(itype);
	ib = mapb%map(ibi);

	! pointers for start index
	allocate(pntr(m3+2));
	pntr(:) = basis(ib)%pntr(1:m3+2) ! only the relevant part
	
	! dimensions of maps for diff k
	allocate(las(m1))
	allocate(lbs(m1))
	las(:) = 0;	lbs(:)=0;
	do k=1,m1
		ntot = pntr(k+2)-pntr(k+1)
		las(k) = k*ntot/n
		lbs(k) = (n-k)*ntot/n
	end do
	! dimensions of full transition matrices
	lat = sum(las)
	lbt = 1 + sum(lbs) ! 1 for k=0 case that is only for chan 2,4 of Dhops; 1,4 of Phihops

	! allocate transition matrix: diagonal format
	allocate(row1(lat)) 
	allocate(row2(lbt))
	allocate(col3(lat)) 
	allocate(col4(lbt))
	
	!	calc the matrix
	!	k=0, chan 2,4 set by hand
	row2(1) = 1;	
	col4(1) = 1 + is; ! LexicoIndex((/ is /),n,1); 1 for k=0 basis
	inda = 1; indb=2;
	do k=1,m1,1
		la = las(k); lb = lbs(k);
		allocate(mapa(la,2));
		allocate(mapbb(lb,2));
		!	calc maps
		call dhopsmap2(ib,n,k,l,mapa,la,mapbb,lb)
		!	assign values to transition matrices
		row1(inda:inda+la-1) = pntr(k+1) + mapa(:,1)
		row2(indb:indb+lb-1) = pntr(k+1) + mapbb(:,1)
		col3(inda:inda+la-1) = pntr(k) + mapa(:,2) !prev sector
		col4(indb:indb+lb-1) = pntr(k+2) + mapbb(:,2)!next sector
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapbb)
	end do


	n3=pntr(m1+2) ! dim of initial hilbert space
	deallocate(pntr,las,lbs)

	! calculate transition amplitudes
	! swap channel 1,2 for PhiHops
	if(dhop) then
		call CalAmp(ih,1,is,row1,lat,n3,'multiplyd') ! ic=1
		call CalAmp(ih,2,is,row2,lbt,n3,'multiplyd') ! ic=2
	else	if(phop) then
		call CalAmp(ih,2,is,row1,lat,n3,'multiplyd') ! ic=2
		call CalAmp(ih,1,is,row2,lbt,n3,'multiplyd') ! ic=1
	else
			write(*,*)"Error(dhops2): stopping!";
			stop
	endif
	call CalAmp0(ih,three,is,row1,lat,n3,col3) ! ic=3
	call CalAmp0(ih,four,is,row2,lbt,n3,col4) ! ic=4


	endif ! m==0

	return
	end subroutine dhops2
!------------------------------------------

	end module dhops
