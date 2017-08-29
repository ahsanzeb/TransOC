

!-----------------------------------
	subroutine dhopsmap1(ib,k,l,mapa,la,mapb,lb)
	use modmain, only: basis
	use lists, only: MemberQ
	implicit none
	integer(kind=1), intent(in) :: ib,k,l
	integer(kind=4), intent(in) :: la,lb
	integer(kind=4), dimension(la), intent(out):: mapa
	integer(kind=4), dimension(lb), intent(out):: mapb
	! local
	integer(kind=1), dimension(k) :: set
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
	use modmain, only: basis,hop,na,nx,dna,dnx,ibs,mapb
	implicit none

	integer(kind=1), intent(in) :: ih, is
	!	local
	integer :: itype=1
	!	itype? n,m, etc....?
	integer(kind=1)::ib,ibi,k,n,m,m1
	integer :: ntot, lat,lbt,inda,indb,la,lb
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:) :: mapa, mapbb
	integer, allocatable, dimension(:):: pntr

	!------------------------------------------	
	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible
	!------------------------------------------
	!	m1==0 case: no excitations: matrices are [[]] and [[1]];
	!------------------------------------------
	if (m1==0) then
		if (allocated(hop(ih)%ht(2,is)%row)) then
			deallocate(hop(ih)%ht(2,is)%row)
		endif
		allocate(hop(ih)%ht(2,is)%row(1))
		! m1=0: empty chan 1 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		hop(ih)%ht(2,is)%row(1) = 1;
		return
	endif
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
	lbt = sum(lbs)
	! allocate transition matrix: diagonal format
	if (allocated(hop(ih)%ht(1,is)%row)) then
		deallocate(hop(ih)%ht(1,is)%row)
		deallocate(hop(ih)%ht(2,is)%row)
	endif
	allocate(hop(ih)%ht(1,is)%row(lat)) 
	allocate(hop(ih)%ht(2,is)%row(lbt))

	!	calc the matrix
	inda = 1; indb=1;
	if (k==0) then
		hop(ih)%ht(2,is)%row(indb) = 1;
		indb = indb+1
	endif

	do k=1,m1,1
		la = las(k); lb = lbs(k); 
		allocate(mapa(la));
		allocate(mapbb(lb))
		!	calc maps
		call dhopsmap1(ib,k,is,mapa,la,mapbb,lb)
		!	assign values to transition matrices
		hop(ih)%ht(1,is)%row(inda:inda+la-1) = pntr(k) + mapa
		hop(ih)%ht(2,is)%row(indb:indb+lb-1) = pntr(k) + mapbb
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapbb)
	end do

	!write(*,*)"row2=", hop(ih)%ht(2,is)%row(1:min(10,lbt))
	
	deallocate(pntr,las,lbs)
	return
	end subroutine dhops1
!------------------------------------------










!-----------------------------------
	subroutine dhopsmap2(ib,n,k,l,mapa,la,mapb,lb)
	use modmain, only: basis
	use basisstates, only: LexicoIndex
	use lists, only: MemberQ,Drop,SortedInsert
	implicit none
	integer(kind=1), intent(in) :: ib,n,k,l
	integer(kind=4), intent(in) :: la,lb
	integer(kind=4), dimension(la,2), intent(out):: mapa
	integer(kind=4), dimension(lb,2), intent(out):: mapb

	! local
	integer(kind=1), dimension(k) :: set
	integer(kind=1), dimension(k-1) :: seta
	integer(kind=1), dimension(k+1) :: setb
	integer:: ntot,i,inda,indb
	integer(kind=1) :: k1

	
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
	use modmain, only: basis,hop,na,nx,dna,dnx,ibs,mapb
	implicit none
	integer(kind=1), intent(in) :: ih, is
	!	local
	integer :: itype=1
	!	itype? n,m, etc....?
	integer(kind=1)::ib,ibi,k,n,m,m1,m3!,k1,m2
	integer :: ntot, lat,lbt,inda,indb,la,lb
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:,:) :: mapa, mapbb
	integer, allocatable, dimension(:):: pntr

	!------------------------------------------	
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
		if (allocated(hop(ih)%ht(2,is)%row)) then
			deallocate(hop(ih)%ht(2,is)%row)
			deallocate(hop(ih)%ht(4,is)%col)
		endif
		allocate(hop(ih)%ht(2,is)%row(1))
		allocate(hop(ih)%ht(4,is)%col(1))
		! m1=0: empty chan 1 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		hop(ih)%ht(2,is)%row(1) = 1;
		hop(ih)%ht(4,is)%col(1) = is; ! LexicoIndex((/ is /),n,1)
		return
	endif
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
	lbt = 1 + sum(lbs) ! 1 for k=0 case that is only for chan 2,4
	! allocate transition matrix: diagonal format
	if (allocated(hop(ih)%ht(1,is)%row)) then
		deallocate(hop(ih)%ht(1,is)%row) ! chan 1,2 diagonal => only row
		deallocate(hop(ih)%ht(2,is)%row)
	endif
	if (allocated(hop(ih)%ht(3,is)%col)) then
		deallocate(hop(ih)%ht(3,is)%col) ! chan 3,4 same row as chan 1,2, keep col only
		deallocate(hop(ih)%ht(4,is)%col)
	endif
	allocate(hop(ih)%ht(1,is)%row(lat)) 
	allocate(hop(ih)%ht(2,is)%row(lbt))
	allocate(hop(ih)%ht(3,is)%col(lat)) 
	allocate(hop(ih)%ht(4,is)%col(lbt))

	! set nnz
	hop(ih)%ht(1,is)%nnz = lat
	hop(ih)%ht(2,is)%nnz = lbt
	hop(ih)%ht(3,is)%nnz = lat
	hop(ih)%ht(4,is)%nnz = lbt

	
	!	calc the matrix
	!	k=0, chan 2,4 set by hand
	hop(ih)%ht(2,is)%row(1) = 1;	
	hop(ih)%ht(4,is)%col(1) = is; ! LexicoIndex((/ is /),n,1)
	inda = 1; indb=2;
	do k=1,m1,1
		la = las(k); lb = lbs(k); 
		allocate(mapa(la,2));
		allocate(mapbb(lb,2));
		!	calc maps
		call dhopsmap2(ib,n,k,is,mapa,la,mapbb,lb)
		!	assign values to transition matrices
		hop(ih)%ht(1,is)%row(inda:inda+la-1) = pntr(k+1) + mapa(:,1)
		hop(ih)%ht(2,is)%row(indb:indb+lb-1) = pntr(k+1) + mapbb(:,1)
		hop(ih)%ht(3,is)%col(inda:inda+la-1) = pntr(k) + mapa(:,2) !prev sector
		hop(ih)%ht(4,is)%col(indb:indb+lb-1) = pntr(k+2) + mapbb(:,2)!next sector
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapbb)
	end do

	!write(*,*)"row2=", hop(ih)%ht(2,is)%row(1:min(10,lbt))	
	deallocate(pntr,las,lbs)
	return
	end subroutine dhops2
!------------------------------------------


