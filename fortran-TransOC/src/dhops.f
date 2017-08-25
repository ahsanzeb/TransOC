

!-----------------------------------
	subroutine dhops12map(ib,k,l,mapa,la,mapb,lb)
	use modmain, only: basis
	use lists, only: MemberQ
	implicit none
	integer(kind=1), intent(in) :: ib,k,l
	integer(kind=4), intent(in) :: la,lb
	integer(kind=4), dimension(la), intent(out):: mapa
	integer(kind=4), dimension(lb), intent(out):: mapb
	! local
	integer(kind=1), dimension(k) :: set
	integer:: ntot

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
	end subroutine dhops12map
!-----------------------------------
	subroutine dhops12(ih,is)
	use modmain, only: basis,hop,na,nx,dna,dnx,nainds
	implicit none

	integer, intent(in) :: ih, is
	!	local
	integer :: itype=1
	!	itype? n,m, etc....?

	integer :: ntot, lat,lbt,inda,indb,n,m,m1
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:) :: mapa, mapb
	integer, allocatable, dimension(:):: pntr

	! NOTE: 	hop(ih)%chan(:)%site(:) dimensions of chan/sites allocated elsewhere 

	!------------------------------------------	
	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible
	!------------------------------------------
	!	m1==0 case: no excitations: matrices are [[]] and [[1]];
	!------------------------------------------
	if (m1==0) then
		if allocated(hop(ih)%chan(2)%site(is)%row) then
			deallocate(hop(ih)%chan(2)%site(is)%row)
		endif
		allocate(hop(ih)%chan(2)%site(is)%row(1))
		! m1=0: empty chan 1 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		hop(ih)%chan(2)%site(is)%row(indb:indb+1) = 1;
	return
	endif
	!------------------------------------------
	! m1 > 0 case
	!------------------------------------------
	! pointers for start index
	allocate(pntr(m1+2));
	! ib: itype ===> which of 5 N case?
	ib = nainds(itype);
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
	if allocated(hop(ih)%chan(1)%site(is)%row) then
		deallocate(hop(ih)%chan(1)%site(is)%row)
		deallocate(hop(ih)%chan(2)%site(is)%row)
	endif
	allocate(hop(ih)%chan(1)%site(is)%row(lat)) 
	allocate(hop(ih)%chan(2)%site(is)%row(lbt))

	!	calc the matrix
	inda = 1; indb=1;
	if (k==0) then
		hop(ih)%chan(1)%site(is)%row(indb:indb+1) = 1;
		indb = indb+1
	endif

	do k=1,m1,1
		la = las(k); lb = lbs(k); 
		allocate(mapa(la));
		allocate(mapb(lb))
		!	calc maps
		call dhops12map(ib,k,l,mapa,la,mapb,lb)
		!	assign values to transition matrices
		hop(ih)%chan(1)%site(is)%row(inda:inda+la) = pntr(k) + mapa(:)
		hop(ih)%chan(2)%site(is)%row(indb:indb+lb) = pntr(k) + mapb(:)
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapb)
	end do

	deallocate(pntr,las,lbs)
	return
	end subroutine dhops12
!------------------------------------------

	
