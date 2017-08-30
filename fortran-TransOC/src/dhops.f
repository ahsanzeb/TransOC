

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
	use modmain, only: basis,na,nx,dna,dnx,ibs,mapb,
     .			qt,hspace, psi,itypes,maph,mapt
	implicit none

	integer(kind=1), intent(in) :: ih, is
	!	local
	integer :: itype=1
	!	itype? n,m, etc....?
	integer(kind=1)::ib,ibi,k,n,m,m1,ia
	integer :: ntot, lat,lbt,inda,indb,la,lb
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:) :: mapa, mapbb
	integer, allocatable, dimension(:):: pntr,row1,row2
	double precision, allocatable:: HtUf(:,:)
	integer:: n1,n2,n3
	integer:: ic,itl

	!------------------------------------------	
	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible
	!------------------------------------------
	!	m1==0 case: no excitations: matrices are [[]] and [[1]];
	!------------------------------------------
	if (m1==0) then
		allocate(row2(1))
		! m1=0: empty chan 1 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		row2(1) = 1;
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
	allocate(row1(lat)) 
	allocate(row2(lbt))

	!	calc the matrix
	inda = 1; indb=1;
	if (k==0) then
		row2(indb) = 1;
		indb = indb+1
	endif

	do k=1,m1,1
		la = las(k); lb = lbs(k); 
		allocate(mapa(la));
		allocate(mapbb(lb))
		!	calc maps
		call dhopsmap1(ib,k,is,mapa,la,mapbb,lb)
		!	assign values to transition matrices
		row1(inda:inda+la-1) = pntr(k) + mapa
		row2(indb:indb+lb-1) = pntr(k) + mapbb
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapbb)
	end do

	! both channel 1,2 same dimensions so loop
	! calculate transition amplitudes
	do ic=1,2
		!itype = itypes(1,ic);
		itl = mapt%map(itype); ! location of final hilber space
		n1=hspace(itl)%n1;
		n2=hspace(itl)%n2;
		n3=pntr(m1+2) ! dim of initial hilbert space
		allocate(HtUf(n3,n2))	
		ia = maph(ih,ic); ! location of amplitudes
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		if(ic==1) then
			call multiplyd(row1,lat,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		else
			call multiplyd(row2,lbt,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		endif
		!if(allocated(qt(ia)%cs(1,is)%amp))deallocate(qt(ia)%cs(1,is)%amp)
		allocate(qt(ia)%cs(ic,is)%amp(n2))
		! set the size for possible future use
		!	DELETE THIS namp VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%namp = n2
		! multiply psi with HtUf to get amplitudes
		! psi should be a row vector; shape = 1 x n3
		qt(ia)%cs(ic,is)%amp = matmul(psi,HtUf); ! both input dense
	deallocate(HtUf)
	end do

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
	use modmain, only: basis,na,nx,dna,maph,mapt,
     .			dnx,ibs,mapb,qt,hspace, psi,itypes
	implicit none
	integer(kind=1), intent(in) :: ih, is
	!	local
	integer :: itype=1
	!	itype? n,m, etc....?
	integer(kind=1)::ib,ibi,k,n,m,m1,m3!,k1,m2
	integer :: ntot, lat,lbt,inda,indb,la,lb,ia,n1,n2,n3
	integer, allocatable, dimension(:) :: las, lbs
	integer, allocatable, dimension(:,:) :: mapa, mapbb
	integer, allocatable, dimension(:):: pntr
	integer, allocatable, dimension(:):: row1,row2,col3,col4
	double precision, allocatable:: HtUf(:,:)
	integer :: ic,itl

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
		! m1=0: empty chan 1 matrix ==> 0 transition matrix,
		! handle it when calc amplitudes??????????
		row2(1) = 1;
		col4(1) = is; ! LexicoIndex((/ is /),n,1)
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
	lbt = 1 + sum(lbs) ! 1 for k=0 case that is only for chan 2,4

	! allocate transition matrix: diagonal format
	allocate(row1(lat)) 
	allocate(row2(lbt))
	allocate(col3(lat)) 
	allocate(col4(lbt))
	
	!	calc the matrix
	!	k=0, chan 2,4 set by hand
	row2(1) = 1;	
	col4(1) = is; ! LexicoIndex((/ is /),n,1)
	inda = 1; indb=2;
	do k=1,m1,1
		la = las(k); lb = lbs(k);
		allocate(mapa(la,2));
		allocate(mapbb(lb,2));
		!	calc maps
		call dhopsmap2(ib,n,k,is,mapa,la,mapbb,lb)
		!	assign values to transition matrices
		row1(inda:inda+la-1) = pntr(k+1) + mapa(:,1)
		row2(indb:indb+lb-1) = pntr(k+1) + mapbb(:,1)
		col3(inda:inda+la-1) = pntr(k) + mapa(:,2) !prev sector
		col4(indb:indb+lb-1) = pntr(k+2) + mapbb(:,2)!next sector
		inda = inda+la; indb = indb+lb
		deallocate(mapa,mapbb)
	end do

	endif ! m==0

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	do ic=1,4
		itl = mapt%map(itypes(1,ic)) ! location of final hilber space
		n1=hspace(itl)%n1 ! dim of final hilbert space
		n2=hspace(itl)%n2
		n3=pntr(m1+2) ! dim of initial hilbert space
	
		allocate(HtUf(n3,n2))
		ia = maph(ih,ic); ! location of amplitudes
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		select case(ic)
			case (1)
				call multiplyd(row1,lat,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
			case (2)
				call multiplyd(row2,lbt,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
			case (3)
				call multiply(row1,col3,lat,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
			case (4)
				call multiply(row2,col4,lbt,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		end select

		!if(allocated(qt(ia)%cs(1,is)%amp))deallocate(qt(ia)%cs(1,is)%amp)
		allocate(qt(ia)%cs(ic,is)%amp(n2))
		! set the size for possible future use
		!	DELETE THIS namp VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%namp = n2
		! multiply psi with HtUf to get amplitudes
		! psi should be a row vector; shape = 1 x n3
		qt(ia)%cs(ic,is)%amp = matmul(psi,HtUf); ! both input dense
		deallocate(HtUf)
	end do
	!---------------------------------

	deallocate(pntr,las,lbs)
	return
	end subroutine dhops2
!------------------------------------------


