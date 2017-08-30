

!------------------------------------------
!	Map for D Phi creation at site l1,l2 for channel 1,2
!------------------------------------------
! starting with (N-2, k-1), add desired up/dn
! or dn/up to make (N, k) 
	subroutine CreatMap1(n,k,l1,l2,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb
	use basisstates, only: LexicoIndex,Shift
	implicit none
	integer(kind=1), intent(in) :: n,k,l1,l2
	integer(kind=4), intent(in) :: ntot
	integer(kind=4), dimension(2,ntot), intent(out):: map
	! local
	integer(kind=1), dimension(k+1) :: set2
	integer(kind=1):: n2,k1
	integer(kind=1):: ibi=1, ib ! ib1=1: dN=-2 for starting basis
	integer(kind=1), dimension(k-1):: set
	
	integer:: i

	ib= mapb%map(ibi);

	n2 = n+2; k1 = k-1;
	do i=1,la
		set = basis(ib)%sec(k1)%sets(i,:)
		! add an up site at n1 position and a down at n2 position;
		! shift n1+ and n2+ labels by 1 or 2 as they should be.
		!	shift by +1 twice for second case means +2 total as desired
		call Shift(set,k1,l1) ! shift all elem by +1 from l1 onwards
		call Shift(set,k1,l2) ! shift all elem by +1 from l2 onwards
		! Adding up/down  at l1/l2
		call SortedInsert(set,k1,l1,set2)
		map(1,i) = LexicoIndex(set2,n2,k)
		! Adding down/up  at l1/l2
		call SortedInsert(set,k1,l2,set2)
		map(2,i) = LexicoIndex(set2,n2,k)
	end do
	return
	end subroutine CreatMap1

!------------------------------------------
!	Map for D Phi creation at site l1,l2 for channel 3
!------------------------------------------
! starting with (N-2, k-2), add up,up at l1,l2 to make (N, k) 
	subroutine CreatMap3(n,k,l1,l2,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb
	use basisstates, only: LexicoIndex,Shift
	implicit none
	integer(kind=1), intent(in) :: n,k,l1,l2
	integer(kind=4), intent(in) :: ntot
	integer(kind=4), dimension(ntot), intent(out):: map
	! local
	integer(kind=1):: n2,k1
	integer(kind=1):: ibi=1,ib ! ib1=1: dN=-2 for starting basis
	integer(kind=1), dimension(k-2):: set
	integer(kind=1), dimension(k-1) :: set2
	integer(kind=1), dimension(k):: set3	
	integer:: i

	ib = mapb%map(ibi);

	n2 = n+2; k1 = k-1; k2=k-2;
	do i=1,ntot
		set = basis(ib)%sec(k2)%sets(i,:)
		! add an up site at n1 position and a down at n2 position;
		! shift n1+ and n2+ labels by 1 or 2 as they should be.
		!	shift by +1 twice for second case means +2 total as desired
		call Shift(set,k2,l1) ! shift all elem by +1 from l1 onwards
		call Shift(set,k2,l2) ! shift all elem by +1 from l2 onwards
		! Adding up  at l1
		call SortedInsert(set,k2,l1,set2)
		! Adding up  at l2
		call SortedInsert(set2,k1,l2,set3)
		map(i) = LexicoIndex(set3,n2,k)
	end do
	return
	end subroutine CreatMap3
!------------------------------------------
!	Map for D Phi creation at site l1,l2 for channel 4
!------------------------------------------
! starting with (N-2, k), add dn,dn at l1,l2 to make (N, k) 
	subroutine CreatMap4(n,k,l1,l2,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb
	use basisstates, only: LexicoIndex,Shift
	implicit none
	integer(kind=1), intent(in) :: n,k,l1,l2
	integer(kind=4), intent(in) :: ntot
	integer(kind=4), dimension(ntot), intent(out):: map
	! local
	integer(kind=1):: n2,k1
	integer(kind=1):: ib=1,ib ! ib1=1: dN=-2 for starting basis
	integer(kind=1), dimension(k):: set	
	integer:: i

	ib = mapb%map(ibi);
	
	n2 = n+2;
	do i=1,ntot
		set = basis(ib)%sec(k)%sets(i,:)
		! add an up site at n1 position and a down at n2 position;
		! shift n1+ and n2+ labels by 1 or 2 as they should be.
		!	shift by +1 twice for second case means +2 total as desired
		call Shift(set,k,l1) ! shift all elem by +1 from l1 onwards
		call Shift(set,k,l2) ! shift all elem by +1 from l2 onwards
		map(i) = LexicoIndex(set,n2,k)
	end do
	return
	end subroutine CreatMap4
!------------------------------------------



































!------------------------------------------
!	D Phi Annihilation for channel 1,2
! up,dn become D,Phi
!------------------------------------------
	subroutine DPhiCreat1(is,l1,l2)
	use modmain, only: basis,hop,na,nx
	implicit none
	integer(kind=1), intent(in):: is,l1,l2 
	!	local
	integer(kind=1):: ih=7, ib1=3, ib2=1 ! see dnalist5 in modmain
	integer(kind=1)::k,n,m,m1,m2,i,m3
	integer :: ntot, ind, ntot1,nnz,n1,n2
	integer, allocatable, dimension(:,:) :: map,row
	integer, allocatable, dimension(:):: pntr1,pntr2,col

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n-2,m-1);
	m3 = min(m1,n-1); ! at least 1 dn required

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2(:) = basis(ib2)%pntr(1:m2+2) ! final
	ntot1 = pntr1(m1+2); ! total number of basis states

	! nnz total number of non zero elements
	nnz = pntr2(m3+2) - 1; ! up to pntr2(m3+2), -1 to exclude k=0 
	! allocate transition matrix: coo format
	!hop(ih)%ht(1,is)%nnz = nnz;
	!hop(ih)%ht(2,is)%nnz = nnz;

	allocate(row(nnz,2)) 
	allocate(col(nnz)) 

	!	calc the matrix
	!	k=1
	ind = 1; ntot=1
	row(1,1) = l1 ! l1 up in k=1 sec
	row(1,2) = l2 ! l2 up in k=1 sec
	! final basis
	col(1,1) = (/ 1 /) ! k=0 sec of final

	!	k>1
	ind = 2;
	do k=2,m3,1
		ntot = pntr2(k+1)-pntr2(k) ! k-1 sec of final
		allocate(map(2,ntot));
		call CreatMap1(n,k,l1,l2,map,ntot)
		! initial basis
		map = pntr1(k+1) + map
		row(ind:ind+ntot-1,1) = map(1,:)
		row(ind:ind+ntot-1,2) = map(2,:)
		! final basis
		col(ind:ind+ntot-1,1) = &
			& (/ (i,i=pntr2(k)+1,pntr2(k+1),1 ) /) ! k-1 sec of final
		ind = ind+ntot
		deallocate(map)
	end do

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	call CalAmp(ih,1,is,row(:,1),nnz,n3,"multiply",col) ! ic=1
	call CalAmp(ih,2,is,row(:,2),nnz,n3,"multiply",col) ! ic=2
	!---------------------------------

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiCreat1
!------------------------------------------










!------------------------------------------
!	D Phi Annihilation for channel 3
! up,up become D,Phi
!------------------------------------------
	subroutine DPhiCreat3(is,l1,l2)
	use modmain, only: basis,hop,na,nx
	implicit none
	integer(kind=1), intent(in):: is,l1,l2 
	!	local
	integer(kind=1):: ih=7, ib1=3, ib2=1 ! see dnalist5 in modmain
	integer(kind=1)::k,n,m,m1,m2,i
	integer :: ntot, ind, nnz,n1,n2
	integer, allocatable, dimension(:) :: map,row,col
	integer, allocatable, dimension(:):: pntr1,pntr2

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n-2,m-2);

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2(:) = basis(ib2)%pntr(1:m2+2) ! final
	!ntot1 = pntr1(m1+2); ! total number of basis states

	! nnz total number fo non zero elements
	nnz = pntr2(m3+2); ! from k=0 sec up to k=m3; pntr2(m3+2) 
	! allocate transition matrix: coo format
	!hop(ih)%ht(3,is)%nnz = nnz;

	allocate(row(nnz)) 
	allocate(col(nnz)) 

	!	2 >= k <= m3
	ind = 1;
	do k=2,m2,1
		ntot = pntr2(k)-pntr2(k-1) ! k-2 sec of final
		allocate(map(ntot));
		call CreatMap3(n,k,l1,l2,map,ntot)
		! initial basis
		map = pntr1(k+1) + map
		row(ind:ind+ntot-1) = map
		! final basis
		col(ind:ind+ntot-1) = &
			& (/ (i,i=pntr2(k-1)+1,pntr2(k),1 ) /) ! k-2 sec of final
		ind = ind+ntot
		deallocate(map)
	end do

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	call CalAmp(ih,3,is,row,nnz,n3,"multiply",col) ! ic=3
	!---------------------------------


	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiCreat3
!------------------------------------------



!------------------------------------------
!	D Phi Annihilation for channel 4
! dn,dn become D,Phi
!------------------------------------------
	subroutine DPhiCreat4(is,l1,l2)
	use modmain, only: basis,hop,na,nx
	implicit none
	integer(kind=1), intent(in):: is,l1,l2 
	!	local
	integer(kind=1):: ih=7, ib1=3, ib2=1 ! see dnalist5 in modmain
	integer(kind=1)::k,n,m,m1,m2,i,m3
	integer :: ntot, ind, nnz,n1,n2
	integer, allocatable, dimension(:) :: map,row,col
	integer, allocatable, dimension(:):: pntr1,pntr2

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n-2,m); ! m3=min(n-2,m)=m2 at least two up required

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2(:) = basis(ib2)%pntr(1:m2+2) ! final
	!ntot1 = pntr1(m1+2); ! total number of basis states

	! nnz total number fo non zero elements
	nnz = pntr2(m2+2);
	! allocate transition matrix: coo format
	!hop(ih)%ht(4,is)%nnz = nnz;

	allocate(row(nnz)) 
	allocate(col(nnz)) 

	! final | initial k=0 case: by hand for efficiency !
	! no actual need to seperate from the k loop below
	hop(ih)%ht(4,is)%row(1) = (/1/)
	hop(ih)%ht(4,is)%col(1) = (/1/)

	!	1 >= k <= m3=m2 cases:
	ind = 2;
	do k=1,m2,1
		ntot = pntr2(k+2)-pntr2(k+1) ! k sec of final
		allocate(map(ntot));
		call CreatMap4(n,k,l1,l2,map,ntot)
		! initial basis
		map = pntr1(k+1) + map
		hop(ih)%ht(4,is)%row(ind:ind+ntot-1) = map
		! final basis
		hop(ih)%ht(4,is)%col(ind:ind+ntot-1) = &
			& (/ (i,i=pntr2(k+1)+1,pntr2(k+2),1 ) /) ! k sec of final
		ind = ind+ntot
		deallocate(map)
	end do

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	call CalAmp(ih,4,is,row,nnz,n3,"multiply",col) ! ic=4
	!---------------------------------

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiCreat4
!------------------------------------------



