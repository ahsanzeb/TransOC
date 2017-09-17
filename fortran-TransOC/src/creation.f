	module Creation
	use amplitudes
	use lists, only: SortedInsert,GetPosition
	use basisstates, only: LexicoIndex,Shift
	implicit none
	integer::one=1,two=2,three=3,four=4

	public::DPhiCreat1, DPhiCreat3, DPhiCreat4
	private::CreatMap1,CreatMap3,CreatMap4
	private:: ActiveOrder
	
	contains
!------------------------------------------
!	Map for D Phi creation at site l1,l2 for channel 1,2
!------------------------------------------
! starting with (N-2, k-1), add desired up/dn
! or dn/up to make (N, k) 
	subroutine CreatMap1(n,k,l1,l2,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb,isk
	implicit none
	integer, intent(in) :: n,k,l1,l2
	integer, intent(in) :: ntot
	integer, dimension(2,ntot), intent(out):: map
	! local
	integer(kind=isk), dimension(k) :: set2
	integer:: n2
	integer:: ibi=1, ib ! ib1=1: dN=-2 for starting basis
	integer(kind=isk), dimension(k-1):: set
	integer:: i,k1

		ib= mapb%map(ibi);
		k1 = k-1;
		do i=1,ntot
			set = basis(ib)%sec(k1)%sets(i,:)
			! add an up site at n1 position and a down at n2 position;
			! shift n1+ and n2+ labels by 1 or 2 as they should be.
			!	shift by +1 twice for second case means +2 total as desired
			call Shift(set,k1,l1) ! shift all elem by +1 from l1 onwards
			call Shift(set,k1,l2) ! shift all elem by +1 from l2 onwards
			! Adding up/down  at l1/l2
			call SortedInsert(set,k1,l1,set2)	
			map(1,i) = LexicoIndex(set2,n,k)
			! Adding down/up  at l1/l2
			call SortedInsert(set,k1,l2,set2)
			map(2,i) = LexicoIndex(set2,n,k)
		end do
		return
	end subroutine CreatMap1

!------------------------------------------
!	Map for D Phi creation at site l1,l2 for channel 3
!------------------------------------------
! starting with (N-2, k-2), add up,up at l1,l2 to make (N, k) 
	subroutine CreatMap3(n,k,l1,l2,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb,isk
	implicit none
	integer, intent(in) :: n,k,l1,l2
	integer, intent(in) :: ntot
	integer, dimension(ntot), intent(out):: map
	! local
	integer:: k1,k2
	integer:: ibi=1,ib ! ib1=1: dN=-2 for starting basis
	integer(kind=isk), dimension(k-2):: set
	integer(kind=isk), dimension(k-1) :: set2
	integer(kind=isk), dimension(k):: set3	
	integer:: i

	if (k==2) then ! only a single final state, k'=k-2=0;
		map = (/ 1 /)
		return
	endif

	ib = mapb%map(ibi);

	k1 = k-1; k2=k-2;
	do i=1,ntot
		!write(*,*) "i, shape(sets(i,:))",i,shape(basis(ib)%sec(k2)%sets)
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
		map(i) = LexicoIndex(set3,n,k)
	end do
	return
	end subroutine CreatMap3
!------------------------------------------
!	Map for D Phi creation at site l1,l2 for channel 4
!------------------------------------------
! starting with (N-2, k), add dn,dn at l1,l2 to make (N, k) 
	subroutine CreatMap4(n,k,l1,l2,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb,isk
	implicit none
	integer, intent(in) :: n,k,l1,l2
	integer, intent(in) :: ntot
	integer, dimension(ntot), intent(out):: map
	! local
	integer:: n2
	integer:: ibi=1, ib ! ib1=1: dN=-2 for starting basis
	integer(kind=isk), dimension(k):: set	
	integer:: i,k1

	ib = mapb%map(ibi);
	
	!n2 = n+2;
	do i=1,ntot
		set = basis(ib)%sec(k)%sets(i,:)
		! add an up site at n1 position and a down at n2 position;
		! shift n1+ and n2+ labels by 1 or 2 as they should be.
		!	shift by +1 twice for second case means +2 total as desired
		call Shift(set,k,l1) ! shift all elem by +1 from l1 onwards
		call Shift(set,k,l2) ! shift all elem by +1 from l2 onwards
		map(i) = LexicoIndex(set,n,k)
	end do
	return
	end subroutine CreatMap4
!------------------------------------------



































!------------------------------------------
!	D Phi Annihilation for channel 1,2
! up,dn become D,Phi
!------------------------------------------
	subroutine DPhiCreat1(is)
	use modmain, only: basis,na,nx,Asites,ways,mapb
	
	implicit none
	integer, intent(in):: is
	!	local
	integer:: ih=7, ib1i=3, ib2i=1,ib1,ib2 ! see dnalist5 in modmain
	integer::k,n,m,m1,m2,m3
	integer :: ntot, ind, ntot1,nnz,n1,n2,n3,i
	integer, allocatable, dimension(:,:) :: map,row
	integer, allocatable, dimension(:):: pntr1,pntr2,col
	integer :: l1,l2,l
	logical :: order

	ib1= mapb%map(ib1i);
	ib2= mapb%map(ib2i);

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n-2,m-1);
	m3 = min(m1,n-1); ! =?min(m,n-1);at least 1 dn required


	l = ways(ih)%active(is); ! localtion of active site in lattice
	! order in list of active sites
	call ActiveOrder(ASites,na,l,l1,l2,order)

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2(:) = basis(ib2)%pntr(1:m2+2) ! final
	ntot1 = pntr1(m1+2); ! total number of basis states

	! nnz total number of non zero elements
	nnz = pntr2(m2+2) !- 1; ! up to pntr2(m3+2), -1 to exclude k=0 
	! allocate transition matrix: coo format
	!hop(ih)%ht(1,is)%nnz = nnz;
	!hop(ih)%ht(2,is)%nnz = nnz;

	allocate(row(nnz,2)) 
	allocate(col(nnz)) 

	col = 0

	!	calc the matrix
	!	init k=1 ===> final k'=1-1=0
	ind = 1; ntot=1
	row(1,1) = 1 + l1 ! l1 up in k=1 sec; pntr2(2) ==> 1
	row(1,2) = 1 + l2 ! l2 up in k=1 sec
	! final basis
	col(1) = 1;! k=0 sec of final

	!	init k>1 ====> final k'=k-1
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
		col(ind:ind+ntot-1) =
     .		(/ (i,i=pntr2(k)+1,pntr2(k+1),1 ) /); ! k-1 sec of final
		!write(*,*)"pntr2(k)+1,pntr2(k+1) = ",pntr2(k)+1,pntr2(k+1)
		!write(*,*)"ind,indf = ",ind,ind+ntot

		ind = ind+ntot;
		
		deallocate(map);
	end do

	! final basis
	!col=(/ (i,i=1,nnz ) /);! same as in k loop
	!write(*,*)"nnz = ",nnz
	!write(*,*)"***********************************************"
	!write(*,*)"pntr1",pntr1
	!write(*,*)"pntr2",pntr2
	!write(*,*) "row",row
	!write(*,*) "col",col

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	if (order) then
		!write(*,*)"creat: order=T"
		!write(*,*) "nnz,n3 = ",nnz,n3
		!write(*,*) row(:,1)
		call CalAmp0(ih,one,is,row(:,1),nnz,n3,col) ! ic=1
		call CalAmp0(ih,two,is,row(:,2),nnz,n3,col) ! ic=2
	else
		!write(*,*)"creat: order=F"
		!write(*,*) "nnz,n3 = ",nnz,n3
		call CalAmp0(ih,one,is,row(:,2),nnz,n3,col) ! ic=1
		call CalAmp0(ih,two,is,row(:,1),nnz,n3,col) ! ic=2
	endif
	!---------------------------------

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiCreat1
!------------------------------------------










!------------------------------------------
!	D Phi Annihilation for channel 3
! up,up become D,Phi
!------------------------------------------
	subroutine DPhiCreat3(is)
	use modmain, only: basis,na,nx,Asites,ways,mapb
	implicit none
	integer, intent(in):: is
	!	local
	integer:: ih=7, ib1i=3, ib2i=1,ib1,ib2 ! see dnalist5 in modmain
	integer::n,k,m,m1,m2
	integer :: ntot, ind, nnz,n1,n2,n3,i1,i2,i
	integer, allocatable, dimension(:) :: map,row,col
	integer, allocatable :: pntr1(:), pntr2(:)
	integer:: l1,l2,l
	logical :: order

	ib1= mapb%map(ib1i);
	ib2= mapb%map(ib2i);

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n-2,m-2);

	l = ways(ih)%active(is); ! localtion of active site in lattice
	! order in list of active sites
	call ActiveOrder(ASites,na,l,l1,l2,order)

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	! ib: itype ===> which of 5 N case?
	pntr1 = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2 = basis(ib2)%pntr(1:m2+2) ! final
	!ntot1 = pntr1(m1+2); ! total number of basis states

	! nnz total number fo non zero elements
	nnz = pntr2(m2+2);
	! allocate transition matrix: coo format
	!hop(ih)%ht(3,is)%nnz = nnz;

	allocate(row(nnz)) 
	allocate(col(nnz)) 
	!write(*,*) "pntr1 = ",pntr1
	!write(*,*) "pntr2 = ",pntr2
	!	2 >= k <= m1; at least two up
	ind = 1;
	do k=2,m1,1
		!write(*,*) " ------------  k = ",k
		!write(*,*) "shape pntr2 = ",shape(pntr2)
		!write(*,*) "pntr2 1= ",pntr2
		ntot = pntr2(k)-pntr2(k-1)
		!write(*,*) "k, ntot = ",k,ntot
		allocate(map(ntot));
		!write(*,*) "pntr2 2= ",pntr2
		!write(*,*) "k, i1,i2 = ",k,i1,i2
		!map = 0
		call CreatMap3(n,k,l1,l2,map,ntot)
		!write(*,*) "pntr2 3= ",pntr2		
		!write(*,*) "pntr1 ?= ",pntr1
		!write(*,*) "map= ",map
		!if (k==3)stop	
		! initial basis
		map = pntr1(k+1) + map
		row(ind:ind+ntot-1) = map
		! final basis
		col(ind:ind+ntot-1) =
     .		(/ (i,i=pntr2(k-1)+1,pntr2(k),1 ) /) ! k-2 sec of final
		ind = ind+ntot
		deallocate(map)
	end do

	!write(*,*) "row = ",row

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space	
	call CalAmp0(ih,three,is,row,nnz,n3,col) ! ic=3
	!---------------------------------



	deallocate(pntr1,pntr2,row,col)
	return	
	end subroutine DPhiCreat3
!------------------------------------------



!------------------------------------------
!	D Phi Annihilation for channel 4
! dn,dn become D,Phi
!------------------------------------------
	subroutine DPhiCreat4(is)
	use modmain, only: basis,na,nx,Asites,ways,mapb
	implicit none
	integer, intent(in):: is
	!	local
	integer:: ih=7, ib1i=3, ib2i=1,ib1,ib2! see dnalist5 in modmain
	integer::k,n,m,m1,m2
	integer :: ntot, ind, nnz,n1,n2,n3,i
	integer, allocatable, dimension(:) :: map,row,col
	integer, allocatable, dimension(:):: pntr1,pntr2
	integer :: l1,l2,l
	logical :: order

	ib1= mapb%map(ib1i);
	ib2= mapb%map(ib2i);

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n-2,m); ! m3=min(n-2,m)=m2 at least two up required

	l = ways(ih)%active(is); ! localtion of active site in lattice
	! order in list of active sites
	call ActiveOrder(ASites,na,l,l1,l2,order)

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2(:) = basis(ib2)%pntr(1:m2+2) ! final
	!ntot1 = pntr1(m1+2); ! total number of basis states

	! nnz total number fo non zero elements
	nnz = pntr2(m2+2);
	! allocate transition matrix: coo format

	allocate(row(nnz)) 
	allocate(col(nnz)) 

	! final | initial k=0 case: by hand for efficiency !
	! no actual need to seperate from the k loop below
	row(1) = 1 
	col(1) = 1 

	!	1 >= k <= m3=m2 cases; at least two dn req
	ind = 2;
	do k=1,m2,1
		ntot = pntr2(k+2)-pntr2(k+1) ! k sec of final
		allocate(map(ntot));
		call CreatMap4(n,k,l1,l2,map,ntot)
		! initial basis
		map = pntr1(k+1) + map
		row(ind:ind+ntot-1) = map
		! final basis
		col(ind:ind+ntot-1) = 
     .			(/ (i,i=pntr2(k+1)+1,pntr2(k+2),1 ) /) ! k sec of final
		ind = ind+ntot
		deallocate(map)
	end do

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3=pntr1(m1+2) ! dim of initial hilbert space
	call CalAmp0(ih,four,is,row,nnz,n3,col) ! ic=4
	!---------------------------------

	deallocate(pntr1,pntr2)
	return
	end subroutine DPhiCreat4
!------------------------------------------

!---------------------------------------
	subroutine ActiveOrder(ASites,na,l,l1,l2,order)
	use modmain, only: sys, periodic,isk
	implicit none
	integer, intent(in) :: na,l
	integer, dimension(na), intent(in):: Asites
	integer, intent(out) :: l1,l2
	logical, intent(out) :: order
	! local
	integer:: x1,x2,lp1

	if(l < sys%nsites) then
		lp1 = l + 1;
	elseif(l==sys%nsites .and. periodic)then
		lp1 = 1;
	else
		write(*,*)"Error(ActiveOrder): ?!"
	endif
	
	! order in the list of active sites
	x1 = GetPosition(ASites,na,l);
	x2 = GetPosition(ASites,na,lp1);
	if (x1<x2) then
		l1=x1; l2=x2; order=.true.;
	else
		l2=x1; l1=x2; order=.false.;
	endif

	return
	end subroutine ActiveOrder
!----------------------------------

	end module Creation
