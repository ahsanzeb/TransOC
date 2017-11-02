
	module Contacts
	use amplitudes
	use lists, only: SortedInsert,GetPosition
	use basisstates, only: LexicoIndex,Shift

	implicit none

	public:: CDAnnihil, CDCreat
	private::CDAnMap, CDCrMap1,CDCrMap2

	contains
!------------------------------------------
!	Map for D Annihilation at a contact
	subroutine CDAnMap(ibl,n,k,map,la)
	! D site appended to active list.
	!? out: index for D,Phi ==> up,dn
	use modmain, only: basis,isk
	use basisstates, only: LexicoIndex
	!use lists, only: MemberQ
	implicit none
	integer, intent(in) :: ibl,n,k
	integer, intent(in) :: la
	integer, dimension(la,2), intent(out):: map
	! local
	integer(kind=isk), dimension(k+1) :: set2
	integer:: n1
	integer:: i,k1

	n1 = n+1;
	k1 = k+1;
	do i=1,la
		set2(1:k) = basis(ibl)%sec(k)%sets(i,:)
		set2(k1) = n1; ! an up at n+1
		map(i,1) = LexicoIndex(set2,n1,k1); !
		map(i,2) = LexicoIndex(set2(1:k),n1,k)	;! dn at n+1:
	end do
	return
	end subroutine CDAnMap
!-----------------------------------




!-----------------------------------
!	D Phi Annihilation at a contact
!------------------------------------------
	subroutine CDAnnihil()
	use modmain, only: basis,na,nx,mapb
	implicit none
	!	local
	integer:: ib1i=3, ib2i=4 ! see dnalist5 in modmain
	integer:: ib1,ib2
	integer::k,n,m,m1,m2,i
	integer :: ntot, ind, n1,n3
	integer, allocatable, dimension(:,:) :: map,col
	integer, allocatable, dimension(:):: pntr1,pntr2

	ib1= mapb%map(ib1i);
	ib2= mapb%map(ib2i);
	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations
	if (n==0) then
		! allocate transition matrix: diagonal format
		allocate(col(1,2))
		col(1,1) = 2 ! Jh
		col(1,2) = 1 ! Jl
		n3=1;! dim of initial hilbert space
	else
		m1 = min(m,n); ! max no of up spins possible
		n1 = n+1;
		m2 = min(m+1,n1);
		allocate(pntr1(m1+2))
		allocate(pntr2(m2+2))
		! ib: itype ===> which of 5 N case?
		pntr1(:) = basis(ib1)%pntr(1:m1+2) ! only the relevant part
		pntr2(:) = basis(ib2)%pntr(1:m2+2) 
		n3 = pntr1(m1+2); ! total number of basis states
		! allocate transition matrix: diagonal format
		allocate(col(n3,2)) 
	
		!	calc the matrix
		! k==0 
		col(1,1) = 1 + n1; ! Jh: final k=1
		col(1,2) = 1; ! Jl: final k=0

		ind = 2;
		do k=1,m1
			ntot = pntr1(k+2)-pntr1(k+1)
			allocate(map(ntot,2));
			call CDAnMap(ib1,n,k,map,ntot)
			! Jh case: shift k by 1; final k+1 up
			col(ind:ind+ntot-1,1) = pntr2(k+2) + map(:,1)
			! Jl case: final k up
			col(ind:ind+ntot-1,2) = pntr2(k+1) + map(:,2) 
			ind = ind+ntot
			deallocate(map)
		end do
		deallocate(pntr1,pntr2)
	endif !n=0

	! map in col is for D annihilation but 
	!	will also be used for Phi annihilation.
	! ih=9 ==> ia = 7; amp for ih=9,11,13,15: final m excitations
	! ih=10 ==> ia = 8; amp for ih=10,12,14,16: final m+1 excitations
	! calculate transition amplitudes
	call CalAmp(9,1,1,col(:,2),n3,n3,"multiplydc");
	call CalAmp(10,1,1,col(:,1),n3,n3,"multiplydc"); 

	return
	end subroutine CDAnnihil
!------------------------------------------




















!------------------------------------------
!	Maps for D creation at a contact
!------------------------------------------
! starting with (N-1, k-1), add an up ==> (N, k) 
	subroutine CDCrMap1(n,k,l1,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb,isk
	implicit none
	integer, intent(in) :: n,k,l1
	integer, intent(in) :: ntot
	integer, dimension(ntot), intent(out):: map
	! local
	integer(kind=isk), dimension(k) :: set2
	integer:: ibi=2, ib
	integer(kind=isk), dimension(k-1):: set
	integer:: i,k1

		ib= mapb%map(ibi);
		k1 = k-1;
		do i=1,ntot
			set = basis(ib)%sec(k1)%sets(i,:)
			! add an up site at n1 position
			! shift l1+ labels by 1;
			call Shift(set,k1,l1) ! shift all elem by +1 from l1 onwards
			! Adding up at l1
			call SortedInsert(set,k1,l1,set2)	
			map(i) = LexicoIndex(set2,n,k)
		end do
		return
	end subroutine CDCrMap1
!-----------------------------------------------
! starting with (N-1, k), add a down ==> (N, k) 
	subroutine CDCrMap2(n,k,l1,map,ntot)
	! only k > 0
	use modmain, only: basis,mapb,isk
	implicit none
	integer, intent(in) :: n,k,l1
	integer, intent(in) :: ntot
	integer, dimension(ntot), intent(out):: map
	! local
	integer:: ibi=2, ib
	integer(kind=isk), dimension(k):: set
	integer:: i

		ib= mapb%map(ibi);
		do i=1,ntot
			set = basis(ib)%sec(k)%sets(i,:)
			! add an up site at n1 position
			! shift l1+ labels by 1;
			call Shift(set,k,l1) ! shift all elem by +1 from l1 onwards
			! Adding a down l1: same set
			map(i) = LexicoIndex(set,n,k)
		end do
		return
	end subroutine CDCrMap2
!-----------------------------------------------














!------------------------------------------
!	D creation at a contact
!------------------------------------------
	subroutine CDCreat(cont)
	use modmain, only: basis,na,nx,Asites,mapb,nsites
	
	implicit none
	character(len=1), intent(in):: cont
	!	local
	integer:: ib1i=3, ib2i=2,ib1,ib2 ! see dnalist5 in modmain
	integer::k,n,m,m1,m2,m3
	integer :: ntot, ind, ntot1,nnz,n1,n2,n3,i
	integer, allocatable, dimension(:) :: map,row,col
	integer, allocatable, dimension(:):: pntr1,pntr2,pntr3
	integer :: is=1, ih,l1

	ib1= mapb%map(ib1i);
	ib2= mapb%map(ib2i);

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = max(0,min(n-1,m-1));! set m2=0 if m=0, needed for Jl case
	m3 = min(m1,n-1); ! =?min(m,n-1);at least 1 dn required

	!write(*,*)"n,m,m2,m3 = ",n,m,m2,m3

	!write(*,*)"cont = ",cont
	!write(*,*)"ASites = ",	ASites

	! position of active contact site in Asites list.
	if(cont=='l') then
		l1 = GetPosition(ASites,n,1);
	elseif(cont=='r') then
		l1 = GetPosition(ASites,n,nsites);
	endif

	allocate(pntr1(m1+2))
	allocate(pntr2(m2+2)) 
	allocate(pntr3(m3+2)) 

	! ib: itype ===> which of 5 N case?
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial
	pntr2(:) = basis(ib2)%pntr(1:m2+2) ! final
	pntr3(:) = basis(ib2)%pntr(1:m3+2) ! final

	n3 = pntr1(m1+2); ! total number of basis states

	!=================== Jh ====================
	!						k ===> k-1
	if (m > 0) then
	! nnz total number of non zero elements
	nnz = pntr2(m2+2); 
	allocate(row(nnz)) 
	allocate(col(nnz)) 
	col = 0

	!	calc the matrix
	!	init k=1 ===> final k'=1-1=0
	! ind = 1;
	row(1) = 1 + l1 ! l1 up in k=1 sec; pntr2(2) ==> 1
	! final basis
	col(1) = 1;! k=0 sec of final

	!	init k>1 ====> final k'=k-1
	ind = 2;
	do k=2,m1
		ntot = pntr2(k+1)-pntr2(k) ! k-1 sec of final
		allocate(map(ntot));
		call CDCrMap1(n,k,l1,map,ntot)
		! initial basis
		row(ind:ind+ntot-1) = pntr1(k+1) + map
		! final basis
		col(ind:ind+ntot-1) =
     .		(/ (i,i=pntr2(k)+1,pntr2(k+1),1 ) /); ! k-1 sec of final
		ind = ind+ntot;		
		deallocate(map);
	end do

	! calculate transition amplitudes
	! ih=19,20; ia = 10; Right conatact
	! ih=23,24; ia = 12; Left contact
	! ih ==> ia in Amplitudes
	if(cont=='l') then ! l=1, Left contact
		ih = 23;
	elseif(cont=='r') then ! l=nsites, Right conatact
		ih = 19;
	endif
	call CalAmp0(ih,1,1,row,nnz,n3,col)
	
	deallocate(row,col)
	endif ! m > 0
	!=================== Jl ====================
	!						k ===> k
	! nnz total number of non zero elements
	nnz = pntr3(m3+2)
	allocate(row(nnz)) 
	allocate(col(nnz)) 
	col = 0

	!	calc the matrix
	!	init k=0
	row(1) = 1; col(1) = 1;

	!write(*,*)"nnz = ",nnz


	!	init k>0
	ind = 2;
	do k=1,m3
		ntot = pntr3(k+2)-pntr3(k+1) ! k sec of final
		allocate(map(ntot));
		call CDCrMap2(n,k,l1,map,ntot)
		! initial basis
		row(ind:ind+ntot-1) = pntr1(k+1) + map
		! final basis
		col(ind:ind+ntot-1) =
     .		(/ (i,i=pntr3(k+1)+1,pntr3(k+2) ) /); ! k sec of final
		ind = ind+ntot;		
		deallocate(map);
	end do

	! calculate transition amplitudes
	! ih=17,18; ia = 9; Right conatact
	! ih=21,22; ia = 11; Left contact
	! ih ==> ia in Amplitudes
	if(cont=='l') then ! l=1, Left contact
		ih = 21;
	elseif(cont=='r') then ! l=nsites, Right conatact
		ih = 17;
	endif
	call CalAmp0(ih,1,1,row,nnz,n3,col)

	deallocate(pntr1,pntr2,row,col)
	
	return
	end subroutine CDCreat
!------------------------------------------









	end module Contacts
