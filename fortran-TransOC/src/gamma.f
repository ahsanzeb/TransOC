
	module gamma
	implicit none

	public:: LossGamma
	private:: GammaMap

	contains
	

	subroutine GammaMap(ib,n,k,l,ntot1,map,ntot)
	use modmain, only: basis, mapb
	use basisstates, only: LexicoIndex
	use lists, only: Drop, MemberQ
	implicit none
	integer(kind=1), intent(in) :: ib,n,k,l
	integer(kind=4), intent(in) :: ntot, ntot1
	integer(kind=4), dimension(2,ntot), intent(out):: map
	! local
	integer(kind=1):: k1
	!integer(kind=1):: ibi=3, ib
	integer(kind=1), dimension(k):: set
	integer(kind=1), dimension(k-1) :: set2
	
	integer:: i, ind

	!ib = mapb%map(ibi);

	k1 = k-1;
	ind = 1
	do i=1,ntot1
		set = basis(ib)%sec(k)%sets(i,:)
		if (MemberQ(set,k,l)) then
			map(1,ind) = i
			call Drop(set,k,l,set2)
			map(1,ind) = LexicoIndex(set2,n,k1)
			ind = ind + 1;
		endif
	end do
	return
	end subroutine GammaMap





	subroutine LossGamma(is)
	use modmain, only: basis,na,nx,mapb
	implicit none
	integer(kind=1), intent(in):: is 
	!	local
	integer(kind=1):: ih=26, ib1i=3 !, ib2i=3 ! see dnalist5 in modmain
	integer(kind=1)::k,n,m,m1,n1,n2,m2,i, ib1 !,ib2
	integer :: ntot, ind, ntot1,nnz,inda,la,lat,n3
	integer, allocatable, dimension(:,:) :: map
	integer, allocatable, dimension(:)::pntr1,las,row,col

	ib1 = mapb%map(ib1i);
	!ib2 = mapb%map(ib2i);
	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n,m-1);

	allocate(pntr1(m1+2))
	pntr1(:) = basis(ib1)%pntr(1:m1+2) ! initial

	! dimensions of maps for diff k
	allocate(las(m1))
	las(:) = 0;
	do k=1,m1
		ntot = pntr1(k+2)-pntr1(k+1)
		las(k) = k*ntot/n ! number of basis with a given spin up
	end do
	! dimensions of full transition matrices
	lat = sum(las)
	! allocate transition matrix: coo format
	allocate(row(lat)) 
	allocate(col(lat))

	!hop(ih)%ht(1,is)%nnz = lat ! nnz 
	
	!	calc the matrix
	! k>0 only, need a spin up to decay it!
	inda = 1;
	do k=1,m1,1
		la = las(k);
		allocate(map(2,la));
		ntot1 = pntr1(k+2) - pntr1(k+1); 
		!	calc maps
		call GammaMap(ib1,n,k,is,ntot1,map,la)
		!	assign values to transition matrices
		row(inda:inda+la-1) = pntr1(k+1) + map(1,:)
		col(inda:inda+la-1) = pntr1(k) + map(2,:)
		inda = inda+la;
		deallocate(map)
	end do

	!-------------------------------------------------------
	! calculate transition amplitudes
	n3 = pntr1(m1+2);
	!-------------------------------------------------------
	call CalAmp(ih,1,is,row,lat,n3,"multiply",col) ! ic=1
	!---------------------------------



	deallocate(pntr1,las)
	return
	end subroutine LossGamma



	end module gamma







