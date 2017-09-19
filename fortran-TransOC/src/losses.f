
	module losses
	use amplitudes
	
	implicit none

	public:: LossGamma,LossKappa
	private:: GammaMap

	contains
	

	subroutine GammaMap(ib,n,k,l,ntot1,map,ntot)
	use modmain, only: basis, mapb,isk
	use basisstates, only: LexicoIndex
	use lists, only: Drop, MemberQ
	implicit none
	integer, intent(in) :: ib,n,k,l
	integer, intent(in) :: ntot, ntot1
	integer, dimension(2,ntot), intent(out):: map
	! local
	integer:: k1
	!integer(kind=1):: ibi=3, ib
	integer(kind=isk), dimension(k):: set
	integer(kind=isk), dimension(k-1) :: set2
	
	integer:: i, ind

	!ib = mapb%map(ibi);

	k1 = k-1;
	ind = 1
	do i=1,ntot1
		set = basis(ib)%sec(k)%sets(i,:)
		if (MemberQ(set,k,l)) then
			map(1,ind) = i
			call Drop(set,k,l,set2)
			map(2,ind) = LexicoIndex(set2,n,k1)
			ind = ind + 1;
		endif
	end do
	return
	end subroutine GammaMap





	subroutine LossGamma(is)
	use modmain, only: basis,na,nx,mapb
	implicit none
	integer, intent(in):: is 
	!	local
	integer:: ih=26, ib1i=3 !, ib2i=3 ! see dnalist5 in modmain
	integer::k,n,m,m1,n1,n2,m2,i, ib1 !,ib2
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

		!write(*,*)"losses: map(2,:) = ",map(2,:) 

		inda = inda+la;
		deallocate(map)
	end do

	!write(*,*)"----===---===--=="
	!-------------------------------------------------------
	! calculate transition amplitudes
	n3 = pntr1(m1+2);
	!write(*,*) "gamma: n3,lat",n3,lat
	!write(*,*) "gamma: col(1:5):",col(1:5)
	!-------------------------------------------------------
	call CalAmp0(ih,1,is,row,lat,n3,col) ! ic=1; "multiply--"
	!---------------------------------

	!write(*,*) "gamma: after callamp"


	deallocate(pntr1,las)
	return
	end subroutine LossGamma





	subroutine LossKappa()
	use modmain, only: basis,na,nx,mapb
	implicit none
	!	local
	integer:: ih=25, is=1, ib1=3, ib2=3 ! see dnalist5 in modmain
	integer::n,m,m1,m2
	integer :: ntot1, ntot2,ibl1,ibl2,i
	integer, allocatable, dimension(:):: col
	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n,m-1);
	ibl1= mapb%map(ib1);
	ibl2=mapb%map(ib2);


	! ib: itype ===> which of 5 N case?
	ntot1 = basis(ibl1)%pntr(m1+2) ! initial
	ntot2 = basis(ibl2)%pntr(m2+2) ! final

	! allocate transition matrix: diag format
	!hop(ih)%ht(1,is)%nnz = ntot2;
	allocate(col(ntot2))

	! spin combinations map to themselves, only a photon is lost
	!	basis with up to m2=min(n,m-1) up spins contribute non-zero
	col(:) = (/ (i,i=1,ntot2,1) /)
	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------

	!write(*,*) "kappa: calling amp"

	!write(*,*) "kappa: ntot2",ntot2
	!write(*,*) "kappa: rowc(1:5)",col(1:min(5,ntot2))

	!write(*,*) "len of multiplyd- =",len("multiplyd-")
	call CalAmp(ih,1,is,col,ntot2,ntot1,'multiplyd') ! ic=1
	!---------------------------------

	!write(*,*) "kappa: done amp"

	return
	end subroutine LossKappa



	end module losses







