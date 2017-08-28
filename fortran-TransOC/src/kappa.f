


	subroutine LossKappa()
	use modmain, only: basis,hop,na,nx
	implicit none
	!	local
	integer(kind=1):: ih=25, is=1, ib1=3, ib2=3 ! see dnalist5 in modmain
	integer(kind=1)::n,m,m1,m2,i
	integer :: ntot1, ntot2

	!------------------------------------------	
	! N,m values of itype:
	n = na; ! no of active sites
	m = nx; ! no of excitations

	m1 = min(n,m);
	m2 = min(n,m-1);

	! ib: itype ===> which of 5 N case?
	ntot1 = basis(ib1)%pntr(m1+2) ! initial
	ntot2 = basis(ib2)%pntr(m2+2) ! final

	! allocate transition matrix: diag format
	hop(ih)%ht(1,is)%nnz = ntot2;
	if (allocated(hop(ih)%ht(1,is)%col))
     .			deallocate(hop(ih)%ht(1,is)%col)
	allocate(hop(ih)%ht(1,is)%col(ntot2)) 

	! spin combinations map to themselves, only a photon is lost
	!	basis with up to m2=min(n,m-1) up spins contribute non-zero
	hop(ih)%ht(1,is)%col(:) = (/ (i,i=1,ntot2,1) /)

	return
	end subroutine LossKappa






