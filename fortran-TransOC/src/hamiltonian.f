
	module hamiltonian
	implicit none
	public :: mkHamilt, MakeHgMulti
	public ::  DegenSectors
	private:: HgMap

	contains

!------------------------------------------
	subroutine mkHamilt()
	use modmain, only: mapt,Hg
	implicit none
	! local
	integer(kind=1):: nt,ib
	
	do ib=1,5
		nt = mapt%ntb(ib)
		if (nt > 0) then ! something to calculate or not?
			call MakeHgMulti(mapt%grouptb(ib,1:nt),nt)
		endif
	end do


	!write(*,*)"hamilt: Hg(:)%xst ============="
	!write(*,*)	Hg(:)%xst
	!write(*,*)	Hg(1)%rowpntr


	
	return
	end	subroutine mkHamilt
!------------------------------------------

!------------------------------------------
	subroutine HgMap(ibl,n,k,ntot,map)
	use modmain, only: basis
	use lists, only: Complement,SortedInsert !SortedInsert=Append+Sort
	use basisstates, only: LexicoIndex
	implicit none
	integer(kind=1), intent(in) :: ibl,n,k
	integer(kind=4), intent(in) :: ntot
	integer(kind=4), dimension(ntot,n-k), intent(out):: map
	! local
	integer:: nk,i,j
	integer(kind=1), dimension(k):: set1
	integer(kind=1), dimension(k+1):: set2
	integer(kind=1), dimension(n-k):: flipsites
	integer(kind=1) k1,flip
	integer(kind=1), allocatable :: sites(:)
	
	if (k==0) then
		map(1,:) = (/ (i,i=1,n,1)/)
	elseif(k==n) then
		! all map to the single state with all up
		map(:,1) = 1;
	else
		! initialise sites array
		if(allocated(sites))deallocate(sites)
		allocate(sites(n))
		sites(:) = (/ (i,i=1,n,1)/)
		
		nk = n-k;
		!write(*,*) "MapHg: k, ntot = ",k, ntot
		do i=1,ntot
		! sec(k) are defined for k>0,
		! so index match with k
			!write(*,*) " k,i = ",k,i
			set1=basis(ibl)%sec(k)%sets(i,:)
			!write(*,*) "sites ",sites
			!write(*,*) "k,i: ",k,i,", set1=",set1

			call Complement(sites,n,set1,k,flipsites)
			do j=1,nk
				flip = flipsites(j);
				!call Append(set1,k,flip,set2)
				k1 = k+1;
				!set2(1:k) = set1; set2(k1) = flip;
				!call Sort(set2,k1);
				call SortedInsert(set1,k1,flip,set2)
				map(i,j) = LexicoIndex(set2,n,k1)
			end do
		end do
	endif
	return
	end subroutine HgMap
!------------------------------------------
	subroutine MakeHgMulti(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	! saves the Hamiltonians in CSR format
	use modmain, only: basis,na,nx,ibs,dna,dnx,
     .								g,dw,detuning,Hg,mapb,mapt
	implicit none
	integer(kind=1), intent(in) :: nt
	integer(kind=1), dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer(kind=4):: indf,i,j,ind,ntot,itype,indi
	integer(kind=1):: ib,k,m,m1,n,ibl,i2,it2
	double precision:: val,kdw
	integer(kind=4), allocatable, dimension(:):: pntr
	integer(kind=4), allocatable, dimension(:):: ind1
	integer(kind=4), allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,nnzg,nnzd,itl
	integer :: ptr,ptrf,ptr0
	!character*1000:: fmt

	!write(*,*)"itlist ",itlist


	! JUST A CHECK: if all types have the same n
	n = dna(itlist(1))
	do i=2,nt,1
		itype = itlist(i)
		if (n .ne. dna(itlist(i))) then
			write(*,*) "makeHg: different N in itlist, aborting."
			stop
		endif
	end do

	! N,m values
	n = na + dna(itlist(1)); ! no of active sites
	do i=1,nt,1
		m	=	 nx + dnx(itlist(i)); ! no of excitations
		ms(i) = m;
		m1s(i)=min(m,n); ! max no of up spins possible
	end do
	m1 = maxval(m1s) ! max m1

	! ib: itype ===> which of 5 N case?
	ib = ibs(itlist(1)); ! same for all cases in itlist
	ibl = mapb%map(ib); ! get the location of data for this ib

	!write(*,*)"n, ib,ibl",n,ib,ibl



	! pointers for start index
	if(allocated(pntr))deallocate(pntr)
	allocate(pntr(m1+2));
	pntr(:) = basis(ibl)%pntr(1:m1+2) ! is this already calculated?

	nnzg=0; nnzd=0;
	do k=0,m1-1
		ntot = pntr(k+2) - pntr(k+1)
		do i=1,nt
			if(k < m1s(i)) nnzg(i) = nnzg(i) + ntot*(n-k);
		end do
	end do

	!	number of diagonal elements
	if (detuning) then
		do i=1,nt
			nnzd(i) = pntr(m1s(i)+2); ! total no of basis states
		end do
	endif
	! total number of non zero elements
	nnz = nnzg + nnzd;

	! hilbert space dimension
	itl = 0;		
	do i=1,nt
		itl(i) = mapt%map(itlist(i));! map for the location of itype
		Hg(itl(i))%xst = .true.
		Hg(itl(i))%n = n;
		Hg(itl(i))%m = ms(i);
		Hg(itl(i))%ntot = pntr(m1s(i)+2);
		Hg(itl(i))%nnz = nnz(i);
		!write(*,*) "ham: it, itl = ",i, itl(i)
		!write(*,*) "ham: n, m, m1 = ",n, ms(i),m1s(i)
		!write(*,*) "ham: nnz, ntot = ",nnz(i),pntr(m1s(i)+2)
	end do
		!write(*,*) "ham: DETUNING = ", detuning

	! allocate memory to Hg(itype)
	do i=1,nt
		if(allocated(Hg(itl(i))%col)) then
			!deallocate(Hg(itl(i))%row)
			deallocate(Hg(itl(i))%col)
			deallocate(Hg(itl(i))%dat)
			deallocate(Hg(itl(i))%rowpntr)
		endif
		!allocate(Hg(itl(i))%row(nnz(i)))
		allocate(Hg(itl(i))%col(nnz(i))); 
		allocate(Hg(itl(i))%dat(nnz(i)));
		! rowpntr size = pntr(m1s(i)+1), one less sectors 
		! k=m1 sector has no higher ksub sector to couple to.
		! +1 for last value =nnz+1
		Hg(itl(i))%srptr= pntr(m1s(i)+1) + 1; ! size of row pntr
		allocate(Hg(itl(i))%rowpntr(pntr(m1s(i)+1) +1))

		Hg(itl(i))%col=-10;
		Hg(itl(i))%dat=0.0d0;
		Hg(itl(i))%rowpntr = -10;

		! set last value of rowptr
		Hg(itl(i))%rowpntr(pntr(m1s(i)+1) +1) = nnz(i) + 1;
		
		!write(*,*)"itl(i),shape(Hg(itl(i))%col)"
		!write(*,*)itl(i),shape(Hg(itl(i))%col)
		!write(*,*)itl(i),shape(Hg(itl(i))%rowpntr)
		!if(detuning) then
		!	allocate(Hg(itl(i))%rowpntr(2*pntr(m1s(i)+2))) ! 2*hilbert space dim???
		!else
		!	allocate(Hg(itl(i))%rowpntr(pntr(m1s(i)+2))) ! hilbert space dim
		!endif
	end do
			
	ind = 1; ptr=1; ptr0=0;
	! calc maps for diff k and combine to get full matrix
	do k=0,m1 ! =m1 for detuning case, to include diag elem
		ntot = pntr(k+2) - pntr(k+1);
		if (k<m1) then
			if(allocated(ind1))deallocate(ind1)
			allocate(ind1(ntot));
			ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);
			! calc map from k to k+1 up spin sector
			if(allocated(map))deallocate(map)
			allocate(map(ntot,n-k));
			call HgMap(ibl,n,k,ntot,map)
		
			! assign values to final matrix			
			indf=ind+ntot*(n-k)-1;
			indi = ind; 
			ptrf = ptr + ntot -1
			do i2=1,nt
				ind = indi; ! save the index for i2 loop
				if(k < m1s(i2)) then
					val = dsqrt((ms(i2)-k)*1.d0)*g;
					it2 = itl(i2);
					Hg(it2)%dat(ind:indf) = val;
					! transpose map to get col first for reshape
					Hg(it2)%col(ind:indf) = 
     .     reshape(transpose(map), (/ntot*(n-k)/))
					Hg(it2)%rowpntr(ptr:ptrf) =
     .    (/ (ptr0 + (i-1)*(n-k)+1, i=1,ntot) /)
     			!if(it2==1) then
					!	write(*,*)"Hg(it2)%rowpntr(ptr:ptrf), m1, k=",m1s(it2),k
					!	write(*,*)Hg(it2)%rowpntr(ptr:ptrf)
					!endif
				endif
			end do
			ind = indf+1; ! advance the index
		endif ! k < m1
		
		! diagonal elements 
		if (detuning) then
			indf=ind+ntot-1;
			indi = ind; ! save the index for i2 loop
			do i2=1,nt
				ind = indi;
				!-----------------------------------------
				! store half of actual values to avoid double counting
				! in mat vec multiplication in module diag, subroutine matvec()
				! ===> put a switch to control this behaviour in case
				! a direct eigensolver is to be used instead of arpack
				!-----------------------------------------
				kdw = (ms(i2) - k)*dw/2.0d0;
				if(k <= m1s(i2)) then
					Hg(itl(i2))%dat(ind:indf) = kdw;
					!Hg(itl(i2))%row(ind:indf) = ind1(:);
					Hg(itl(i2))%col(ind:indf) = ind1(:);
					ind = indf + 1;
					! shift pntrs, by one for each prev row
					Hg(itl(i2))%rowpntr(ptr:ptrf) = 
     .   		Hg(itl(i2))%rowpntr(ptr:ptrf) + 
     .   			(/(i,i=0,ntot-1)/)
				endif
			end do
			ptr = ptrf+1	
			ptr0 = ptr0 + ntot*(n-k+1); ! shift for nxt k-blocks		
		else
			ptr = ptrf+1		
			ptr0 = ptr0 + ntot*(n-k); ! shift for nxt k-blocks		
		end if ! detuning

	end do ! k

	return
	end subroutine MakeHgMulti
!------------------------------------------



!---------------------------------------
	subroutine DegenSectors(es,ne,nsec,esec,ind)
	! importnat:only if es are in ascending order
	integer, intent(in) :: ne
	double precision,dimension(ne),intent(in):: es
	double precision,dimension(ne),intent(out):: esec ! 1:nsec contains esec
	integer, intent(out) :: nsec
	integer,dimension(ne),intent(out):: ind !start index of sectors; 1:nsec
	! local
	double precision:: tol = 1.0d-6 ! tolerance
	integer:: i,j

	esec(1) = es(1); j = 1; ! start with the lowest energy
	ind(1) = 1;
	do i=1,ne
		if (es(i) > esec(j) + tol ) then
			j = j + 1; 
			esec(j) = es(i);
			ind(j) = i;
		endif		
	end do

	nsec = j; ! total number of sectors
	
	return
	end subroutine DegenSectors
!---------------------------------------

	
	end module


	
