
	module hamiltonian
	implicit none
	public :: mkHamilt
	private:: MakeHgMulti, MakeHgMulti1, HgMap

	contains

!------------------------------------------
	subroutine mkHamilt()
	use modmain, only: mapt,Hg,detuning,sameg
	implicit none
	! local
	integer(kind=1):: nt,ib

		do ib=1,5
			nt = mapt%ntb(ib)
			if (nt > 0) then ! something to calculate or not?
				if ((.not. detuning) .and. sameg) then
					! use better storage format (and matvec routines for iter diag)
					write(*,*) "hamilt: better = T"
					call MakeHgMulti1(mapt%grouptb(ib,1:nt),nt)
				else
					write(*,*) "hamilt: better = F"
					! simple CSR format to store sparse
					call MakeHgMulti(mapt%grouptb(ib,1:nt),nt)
				endif
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
     .	 g,dw,detuning,Hg,eig,mapb,mapt,smalln,sameg,ndsec
	implicit none
	integer(kind=1), intent(in) :: nt
	integer(kind=1), dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer(kind=4):: indf,i,j,ind,ntot,itype,indi,irow,icol
	integer(kind=1):: ib,k,m,m1,n,ibl,i2,it2
	double precision:: val,kdw
	integer(kind=4), allocatable, dimension(:):: pntr
	integer(kind=4), allocatable, dimension(:):: ind1
	integer(kind=4), allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,nnzg,nnzd,itl
	integer :: ptr,ptrf,ptr0,nev
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
		Hg(itl(i))%m1 = m1s(i);
		Hg(itl(i))%ntot = pntr(m1s(i)+2);
		Hg(itl(i))%nnz = nnz(i);
		if (pntr(m1s(i)+2) .le. smalln) then
			Hg(itl(i))%dense = .true.
			Hg(itl(i))%nev = ntot; ! ncv irrelevant for direct
			! if true, save hamiltonian in eig(itl(i))%evec
		else
			! set this in case it is true from a prev run
			Hg(itl(i))%dense = .false.
			nev = pntr(min(ndsec,m1s(i))+2); ! ~ nsnev degen sectors
			Hg(itl(i))%nev = nev;
			Hg(itl(i))%ncv = min(2*nev, ntot);
		endif		
	end do

	! allocate memory to Hg(itype)
	do i=1,nt
		if(allocated(Hg(itl(i))%col)) then
			deallocate(Hg(itl(i))%col)
			deallocate(Hg(itl(i))%rowpntr)
			deallocate(Hg(itl(i))%dat)
		endif
		if (Hg(itl(i))%dense) then ! dense
			ntot = pntr(m1s(i)+2)
			eig(itl(i))%ntot = ntot;
			eig(itl(i))%n1 = ntot;
			if(allocated(eig(itl(i))%evec))deallocate(eig(itl(i))%evec)
			if(allocated(eig(itl(i))%eval))deallocate(eig(itl(i))%eval)
			allocate(eig(itl(i))%evec(ntot,ntot))
			eig(itl(i))%evec(:,:) = 0.0d0 ! initialise = 0
			allocate(eig(itl(i))%eval(ntot))
		else	 ! sparse
			allocate(Hg(itl(i))%dat(nnz(i)));
			!allocate(Hg(itl(i))%row(nnz(i)))
			allocate(Hg(itl(i))%col(nnz(i))); 
			! rowpntr size = pntr(m1s(i)+1), one less sectors 
			! k=m1 sector has no higher ksub sector to couple to.
			! +1 for last value =nnz+1
			Hg(itl(i))%srptr= pntr(m1s(i)+1) + 1; ! size of row pntr
			allocate(Hg(itl(i))%rowpntr(pntr(m1s(i)+1) +1))
			! set last value of rowptr
			Hg(itl(i))%rowpntr(pntr(m1s(i)+1) +1) = nnz(i) + 1;
		endif
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
					if (Hg(it2)%dense) then
						do i=1,ntot
							do j=1,n-k
								irow = pntr(k+1) + i - 1;
								icol = map(i,j)! col indx
								eig(it2)%evec(irow,icol) = val
							enddo
						enddo
					else
						Hg(it2)%dat(ind:indf) = val;
						! transpose map to get col first for reshape
						Hg(it2)%col(ind:indf) = 
     .    			reshape(transpose(map), (/ntot*(n-k)/))
						Hg(it2)%rowpntr(ptr:ptrf) =
     .       	(/ (ptr0 + (i-1)*(n-k)+1, i=1,ntot) /)
					endif
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
				if(k <= m1s(i2)) then
					if (Hg(itl(i2))%dense)then
						kdw = (ms(i2) - k)*dw; ! full value
						do i=1,ntot
							irow = ind1(i)
							eig(it2)%evec(irow,irow) = kdw
						enddo
					else
						kdw = (ms(i2) - k)*dw/2.0d0; ! half value
						Hg(itl(i2))%dat(ind:indf) = kdw;
						!Hg(itl(i2))%row(ind:indf) = ind1(:);
						Hg(itl(i2))%col(ind:indf) = ind1(:);
						ind = indf + 1;
						! shift pntrs, by one for each prev row
						Hg(itl(i2))%rowpntr(ptr:ptrf) = 
     .   		Hg(itl(i2))%rowpntr(ptr:ptrf) + 
     .   			(/(i,i=0,ntot-1)/)
     			endif
				endif
			end do
			ptr = ptrf+1	
			ptr0 = ptr0 + ntot*(n-k+1); ! shift for nxt k-blocks		
		else
			ptr = ptrf+1		
			ptr0 = ptr0 + ntot*(n-k); ! shift for nxt k-blocks		
		end if ! detuning

	write(*,*)	"--- 1 "		

	end do ! k

	write(*,*)	"--- 2 "		


	return
	end subroutine MakeHgMulti
!------------------------------------------































!------------------------------------------
	! when no detuning ( and all w0i=w0)
	subroutine MakeHgMulti1(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	! saves the Hamiltonians in a format similar to CSR,
	!	pointers to k-sub sectors in rowpntr, 
	!	and single val per sector in %dat
	use modmain, only: basis,na,nx,ibs,dna,dnx,
     .	 g,dw,detuning,Hg,eig,mapb,mapt,smalln,sameg
	implicit none
	integer(kind=1), intent(in) :: nt
	integer(kind=1), dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer(kind=4):: indf,i,j,ind,ntot,itype,indi,irow,icol
	integer(kind=1):: ib,k,m,m1,n,ibl,i2,it2
	double precision:: val,kdw
	integer(kind=4), allocatable, dimension(:):: pntr
	integer(kind=4), allocatable, dimension(:):: ind1
	integer(kind=4), allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,itl
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

	nnz=0;
	do k=0,m1-1
		ntot = pntr(k+2) - pntr(k+1)
		do i=1,nt
			if(k < m1s(i)) nnz(i) = nnz(i) + ntot*(n-k);
		end do
	end do

	! hilbert space dimension
	itl = 0;		
	do i=1,nt
		itl(i) = mapt%map(itlist(i));! map for the location of itype
		Hg(itl(i))%xst = .true.
		Hg(itl(i))%n = n;
		Hg(itl(i))%m = ms(i);
		Hg(itl(i))%m1 = m1s(i);
		Hg(itl(i))%ntot = pntr(m1s(i)+2);
		Hg(itl(i))%nnz = nnz(i);
		if (pntr(m1s(i)+2) .le. smalln) then
			Hg(itl(i))%dense = .true.
			! if true, save hamiltonian in eig(itl(i))%evec
		else
			! set this in case it is true from a prev run
			Hg(itl(i))%dense = .false. 
		endif		
		!write(*,*)"it, Hg(itl(i))%dense",itl(i),Hg(itl(i))%dense
	end do

	! allocate memory to Hg(itype)
	do i=1,nt
		if(allocated(Hg(itl(i))%col)) then
			deallocate(Hg(itl(i))%col)
			deallocate(Hg(itl(i))%rowpntr)
			deallocate(Hg(itl(i))%dat)
			deallocate(Hg(itl(i))%spntr)
		endif
		if (Hg(itl(i))%dense) then ! dense
			ntot = pntr(m1s(i)+2)
			if(allocated(eig(itl(i))%evec))deallocate(eig(itl(i))%evec)
			if(allocated(eig(itl(i))%eval))deallocate(eig(itl(i))%eval)
			allocate(eig(itl(i))%evec(ntot,ntot))
			eig(itl(i))%evec(:,:) = 0.0d0 ! initialise = 0
			allocate(eig(itl(i))%eval(ntot))
		else	 ! sparse
			allocate(Hg(itl(i))%spntr(m1s(i)+1));
			allocate(Hg(itl(i))%dat(m1s(i)));			
			allocate(Hg(itl(i))%col(nnz(i))); 
			! rowpntr size = pntr(m1s(i)+1), one less sectors 
			! k=m1 sector has no higher ksub sector to couple to.
			! +1 for last value =nnz+1
			Hg(itl(i))%srptr= pntr(m1s(i)+1) + 1; ! size of row pntr
			allocate(Hg(itl(i))%rowpntr(pntr(m1s(i)+1) +1))
			! set last value of rowptr
			Hg(itl(i))%rowpntr(pntr(m1s(i)+1) +1) = nnz(i) + 1;
			! set last value of pointer to (k-sub) sectors in rowpntr
			!write(*,*)"shape(Hg(it2)%spntr)=",shape(Hg(itl(i))%spntr)
			!write(*,*)"shape(pntr(1:m1s(i)+1))=",shape(pntr(1:m1s(i)+1))
			Hg(itl(i))%spntr(:) = pntr(1:m1s(i)+1) +1; ! upto k=m1-1 sectors
		endif
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
					if (Hg(it2)%dense) then
						do i=1,ntot
							do j=1,n-k
								irow = pntr(k+1) + i - 1;
								icol = map(i,j)! col indx
								eig(it2)%evec(irow,icol) = val
							enddo
						enddo
					else
						Hg(it2)%dat(k+1) = val;
						! transpose map to get col first for reshape
						Hg(it2)%col(ind:indf) = 
     .    			reshape(transpose(map), (/ntot*(n-k)/))
						Hg(it2)%rowpntr(ptr:ptrf) =
     .       	(/ (ptr0 + (i-1)*(n-k)+1, i=1,ntot) /)
					endif
				endif
			end do
			ind = indf+1; ! advance the index
		endif ! k < m1
		
		ptr = ptrf+1		
		ptr0 = ptr0 + ntot*(n-k); ! shift for nxt k-blocks		

	end do ! k

	

	return
	end subroutine MakeHgMulti1
!------------------------------------------



	
	end module


	
