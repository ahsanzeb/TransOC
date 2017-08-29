
	module hamiltonian
	implicit none

	contains
!------------------------------------------
	subroutine HgMap(ib,n,k,ntot,map)
	use modmain, only: basis
	use lists, only: Complement,SortedInsert !SortedInsert=Append+Sort
	use basisstates, only: LexicoIndex
	implicit none
	integer(kind=1), intent(in) :: ib,n,k
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
			set1=basis(ib)%sec(k)%sets(i,:)
			!write(*,*) "ib,k,i: ",ib,k,i,", set1=",set1
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
	subroutine makeHg(itype,nnz)
	! makes Hg of type index itype \in [1,13]
	use modmain, only: basis,na,nx,ibs,dna,dnx,
     .								g,dw,detuning,Hg,mapb,mapt
	implicit none
	integer, intent(in) :: itype
	integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer(kind=4) nnzg,nnzd,indf,i,j,ind,ntot
	integer(kind=1) ib,k,m,m1,n,ibl,itl
	double precision:: val,kdw
	integer(kind=4), allocatable, dimension(:):: pntr
	integer(kind=4), allocatable, dimension(:):: ind1
	integer(kind=4), allocatable, dimension(:,:):: map
	!character*1000:: fmt

	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible

	! ib: itype ===> which of 5 N case?
	ib = ibs(itype);
	ibl = mapb%map(ib); ! ib location
	itl = mapt%map(itype); ! itype location
	! pointers for start index
	if(allocated(pntr))deallocate(pntr)
	allocate(pntr(m1+2));
	pntr(:) = basis(ibl)%pntr(1:m1+2) ! only the relevant part
	!fmt = '(i5' // repeat (', 1x, i5', m1+2 - 1) // ')'
	!write(*,fmt),itype, pntr

	! find nnz: no of non-zero elem in Hg
	!	light-matter coupling terms:
	!	 nnzg = sumation_{k=0,m1-1}{^nC_k \times (n-k)}
	!	diagonal elem:
	!		if no site energy disorder or detuning:
	!			nnzd = tot no of basis
	!		else, nnzd=0	 
	nnzg = 0; nnzd=0;
	do k=0,m1-1
		ntot = pntr(k+2) - pntr(k+1)
		nnzg = nnzg + ntot*(n-k);
	end do
	if (detuning) nnzd = pntr(m1+2) ! total no of basis states
	nnz = nnzg + nnzd;

	Hg(itl)%ntot = pntr(m1+2); 	! ntot = hilbert space dimension
	write(*,*) "itype, ntot ",itype, pntr(m1+2)

	if(allocated(Hg(itl)%row)) then
		deallocate(Hg(itl)%row)
		deallocate(Hg(itl)%col)
		deallocate(Hg(itl)%dat)
		!call memory('d','i4',100,'makeHg')
	endif
	allocate(Hg(itl)%row(nnz))
	allocate(Hg(itl)%col(nnz))
	allocate(Hg(itl)%dat(nnz))
	!call memory('a','i',100,'makeHg')

	ind = 1;
	! calc maps for diff k and combine to get full matrix
	do k=0,m1-1,1
		ntot = pntr(k+2) - pntr(k+1);
		!write(*,*)"makeHg: n,m,k,m1,ntot = ",n,m,k,m1,ntot
		if(allocated(ind1))deallocate(ind1)
		allocate(ind1(ntot));
		ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);
		! calc map from k to k+1 up spin sector
		if(allocated(map))deallocate(map)
		allocate(map(ntot,n-k));
		call HgMap(ib,n,k,ntot,map)
		! assign values to final matrix		
		val = dsqrt((m-k)*1.d0)*g;
		indf=ind+ntot*(n-k);
		Hg(itl)%dat(ind:indf) = val;
		do i=1,ntot
			do j=1,n-k
				Hg(itl)%row(ind) = ind1(i);
				Hg(itl)%col(ind) = map(i,j);
				!Hg(itype)%dat(ind) = val;
				ind = ind + 1;
			end do
		end do
		
		! diagonal elements 
		kdw = (m - k)*dw;
		if (detuning) then
			indf=ind+ntot-1;
			Hg(itl)%dat(ind:indf) = kdw;
			Hg(itl)%row(ind:indf) = ind1(:);
			Hg(itl)%col(ind:indf) = ind1(:);
			ind = indf + 1;
		end if
	end do

	return
	end subroutine makeHg
!------------------------------------------

!------------------------------------------
	subroutine MakeHgMulti(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	use modmain, only: basis,na,nx,ibs,dna,dnx,
     .								g,dw,detuning,Hg,mapb,mapt
	implicit none
	integer, intent(in) :: nt
	integer, dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer(kind=4) indf,i,j,ind,ntot,itype,indi
	integer(kind=1) ib,k,m,m1,n,ibl,i2
	double precision:: val,kdw
	integer(kind=4), allocatable, dimension(:):: pntr
	integer(kind=4), allocatable, dimension(:):: ind1
	integer(kind=4), allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,nnzg,nnzd,itl
	
	!character*1000:: fmt


	! check if all types have the same n
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
	do i=2,nt,1
		m	=	 nx + dnx(itlist(i)); ! no of excitations
		ms(i) = m;
		m1s(i)=min(m,n); ! max no of up spins possible
	end do
	m1 = maxval(m1s) ! max m1

	! ib: itype ===> which of 5 N case?
	ib = ibs(itlist(1)); ! same for all cases in itlist
	ibl = mapb%map(ib); ! get the location of data for this ib

	! pointers for start index
	if(allocated(pntr))deallocate(pntr)
	allocate(pntr(m1+2));
	pntr(:) = basis(ibl)%pntr(1:m1+2) ! is this already calculated?



	nnzg=0; nnzd=0;
	do k=0,m1-1
		ntot = pntr(k+2) - pntr(k+1)
		do i=1,nt
			if(k<m1s(i)) nnzg(i) = nnzg(i) + ntot*(n-k);
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
		Hg(itl(i))%ntot = pntr(m1s(i)+2);
	end do

	! allocate memory to Hg(itype)
	do i=1,nt
		if(allocated(Hg(itl(i))%row)) then
			deallocate(Hg(itl(i))%row)
			deallocate(Hg(itl(i))%col)
			deallocate(Hg(itl(i))%dat)
		endif
		allocate(Hg(itl(i))%row(nnz(i)))
		allocate(Hg(itl(i))%col(nnz(i)))
		allocate(Hg(itl(i))%dat(nnz(i)))
	end do
	
	ind = 1;
	! calc maps for diff k and combine to get full matrix
	do k=0,m1-1,1
		ntot = pntr(k+2) - pntr(k+1);
		if(allocated(ind1))deallocate(ind1)
		allocate(ind1(ntot));
		ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);
		! calc map from k to k+1 up spin sector
		if(allocated(map))deallocate(map)
		allocate(map(ntot,n-k));
		call HgMap(ib,n,k,ntot,map)
		
		! assign values to final matrix			
		indf=ind+ntot*(n-k);
		indi = ind;
		do i2=1,nt
			ind = indi; ! save the index for i2 loop
			if(k < m1s(i2)) then
				val = dsqrt((ms(i2)-k)*1.d0)*g;
				Hg(itl(i2))%dat(ind:indf) = val;
				do i=1,ntot
					do j=1,n-k
						Hg(itl(i2))%row(ind) = ind1(i);
						Hg(itl(i2))%col(ind) = map(i,j);
						ind = ind + 1;
					end do
				end do
			endif
		end do
		
		! diagonal elements 
		if (detuning) then
			indf=ind+ntot-1;
			indi = ind; ! save the index for i2 loop
			do i2=1,nt
				ind = indi;
				kdw = (ms(i2) - k)*dw;
				if(k < m1s(i2)) then
					Hg(itl(i2))%dat(ind:indf) = kdw;
					Hg(itl(i2))%row(ind:indf) = ind1(:);
					Hg(itl(i2))%col(ind:indf) = ind1(:);
					ind = indf + 1;
				endif
			end do
		end if ! detuning

	end do ! k


	return
	end subroutine MakeHgMulti
!------------------------------------------

	
	end module


	
