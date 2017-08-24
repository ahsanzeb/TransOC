
	module hamiltonian
	use modmain
	implicit none




!------------------------------------------
	subroutine HgMap(ib,n,k,ntot,map)
	use modmain, only: basis, sites,na
	use lists, only: Complement,Append,Sort,LexicoIndex
	implicit none
	integer, intent(in) :: ib,n,k,ntot
	integer(kind=4), dimension(ntot,n-k), intent(out):: map
	! local
	integer:: nak,flip,i,j
	integer(kind=1), dimension(k):: set1
	integer(kind=1), dimension(k+1):: set2
	integer(kind=1), dimension(n-k):: flipsites
		
	if (k==0) then
		map(1,:) = (/ (i,i=1,n,1)/)
	elseif(k==n) then
		! all map to the single state with all up
		map(:,1) = 1;
	else
		nak = na-k;
		do i=1,ntot
		! sec(k) are defined for k>0,
		! so index match with k
			set1=basis(ib)%sec(k)%sets(i,:)
			call Complement(sites,na,set1,k,flipsites)
			do j=1,nak
				flip = flipsites(j);
				call Append(set1,k,flip,set2)
				call Sort(set2,k+1);
				map(i,j) = LexicoIndex(set2,n,k+1)
			end do
		end do
	endif
	return
	end subroutine HgMap
!------------------------------------------
	subroutine makeHg(itype)
	! makes Hg of type index itype \in [1,13]
	use modmain, only: basis,na,nx,nainds,dna,dnx,g,dw
	implicit none
	integer nnz,nnzg,nnzd,indf
	double precision:: val

	! N,m values of itype:
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible

	! ib: itype ===> which of 5 N case?
	ib = nainds(itype);
	! pointers for start index
	pntr(:) = basis(ib)%pntr(1:m1+2) ! only the relevant part

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
	! allocate memory to Hg(itype)
	if(allocated(Hg(itype)%row))deallocate(Hg(itype)%row)
	if(allocated(Hg(itype)%row))deallocate(Hg(itype)%col)
	if(allocated(Hg(itype)%row))deallocate(Hg(itype)%dat)
	allocate(Hg(itype)%row(nnz))
	allocate(Hg(itype)%row(nnz))
	allocate(Hg(itype)%row(nnz))

	ind = 0;
	! calc maps for diff k and combine to get full matrix
	do k=0,m1-1,1
		ntot = pntr(k+2) - pntr(k+1);
		ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);
		! calc map from k to k+1 up spin sector
		call HgMap(ib,n,k,ntot,map)
		! assign values to final matrix		
		val = dsqrt((m-k)*1.d0)*g;
		indf=ind+ntot*(n-k);
		Hg(itype)%dat(ind:indf) = val;
		do i=1,ntot
			do j=1,n-k
				Hg(itype)%row(ind) = ind1(i);
				Hg(itype)%col(ind) = map(i,j);
				!Hg(itype)%dat(ind) = val;
				ind = ind + 1;
			end do
		end do
		
		! diagonal elements 
		kdw = (m - k)*dw;
		if (detuning) then
			indf=ind+ntot;
			Hg(itype)%dat(ind:indf) = kdw;
			Hg(itype)%row(ind:indf) = ind1(:);
			Hg(itype)%col(ind:indf) = ind1(:);
			ind = indf + 1;
		end if
		
	end do
!------------------------------------------


	



	





	end subroutine makeHg
!------------------------------------------

	
	end module


	
