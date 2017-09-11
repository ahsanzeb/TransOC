
	module hamiltonian
	implicit none
	public :: mkHamilt
	private:: MakeHgMulti, MakeHgMulti1, HgMap

	contains

!------------------------------------------
	subroutine mkHamilt()
	use modmain, only: mapt,Hg,detuning,sameg,debug,nog
	implicit none
	! local
	integer:: nt,ib
	integer:: grouptb(5,13)
	
		grouptb = mapt%grouptb;
		do ib=1,5
			nt = mapt%ntb(ib)
			if (nt > 0) then ! something to calculate or not?
				if(  nog    ) then 
					! H is diagonal, store evec sparse, later Ht.Uf sparse
					call DiagHgMulti(grouptb(ib,1:nt),nt)
				elseif ((.not. detuning) .and. sameg) then
					! use better storage format (and matvec routines for iter diag)
					!if(debug) write(*,*) "hamilt: better = T"
					call MakeHgMulti1(grouptb(ib,1:nt),nt)
				elseif (detuning) then
					!if(debug) write(*,*) "hamilt: better = F"
					write(*,*) nt, grouptb(ib,1:nt)					
					! simple CSR format to store sparse
					call MakeHgMulti(grouptb(ib,1:nt),nt)
					! CSR format, no detuning.
					!	call MakeHgMulti0(grouptb(ib,1:nt),nt)
					!write(*,*) 'done... '	
				endif
			endif
		end do

	!if(debug) write(*,*)"hamilt: Hg(:)%xst ============="
	!if(debug) write(*,*)	Hg(:)%xst
	!if(debug) write(*,*)	Hg(1)%rowpntr
	
	return
	end	subroutine mkHamilt
!------------------------------------------

!------------------------------------------
	subroutine HgMap(ibl,n,k,ntot,map)
	use modmain, only: basis,isk
	use lists, only: Complement,SortedInsert !SortedInsert=Append+Sort
	use basisstates, only: LexicoIndex
	implicit none
	integer, intent(in) ::n,k
	integer, intent(in) :: ibl
	integer, intent(in) :: ntot
	integer, dimension(ntot,n-k), intent(out):: map
	! local
	integer:: nk,i,j
	integer(kind=isk), dimension(k):: set1
	integer(kind=isk), dimension(k+1):: set2
	integer(kind=isk), dimension(n-k):: flipsites
	integer:: flip,k1
	integer(kind=isk), allocatable :: sites(:)
	
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
				!write(*,*)"ham: flip=",flip

				call SortedInsert(set1,k,flip,set2)
				!write(*,*)"ham: set2=",set2
				map(i,j) = LexicoIndex(set2,n,k1)
				!write(*,*)"ham:k,n, i,j,map(i,j)",k,n,i,j,map(i,j)
			end do
		end do
	endif
	return
	end subroutine HgMap
!------------------------------------------
	subroutine MakeHgMulti(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	! saves the Hamiltonians in CSR format
	use modmain, only: basis,na,nx,ibs,dna,dnx,isk,debug,
     .	 g,dw,Hg,eig,mapb,mapt,smalln,sameg,ndsec
	implicit none
	integer, intent(in) :: nt
	integer, dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer:: n,k

	integer:: indf,i,j,ind,ntot,itype,indi,irow,icol
	integer:: ib,m,m1,ibl,i2,it2
	double precision:: val,kdw
	integer, allocatable, dimension(:):: pntr
	integer, allocatable, dimension(:):: ind1
	integer, allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,nnzg,nnzd,itl
	integer :: ptr,ptrf,ptr0,nev,indx1,indx2
	!character*1000:: fmt

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

	write(*,*) " n = ",n
	
	do i=1,nt,1
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

	!if(debug) write(*,*)"ham: pntr=",pntr

	nnzg=0; nnzd=0;
	do k=0,m1-1 ! exclude the last k-sub sector
		ntot = pntr(k+2) - pntr(k+1)
		do i=1,nt
			if(k < m1s(i)) nnzg(i) = nnzg(i) + ntot*(n-k);
		end do
	end do

	!	number of diagonal elements
	do i=1,nt
		if(m1s(i) < ms(i)) then
			nnzd(i) = pntr(m1s(i)+2); ! total no of basis states
		else ! m1s(i)=ms(i)
			nnzd(i) = pntr(m1s(i)+1); ! last ksub sector dig elem = (m-k)*dw=0
		endif
	end do
	! total number of non zero elements
	nnz = nnzg + nnzd;

	! hilbert space dimension
	itl = 0;		
	do i=1,nt
		ntot = pntr(m1s(i)+2)
		it2 = mapt%map(itlist(i));! map for the location of itype
		itl(i) = it2;	
		Hg(it2)%xst = .true.
		Hg(it2)%n = n;
		Hg(it2)%m = ms(i);
		Hg(it2)%m1 = m1s(i);
		Hg(it2)%ntot = pntr(m1s(i)+2);
		Hg(it2)%nnz = nnz(i);
		if (pntr(m1s(i)+2) .le. smalln) then
			Hg(it2)%dense = .true.
			Hg(it2)%nev = ntot; ! ncv irrelevant for direct
			! if true, save hamiltonian in eig(itl(i))%evec
		else
			! set this in case it is true from a prev run
			Hg(it2)%dense = .false.
			nev = pntr(min(ndsec,m1s(i))+2); ! ~ nsnev degen sectors
			Hg(it2)%nev = nev;
			Hg(it2)%ncv = min(2*nev, ntot);
		endif		
	end do

	! allocate memory to Hg(itype)
	do i=1,nt
		it2 = itl(i);
		if(allocated(Hg(it2)%col)) then
			deallocate(Hg(it2)%col)
			deallocate(Hg(it2)%rowpntr)
			deallocate(Hg(it2)%dat)
		endif
		if (Hg(it2)%dense) then ! dense
			ntot = pntr(m1s(i)+2)
			eig(it2)%ntot = ntot;
			eig(it2)%n1 = ntot;
			if(allocated(eig(it2)%evec))deallocate(eig(it2)%evec)
			if(allocated(eig(it2)%eval))deallocate(eig(it2)%eval)
			allocate(eig(it2)%evec(ntot,ntot))
			eig(it2)%evec(:,:) = 0.0d0 ! initialise = 0
			allocate(eig(it2)%eval(ntot))
		else	 ! sparse
			allocate(Hg(it2)%dat(nnz(i)));
			allocate(Hg(it2)%col(nnz(i))); 
			! rowpntr size = pntr(m1s(i)+2), all ksub sectors 
			Hg(it2)%srptr= nnzd(i) + 1; ! size of row pntr
			allocate(Hg(it2)%rowpntr(nnzd(i) +1))
			! set last value of rowptr
			Hg(it2)%rowpntr(nnzd(i) +1) = nnz(i) + 1;
			!write(*,*) "i2, ntot, shape(Hg(it2)%rowpntr) = "
			!write(*,*) i, Hg(it2)%ntot, shape(Hg(it2)%rowpntr)
			!write(*,*) "nnz = ",nnz(i)
			!write(*,*) "shape(Hg(it2)%col)",shape(Hg(it2)%col)

		endif
	end do



			
	ind = 1; ptr=1; ptr0=0;
	! calc maps for diff k and combine to get full matrix
	do k=0,m1 ! =m1 for detuning case, to include diag elem
		ntot = pntr(k+2) - pntr(k+1);
		if(allocated(ind1))deallocate(ind1)
		allocate(ind1(ntot));
		ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);

		if (k<m1) then
			!write(*,*)" =====> k = ",k
			! calc map from k to k+1 up spin sector
			if(allocated(map))deallocate(map)
			allocate(map(ntot,n-k));
			call HgMap(ibl,n,k,ntot,map)
		endif

		!****************************************		
			! assign values to final matrix
			! to use if k < m1s()		
			indf=ind+ntot*(n-k+1)-1;
			ptrf = ptr + ntot -1
			do i2=1,nt
				it2 = itl(i2);
				!==================================
				if(k < m1s(i2)) then ! diagonal + off-diagonal
					!write(*,*)" =====> i2, k<m1s(i2)",i2, k,m1s(i2)
					!if(debug)write(*,*)"ham: k < m1s(i2)",k, m1s(i2)
					val = dsqrt((ms(i2)-k)*1.d0)*g;
					if (Hg(it2)%dense) then
						!---------------------------------
						!if(debug)write(*,*)"ham: dense"
						do i=1,ntot
							do j=1,n-k
								irow = pntr(k+1) + i;
								icol = pntr(k+2) + map(i,j)! col indx
								eig(it2)%evec(irow,icol) = val
							enddo
						enddo					
						kdw = (ms(i2) - k)*dw; ! full value
						do i=1,ntot
							irow = ind1(i)
							eig(it2)%evec(irow,irow) = kdw
						enddo
						!---------------------------------
					else ! sparse
						kdw = (ms(i2) - k)*dw/2.0d0; ! half value
						indx1 = ind;
						!write(*,*)"=====> col ========="
						!write(*,*)	map	
						!write(*,*)"=====> ind ========="
						!write(*,*)ind1		
						do i=1,ntot
							indx2 = indx1 + n-k - 1;
							Hg(it2)%dat(indx1:indx2) = val;
							Hg(it2)%col(indx1:indx2) = pntr(k+2) + map(i,:)
							indx2 = indx2 + 1; ! for diagonal element
							Hg(it2)%dat(indx2) = kdw;
							Hg(it2)%col(indx2) = ind1(i)
							!write(*,*)"indx2, ind1(i) = ",indx2, ind1(i)
							indx1 = indx2 + 1;
						enddo
						Hg(it2)%rowpntr(ptr:ptrf) =
     .         (/ (ptr0 + (i-1)*(n-k+1) + 1, i=1,ntot) /)
						!-----------------------------------     	
					endif
				! k==m1s(i2) .and. k < ms(i2), otherwise: (ms(i2) - k)*dw = 0
				elseif(k==m1s(i2) .and. k < ms(i2)) then ! only diagonal
					!write(*,*)" =====> i2, k==m1s(i2)",i2, k,m1s(i2)
				!-----------------------------------------
				! store half of actual values to avoid double counting
				! in mat vec multiplication in module diag, subroutine matvec()
				! ===> put a switch to control this behaviour in case
				! a direct eigensolver is to be used instead of arpack
				!-----------------------------------------
				! this might apply to only some i2 cases,
				! so dont affect ptr0 etc.
					if (Hg(it2)%dense)then
						kdw = (ms(i2) - k)*dw; ! full value
						do i=1,ntot
							irow = ind1(i)
							eig(it2)%evec(irow,irow) = kdw
						enddo
					else ! sparse
						!write(*,*) " ind1 = ",ind1
						!write(*,*) " ptr0+ntot = ",ptr0+ntot
						kdw = (ms(i2) - k)*dw/2.0d0; ! half value
						Hg(it2)%dat(ind:ind+ntot-1) = kdw;
						Hg(it2)%col(ind:ind+ntot-1) = ind1(:);
						Hg(it2)%rowpntr(ptr:ptr+ntot-1) = 
     .      (/(i,i=ptr0+1, ptr0+ntot)/)
     			endif
				endif ! k < m1s()
				!==================================			
			end do ! i2
			!****************************************		
		
			! advance the indexes
			ind = indf+1; 
			ptr = ptrf+1	;
			ptr0 = ptr0 + ntot*(n-k+1); ! shift for nxt k-blocks					
		
	end do ! k

	!write(*,*)"ham: ccccc  666666666666666"

	return
	end subroutine MakeHgMulti
!------------------------------------------

























!------------------------------------------
	! if no detuning
	subroutine MakeHgMulti0(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	! saves the Hamiltonians in CSR format
	use modmain, only: basis,na,nx,ibs,dna,dnx,isk,debug,
     .	 g,dw,Hg,eig,mapb,mapt,smalln,sameg,ndsec
	implicit none
	integer, intent(in) :: nt
	integer, dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer:: n,k

	integer:: indf,i,j,ind,ntot,itype,indi,irow,icol
	integer:: ib,m,m1,ibl,i2,it2
	double precision:: val,kdw
	integer, allocatable, dimension(:):: pntr
	integer, allocatable, dimension(:):: ind1
	integer, allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,nnzg,nnzd,itl
	integer :: ptr,ptrf,ptr0,nev
	!character*1000:: fmt

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

	! pointers for start index
	if(allocated(pntr))deallocate(pntr)
	allocate(pntr(m1+2));
	pntr(:) = basis(ibl)%pntr(1:m1+2) ! is this already calculated?

	!if(debug) write(*,*)"ham: pntr=",pntr

	nnz = 0
	do k=0,m1-1
		ntot = pntr(k+2) - pntr(k+1)
		do i=1,nt
			if(k < m1s(i)) nnz(i) = nnz(i) + ntot*(n-k);
		end do
	end do

	! hilbert space dimension
	itl = 0;		
	do i=1,nt
		ntot = pntr(m1s(i)+2)
		itl(i) = mapt%map(itlist(i));! map for the location of itype
		Hg(itl(i))%xst = .true.
		Hg(itl(i))%n = n;
		Hg(itl(i))%m = ms(i);
		Hg(itl(i))%m1 = m1s(i);
		Hg(itl(i))%ntot = pntr(m1s(i)+2);
		Hg(itl(i))%nnz = nnz(i);
		if (pntr(m1s(i)+2) .le. smalln ) then
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
	do k=0,m1-1
		ntot = pntr(k+2) - pntr(k+1);
		if(allocated(ind1))deallocate(ind1)
		allocate(ind1(ntot));
		ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);
		! calc map from k to k+1 up spin sector
		if(allocated(map))deallocate(map)
		allocate(map(ntot,n-k));
		call HgMap(ibl,n,k,ntot,map)
		!****************************************		
			! assign values to final matrix
			! to use if k < m1s()		
			indf=ind+ntot*(n-k)-1; 
			ptrf = ptr + ntot -1
			do i2=1,nt
				it2 = itl(i2);
				!==================================
				if(k < m1s(i2)) then
					if (Hg(it2)%dense) then
						do i=1,ntot
							do j=1,n-k
								irow = pntr(k+1) + i;
								icol = pntr(k+2) + map(i,j)! col indx
								eig(it2)%evec(irow,icol) = val
							enddo
						enddo					
					else ! sparse
						Hg(it2)%dat(ind:indf) = val;
						! transpose map to get col first for reshape
						Hg(it2)%col(ind:indf) = pntr(k+2) + 
     .    			reshape(transpose(map), (/ntot*(n-k)/))
						Hg(it2)%rowpntr(ptr:ptrf) =
     .       	(/ (ptr0 + (i-1)*(n-k)+1, i=1,ntot) /)
					endif						
				endif ! k < m1s()
				!==================================			
			end do ! i2
			!****************************************		
	
			! advance the indexes
			ind = indf+1; 
			ptr = ptrf+1	;
			ptr0 = ptr0 + ntot*(n-k); ! shift for nxt k-blocks					
		
	end do ! k

	return
	end subroutine MakeHgMulti0


























!------------------------------------------
	! when no detuning ( and all w0i=w0)
	subroutine MakeHgMulti1(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	! saves the Hamiltonians in a format similar to CSR,
	!	pointers to k-sub sectors in rowpntr, 
	!	and single val per sector in %dat
	use modmain, only: basis,na,nx,ibs,dna,dnx,ndsec,isk,
     .	 g,dw,detuning,Hg,eig,mapb,mapt,smalln,sameg
	implicit none
	integer, intent(in) :: nt
	integer, dimension(nt), intent(in) :: itlist
	!integer(kind=4), intent(out) :: nnz ! no of nonzero elements
	! local
	integer:: n,k
	integer:: indf,i,j,ind,ntot,itype,indi,irow,icol
	integer:: ib,m,m1,ibl,i2,it2
	double precision:: val,kdw
	integer, allocatable, dimension(:):: pntr
	integer, allocatable, dimension(:):: ind1
	integer, allocatable, dimension(:,:):: map
	integer, dimension(nt) :: ms,m1s,nnz,itl
	integer :: ptr,ptrf,ptr0, sr,nev
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
		ntot = pntr(m1s(i)+2)
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
			nev = pntr(min(ndsec,m1s(i))+2); ! ~ nsnev degen sectors
			if(nev > int(ntot/2)) nev = int(ntot/2)
			Hg(itl(i))%nev = nev;
			Hg(itl(i))%ncv = min(2*nev, ntot);

			!write(*,*) "	ham:	 n, nev,ncv = ",ntot,nev,min(2*nev, ntot)

		endif		
		!write(*,*)"it, Hg(itl(i))%dense",itl(i),Hg(itl(i))%dense
	end do


	!write(*,*) "	ham: itl(:)=",		itl



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
			sr = pntr(m1s(i)+1) +1;
			!write(*,*)"ham: sr = ",sr
			allocate(Hg(itl(i))%rowpntr(sr))
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
		if(allocated(ind1))deallocate(ind1)
		allocate(ind1(ntot));
		ind1(:) = pntr(k+1) + (/ (i, i=1,ntot,1) /);
		if (k<m1) then
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
				it2 = itl(i2);
				if(k < m1s(i2)) then
					val = dsqrt((ms(i2)-k)*1.d0)*g;
					if (Hg(it2)%dense) then
						do i=1,ntot
							do j=1,n-k
								irow = pntr(k+1) + i;
								icol = pntr(k+2) + map(i,j)! col indx
								!if(irow<=0)write(*,*)"**************row*************"
								!if(icol<=0)then
								!	write(*,*)"**************col*************"
								!	write(*,*)map
								!	write(*,*)"ibl,n,k,ntot",ibl,n,k,ntot,map
								!	stop			
								!endif
								eig(it2)%evec(irow,icol) = val
							enddo
						enddo
					else
						Hg(it2)%dat(k+1) = val;
						! transpose map to get col first for reshape
						Hg(it2)%col(ind:indf) = pntr(k+2) + 
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






!------------------------------------------
	subroutine DiagHgMulti(itlist,nt)
	! makes Hg of types in itlist; all for same n, diff m.
	! Hamiltonians are diagonal (nog=T)
	!	save evec sparse, just numbers.
	! eval =  diag elem of H
	use modmain, only: basis,na,nx,ibs,dna,dnx,
     .	 dw,Hg,eig,mapb,mapt
	implicit none
	integer, intent(in) :: nt
	integer, dimension(nt), intent(in) :: itlist
	! local
	integer:: n,k

	integer:: indf,i,j,ind,ntot,itype
	integer:: ib,m,m1,ibl,i2,it2
	double precision:: kdw
	integer, allocatable, dimension(:):: pntr
	integer, dimension(nt) :: ms,m1s,nnz,itl

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

	! pointers for start index
	if(allocated(pntr))deallocate(pntr)
	allocate(pntr(m1+2));
	pntr(:) = basis(ibl)%pntr(1:m1+2) ! is this already calculated?

	nnz=0;
	!	total number of (diagonal) elements
	do i=1,nt
		if(m1s(i) < ms(i)) then
			nnz(i) = pntr(m1s(i)+2); ! total no of basis states
		else ! m1s(i)=ms(i)
			nnz(i) = pntr(m1s(i)+1); ! last ksub sector dig elem = (m-k)*dw=0
		endif
	end do

	! hilbert space dimension
	itl = 0;		
	do i=1,nt
		ntot = pntr(m1s(i)+2)
		it2 = mapt%map(itlist(i));! map for the location of itype
		itl(i) = it2;	
		Hg(it2)%xst = .true.
		Hg(it2)%n = n;
		Hg(it2)%m = ms(i);
		Hg(it2)%m1 = m1s(i);
		Hg(it2)%ntot = pntr(m1s(i)+2);
		Hg(it2)%nnz = nnz(i);
		Hg(it2)%dense = .true.
		Hg(it2)%nev = ntot;
	end do

	! allocate memory
	do i=1,nt
		it2 = itl(i);
		ntot = pntr(m1s(i)+2)
		eig(it2)%ntot = ntot;
		eig(it2)%n1 = ntot;
		eig(it2)%nsec = ntot;
		eig(it2)%n2 = ntot;
		if(allocated(eig(it2)%eval))deallocate(eig(it2)%eval)
		allocate(eig(it2)%eval(ntot))
		eig(it2)%eval = 0.0d0 
	end do

	ind = 1;
	! calc maps for diff k and combine to get full matrix
	do k=0,m1 ! =m1 for detuning case, to include diag elem
			ntot = pntr(k+2) - pntr(k+1);
			do i2=1,nt
				it2 = itl(i2);
				if(k <= m1s(i2) .and. k < ms(i2)) then ! only diagonal
					kdw = (ms(i2) - k)*dw;! full value
					eig(it2)%eval(ind:ind+ntot-1) = kdw;
				endif
			end do ! i2
			ind = ind+ntot; 		
	end do ! k

	return
	end subroutine DiagHgMulti
!------------------------------------------


	
	end module


	
