
! transition matrices and amplitudes
	subroutine CalAmp(ih,ic,is,rowc,nnz,n3,routine,col)
	use modmain, only: qt,mapt,maph,hspace,itypes,psi
	implicit none
	integer(kind=1), intent(in) :: ih,is,ic,nnz
	integer, dimension(nnz), intent(in)  :: rowc
	integer, dimension(nnz), intent(in), optional :: col
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,n3,itl
	character(len=*) :: routine
	
	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
		itl = mapt%map(itypes(ih,ic)) !???????! location of final hilber space
		n1=hspace(itl)%n1 ! dim of final hilbert space
		n2=hspace(itl)%n2
!		n3=pntr1(m1+2) ! dim of initial hilbert space
	
		allocate(HtUf(n3,n2))
		ia = maph(ih,ic); ! location of amplitudes
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		if (present(col)) then
			call multiply(rowc,col,nnz,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		elseif(routine=="multiplyd") then
			call multiplyd(rowc,nnz,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		elseif(routine=="multiplydc") then
			call multiplydc(rowc,nnz,hspace(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		else
			write(*,*) "amplitudes: something wrong....!"
			stop
		endif
		!if(allocated(qt(ia)%cs(1,is)%amp))deallocate(qt(ia)%cs(1,is)%amp)
		allocate(qt(ia)%cs(ic,is)%amp(n2))
		! set the size for possible future use
		!	DELETE THIS namp VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%namp = n2
		! multiply psi with HtUf to get amplitudes
		! psi should be a row vector; shape = 1 x n3
		qt(ia)%cs(ic,is)%amp = matmul(psi,HtUf); ! both input dense
		deallocate(HtUf)

	end subroutine CalAmp




