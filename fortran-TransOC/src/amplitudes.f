
	module amplitudes
	use modmultiply
	
	implicit none

	contains
! transition matrices and amplitudes
	subroutine CalAmp(ih,ic,is,rowc,nnz,n3,routine)
	use modmain, only: qt,mapt,maph,eig,itypes,psi
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ic
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	character(len=*), intent(in) :: routine
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i

	!write(*,*)"amp: ih,ic,is = ",ih,ic,is
	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
		itl = mapt%map(itypes(ih,ic)) !???????! location of final hilber space
		!write(*,*)"ih,ic,it,itl= ",ih,ic,itypes(ih,ic),itl
		
		n1=eig(itl)%n1 ! dim of final hilbert space
		n2=eig(itl)%n2

		!write(*,*) "n1,n2 = ", n1,n2
		
		allocate(HtUf(n3,n2))
		ia = maph(ih,ic); ! location of amplitudes
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		if(routine=='multiplyd') then
			call multiplyd(rowc,nnz,eig(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		elseif(routine=='multiplydc') then
			call multiplydc(rowc,nnz,eig(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		else
			write(*,*) "amplitudes: something wrong....!"
			stop
		endif
		if(allocated(qt(ia)%cs(ic,is)%amp))
     .						deallocate(qt(ia)%cs(ic,is)%amp)
		allocate(qt(ia)%cs(ic,is)%amp(n2))
		! set the size for possible future use
		!	DELETE THIS namp VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%namp = n2
		! multiply psi with HtUf to get amplitudes
		! psi should be a row vector; shape = 1 x n3

		!HtUf = 1.0d0;
		!write(*,*)"amp: n3,n1,n2=",n3,n1,n2
		!write(*,*)"amp: shape(psi),shape(HtUf) ",shape(psi),shape(HtUf) 
		!if(ih==0) then
		!	write(*,*)"==================="
		!	write(*,*)"HtUf=",HtUf
		!	write(*,*)"==================="
		!endif
		! testing HtUf(1,:);!
		qt(ia)%cs(ic,is)%amp=reshape(matmul(psi,HtUf),(/n2/)); ! both input dense
		deallocate(HtUf)
		! calculate amp^2 for degenerate sectors
		if(allocated(qt(ia)%cs(ic,is)%amp2))
     .					deallocate(qt(ia)%cs(ic,is)%amp2)
		allocate(qt(ia)%cs(ic,is)%amp2(eig(itl)%nsec))
		!	DELETE THIS nsec VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%nsec = eig(itl)%nsec !

		!write(*,*) "amp: ih,ic, ia: ic,is,itl=",ih,ic, ia, ic,is,itl
		!write(*,*)"eig(itl)%nsec,eig(itl)%ind",eig(itl)%nsec,eig(itl)%ind
		!write(*,*)"eval=",eig(itl)%eval
		!write(*,*)"............................"
		!write(*,*)"amp =",	qt(ia)%cs(ic,is)%amp
		call GetAmp2(qt(ia)%cs(ic,is)%amp, n2,
     .		qt(ia)%cs(ic,is)%amp2, eig(itl)%nsec, eig(itl)%ind)


	return
	end subroutine CalAmp
!---------------------------------------------


	subroutine CalAmp0(ih,ic,is,rowc,nnz,n3,col)
	use modmain, only: qt,mapt,maph,eig,itypes,psi
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ic
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	integer, dimension(nnz), intent(in) :: col
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
		itl = mapt%map(itypes(ih,ic)) !???????! location of final hilber space
		!write(*,*)"itypes(ih,ic), itl = ",itypes(ih,ic),itl

		
		n1=eig(itl)%n1 ! dim of final hilbert space
		n2=eig(itl)%n2
		allocate(HtUf(n3,n2))
		ia = maph(ih,ic); ! location of amplitudes
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		call multiply(rowc,col,nnz,eig(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
		if(allocated(qt(ia)%cs(ic,is)%amp))
     .						deallocate(qt(ia)%cs(ic,is)%amp)
		allocate(qt(ia)%cs(ic,is)%amp(n2))
		! set the size for possible future use
		!	DELETE THIS namp VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%namp = n2
		! multiply psi with HtUf to get amplitudes
		! psi should be a row vector; shape = 1 x n3
		!HtUf = 1.0d0;
		!write(*,*)"amp0: n3,n1,n2=",n3,n1,n2
		!write(*,*)"amp0: shape(psi),shape(HtUf) ",shape(psi),shape(HtUf) 
		qt(ia)%cs(ic,is)%amp=reshape(matmul(psi,HtUf),(/n2/)); ! both input dense
		deallocate(HtUf)
		! calculate amp^2 for degenerate sectors
		if(allocated(qt(ia)%cs(ic,is)%amp2))
     .					deallocate(qt(ia)%cs(ic,is)%amp2)
		allocate(qt(ia)%cs(ic,is)%amp2(eig(itl)%nsec))
		!	DELETE THIS nsec VARIABLE IF NOT NEEDED/USED.
		qt(ia)%cs(ic,is)%nsec = eig(itl)%nsec ! 


		!write(*,*) "amp: ih,ic, ia: ic,is,itl=",ih,ic, ia, ic,is,itl

		call GetAmp2(qt(ia)%cs(ic,is)%amp, n2,
     .		qt(ia)%cs(ic,is)%amp2,
     .		eig(itl)%nsec, eig(itl)%ind)


	return
	end subroutine CalAmp0


!----------------------------------------------
	! if no disorder in site w0 and g, and AlwaysLP,
	!	then permutation symmetry holds for all sites for a given hop,
	!	so calculate amplitudes for a single site and use it for all sites; 
	! similarly, calcualte amp2 once.
	subroutine GetAmp2(amp,ne,amp2,nsec,ind)
	! calculates total amplitude squared for degenerate sectors
	integer, intent(in) :: ne, nsec
	double precision,dimension(ne),intent(in):: amp
	integer,dimension(nsec+1),intent(in):: ind !start index of sectors
	double precision,dimension(nsec),intent(out):: amp2 !tot amp^2 of sec
	! local
	integer:: i,i1,i2,j

	amp2 = 0.0d0;
	do i=1,nsec !-1
		i1 = ind(i); i2=ind(i+1)-1
		!write(*,*) "ne, i,i1,i2 = ",ne, i,i1,i2
		!do j = i1,i2
		!	amp2(i) = amp2(i) + amp(j)
		!end do 
		amp2(i) = sum( amp(i1:i2)**2 )	
	end do

	if(sum(amp2)> 10000) then
		write(*,*) "amp:   ne, nsec,shape(ind),shape(amp) = "
		write(*,*) ne,nsec, shape(ind),shape(amp)
		write(*,*) "amp:   ind = ",ind
		write(*,*) "amp:   amp = ",amp
		write(*,*) "....................."
		write(*,*) "amp:   amp2 = ",amp2
		stop
	endif
	
	return
	end subroutine GetAmp2
!---------------------------------------

	end 	module amplitudes
