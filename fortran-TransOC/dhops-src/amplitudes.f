
	module amplitudes
	use modmultiply
	
	implicit none

! CalAmp and CalAmp0 use different routines for Ht.Uf multiplication

	contains
! transition matrices and amplitudes
	subroutine CalAmp(ih,ici,is,rowc,nnz,n3,routine)
	use modmain, only: qt,mapt,maph,eig,itypes,psi,ipsi,nog,mapc
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ici
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	character(len=*), intent(in) :: routine
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i,it1,ic
	double precision, allocatable, dimension(:,:):: Hte ! Ht in eigenbasis

	!write(*,*)"amp: ih,ic,is = ",ih,ic,is
	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	itl = mapt%map(itypes(ih,ici)) !???????! location of final hilber space
	!write(*,*)"ih,ic,it,itl= ",ih,ic,itypes(ih,ic),itl
	ia = maph(ih,ici); ! location of amplitudes
	ic = mapc(ih,ici); ! location of ici 
	! 08Dec2017: just using ici for ic and ic for icl becasue otherwise I will have to change ic to icl in too many places in these amp,amp0 routines... this is becuase for PermSym case, I now want to use Dhops amp for Phihops amp, so need to give location for storing the amp that both Dhop and phihops processes know later when retrieving qt%cs%amp for rate calculations. 
	n1=eig(itl)%n1 ! dim of final hilbert space
	n2=eig(itl)%n2
	
	call allocqt(ia,ic,is,itl,n2)

	if(nog) then
		!write(*,*)"========> ih,ic,is = ",ih,ic,is 
		! use efficient mat vec multiplications for diagonal Uf
		qt(ia)%cs(ic,is)%amp= 0.0d0; ! psi.Ht.Uf
		! just index: psi = ipsi; Uf = Identity
		if(routine=='multiplyd') then
			do i=1,nnz ! diagonal but not full diagonal, selected rows/cols
				if (rowc(i) == ipsi) then ! if ipsi is included in rowc
					!write(*,*)"========> amp: i, rowc(i)=psi=",i,ipsi
					qt(ia)%cs(ic,is)%amp(rowc(i)) = 1.0d0
					exit
				endif
			end do
		elseif(routine=='multiplydc') then
			!write(*,*)"amp: shape(rowc), ipsi = ",shape(rowc), ipsi
			qt(ia)%cs(ic,is)%amp(rowc(ipsi)) = 1.0d0
		else
			write(*,*) "amplitudes: something wrong....!"
			stop
		endif
		qt(ia)%cs(ic,is)%amp2 = qt(ia)%cs(ic,is)%amp; ! 0.0, or +1.0 
	else
		allocate(HtUf(n3,n2))
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
			! multiply psi with HtUf to get amplitudes
			! psi should be a row vector; shape = 1 x n3
			! testing HtUf(1,:);!
			qt(ia)%cs(ic,is)%amp=reshape(matmul(psi,HtUf),(/n2/)); ! both input dense
			!write(*,*) "amp: ih,ic, ia, is,itl=",ih,ic, ia, is,itl
			call GetAmp2(qt(ia)%cs(ic,is)%amp, n2,
     .		qt(ia)%cs(ic,is)%amp2, eig(itl)%nsec, eig(itl)%ind)
		deallocate(HtUf)
	endif

	return
	end subroutine CalAmp
!---------------------------------------------


	subroutine CalAmp0(ih,ici,is,rowc,nnz,n3,col)
	use modmain, only: qt,mapt,maph,eig,itypes,psi,ipsi,nog,mapc
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ici
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	integer, dimension(nnz), intent(in) :: col
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i,it1,ic
	double precision, allocatable, dimension(:,:):: Hte ! Ht in eigenbasis

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	itl = mapt%map(itypes(ih,ici)) !???????! location of final hilber space
	!write(*,*)"itypes(ih,ic), itl = ",itypes(ih,ic),itl
	ia = maph(ih,ici); ! location of amplitudes
	ic = mapc(ih,ici); ! confused about ic and ici? see comments in calamp above...
	n1=eig(itl)%n1 ! dim of final hilbert space
	n2=eig(itl)%n2
	call allocqt(ia,ic,is,itl,n2)

	if(nog) then
		qt(ia)%cs(ic,is)%amp = 0.0d0
		do i=1,nnz
			if (rowc(i) == ipsi) then
				!write(*,*)"ih, ic,is = ",ih, ic,is
				!write(*,*)"ipsi, i, ir, ic = ",ipsi,i,rowc(i),col(i)
				qt(ia)%cs(ic,is)%amp(col(i)) = 1.0d0
				exit ! just a single entry
			endif
		end do
		qt(ia)%cs(ic,is)%amp2 = qt(ia)%cs(ic,is)%amp; ! 0.0, +1.0
	else
		allocate(HtUf(n3,n2))
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		call multiply(rowc,col,nnz,eig(itl)%evec,n1,n2,
     .							HtUf,n3,n2)

			! multiply psi with HtUf to get amplitudes
			! psi should be a row vector; shape = 1 x n3
			!HtUf = 1.0d0;
			!write(*,*)"amp0: n3,n1,n2=",n3,n1,n2
			!write(*,*)"amp0: shape(psi),shape(HtUf) ",shape(psi),shape(HtUf) 
			qt(ia)%cs(ic,is)%amp=reshape(matmul(psi,HtUf),(/n2/)); ! both input dense
			!write(*,*) "amp0: ih,ic, ia, is,itl=",ih,ic, ia, is,itl
			call GetAmp2(qt(ia)%cs(ic,is)%amp, n2,
     .		qt(ia)%cs(ic,is)%amp2,
     .		eig(itl)%nsec, eig(itl)%ind)
		deallocate(HtUf)  
	endif

	
	return
	end subroutine CalAmp0
!----------------------------------------------------------
	! cavity losses case
	subroutine CalAmpk(ih,ic,is,rowc,nnz,n3)
	use modmain, only: qt,mapt,maph,eig,itypes,psi,ipsi,nog
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ic
	integer, intent(in) :: nnz,n3
	double precision, dimension(nnz), intent(in)  :: rowc
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i,it1
	double precision, allocatable, dimension(:,:):: Hte ! Ht in eigenbasis

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	itl = mapt%map(itypes(ih,ic)) !???????! location of final hilber space
	ia = maph(ih,ic); ! location of amplitudes
		
	n1=eig(itl)%n1 ! dim of final hilbert space
	n2=eig(itl)%n2
	
	call allocqt(ia,ic,is,itl,n2)

	if(nog) then
		! use efficient mat vec multiplications for diagonal Uf
		qt(ia)%cs(ic,is)%amp= 0.0d0; ! psi.Ht.Uf
		! just index: psi = ipsi; Uf = Identity
		!qt(ia)%cs(ic,is)%amp(ipsi) = rowc(ipsi)
		if(ipsi<= n2) then ! only if ipsi has a poton, otherwise amp=0
			qt(ia)%cs(ic,is)%amp(ipsi) = rowc(ipsi)
		endif
		qt(ia)%cs(ic,is)%amp2 = qt(ia)%cs(ic,is)%amp**2;
	else
		allocate(HtUf(n3,n2))
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
		call multiplydk(rowc,nnz,eig(itl)%evec,n1,n2,
     .							HtUf,n3,n2)
			! multiply psi with HtUf to get amplitudes
			! psi should be a row vector; shape = 1 x n3
			! testing HtUf(1,:);!
			qt(ia)%cs(ic,is)%amp=reshape(matmul(psi,HtUf),(/n2/)); ! both input dense
			!write(*,*) "amp: ih,ic, ia, is,itl=",ih,ic, ia, is,itl
			call GetAmp2(qt(ia)%cs(ic,is)%amp, n2,
     .		qt(ia)%cs(ic,is)%amp2, eig(itl)%nsec, eig(itl)%ind)
		deallocate(HtUf)
	endif

	return
	end subroutine CalAmpk
!---------------------------------------------


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

	!-------------------------------------
	subroutine allocqt(ia,ic,is,itl,n2)
	use modmain, only: qt,eig
	implicit none
	integer, intent(in):: ia,ic,is,itl,n2
	if(allocated(qt(ia)%cs(ic,is)%amp))
     .						deallocate(qt(ia)%cs(ic,is)%amp)
	allocate(qt(ia)%cs(ic,is)%amp(n2))
	qt(ia)%cs(ic,is)%namp = n2
	! calculate amp^2 for degenerate sectors
	if(allocated(qt(ia)%cs(ic,is)%amp2))
     .					deallocate(qt(ia)%cs(ic,is)%amp2)
	allocate(qt(ia)%cs(ic,is)%amp2(eig(itl)%nsec))
	qt(ia)%cs(ic,is)%nsec = eig(itl)%nsec !

	return
	end subroutine allocqt
	!-------------------------------------

! transition amplitudes for bosonic states, complex in this case
	subroutine CalAmpBose(ih,ici,is,rowc,nnz,n3)
	use modmain, only: qt,mapt,maph,itypes,eig,psi,ipsi,nog,mapc,eigb
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ici
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	character(len=*), intent(in) :: routine
	! local
	double complex, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i,it1,ic
	!double complex, allocatable, dimension(:,:):: Hte ! Ht in eigenbasis

	!write(*,*)"amp: ih,ic,is = ",ih,ic,is
	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	itl = mapt%map(itypes(ih,ici)) !???????! location of final hilber space
	!write(*,*)"ih,ic,it,itl= ",ih,ic,itypes(ih,ic),itl
	ia = maph(ih,ici); ! location of amplitudes
	ic = mapc(ih,ici); ! location of ici 
	! 08Dec2017: just using ici for ic and ic for icl becasue otherwise I will have to change ic to icl in too many places in these amp,amp0 routines... this is becuase for PermSym case, I now want to use Dhops amp for Phihops amp, so need to give location for storing the amp that both Dhop and phihops processes know later when retrieving qt%cs%amp for rate calculations. 
	n1=eig(itl)%n1 ! dim of final hilbert space
	n2=eig(itl)%n2
	
	if(allocated(qtbc(ic)%amp)) deallocate(qtbc(ic)%amp)
	allocate(qtbc(ic)%amp(n2))
	qtbc(ic)%namp = n2
	! calculate amp^2 for degenerate sectors
	if(allocated(qtbc(ic)%amp2)) deallocate(qtbc(ic)%amp2)
	allocate(qtbc(ic)%amp2(eig(itl)%nsec))
	qtbc(ic)%nsec = eig(itl)%nsec ! nsec???? 

	allocate(HtUf(n3,n2))
		! NOTE: allocate qt(ih)%cs(:,:) in calling routine
		! sizes: Ht(n3 x n1) . Uf(n1 x n2) = HtUf(n3 x n2)
		! out: HtUf
	call multiplydz(rowc,nnz,eigb%evec,n1,n2,HtUf,n3,n2)

			! multiply psi with HtUf to get amplitudes
			! psi should be a row vector; shape = 1 x n3
			! testing HtUf(1,:);!
			qtbc(ic)%amp=reshape(matmul(psib,HtUf),(/n2/)); ! both input dense
			!write(*,*) "amp: ih,ic, ia, is,itl=",ih,ic, ia, is,itl
			call GetAmp2(qtbc(ic)%amp, n2,
     .		qtbc(ic)%amp2, eigb%nsec, eigb%ind)
		deallocate(HtUf)

	return
	end subroutine CalAmpBose
!---------------------------------------------


!------------------------------------------
	end 	module amplitudes
