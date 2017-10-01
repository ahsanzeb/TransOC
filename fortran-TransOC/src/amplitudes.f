
	module amplitudes
	use modmultiply
	
	implicit none

	contains
! transition matrices and amplitudes
	subroutine CalAmp(ih,ic,is,rowc,nnz,n3,routine)
	use modmain, only: qt,mapt,maph,eig,itypes,psi,ipsi,nog,master
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ic
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	character(len=*), intent(in) :: routine
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i,it1
	double precision, allocatable, dimension(:,:):: Hte ! Ht in eigenbasis

	!write(*,*)"amp: ih,ic,is = ",ih,ic,is
	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	itl = mapt%map(itypes(ih,ic)) !???????! location of final hilber space
	!write(*,*)"ih,ic,it,itl= ",ih,ic,itypes(ih,ic),itl
	ia = maph(ih,ic); ! location of amplitudes
		
	n1=eig(itl)%n1 ! dim of final hilbert space
	n2=eig(itl)%n2
	
	if(.not. master) call allocqt(ia,ic,is,itl,n2)

	if(nog) then
		!write(*,*)"========> ih,ic,is = ",ih,ic,is 
		! use efficient mat vec multiplications for diagonal Uf
		qt(ia)%cs(ic,is)%amp= 0.0d0; ! psi.Ht.Uf
		! just index: psi = ipsi; Uf = Identity
		if(routine=='multiplyd') then
			do i=1,nnz ! diagonal but not full diagonal, selected rows/cols
				if (rowc(i) == ipsi) then
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

		if (master) then
			it1 = mapt%map(1); ! inition/current hilbert space location
			allocate(Hte(eig(it1)%n2,n2))
			Hte = matmul(eig(it1)%evec,HtUf);
			! Hte is Ht in eigenbasis.
			!	split rhovt into degenerate sectors
			!	and do: Hte^T.rhovtsec.Hte OR Hte.rhovtsec.Hte^T ?
			! for all t-points on coarse grid.
			! get amp2. Multiply penalty function and get rates.
			call merates(ih,ic,is,Hte,eig(it1)%n2,n2,it1,itl)
			deallocate(Hte)
		else
			! multiply psi with HtUf to get amplitudes
			! psi should be a row vector; shape = 1 x n3
			! testing HtUf(1,:);!
			qt(ia)%cs(ic,is)%amp=reshape(matmul(psi,HtUf),(/n2/)); ! both input dense
			!write(*,*) "amp: ih,ic, ia, is,itl=",ih,ic, ia, is,itl
			call GetAmp2(qt(ia)%cs(ic,is)%amp, n2,
     .		qt(ia)%cs(ic,is)%amp2, eig(itl)%nsec, eig(itl)%ind)
		endif
		deallocate(HtUf)
	endif

	return
	end subroutine CalAmp
!---------------------------------------------


	subroutine CalAmp0(ih,ic,is,rowc,nnz,n3,col)
	use modmain, only: qt,mapt,maph,eig,itypes,psi,ipsi,nog,master
	implicit none
	integer, intent(in) :: is
	integer, intent(in) :: ih,ic
	integer, intent(in) :: nnz,n3
	integer, dimension(nnz), intent(in)  :: rowc
	integer, dimension(nnz), intent(in) :: col
	! local
	double precision, allocatable:: HtUf(:,:)
	integer :: ia,n1,n2,itl,i,it1
	double precision, allocatable, dimension(:,:):: Hte ! Ht in eigenbasis

	!-------------------------------------------------------
	! calculate transition amplitudes
	!-------------------------------------------------------
	itl = mapt%map(itypes(ih,ic)) !???????! location of final hilber space
	!write(*,*)"itypes(ih,ic), itl = ",itypes(ih,ic),itl
	ia = maph(ih,ic); ! location of amplitudes

		
	n1=eig(itl)%n1 ! dim of final hilbert space
	n2=eig(itl)%n2
	if(.not. master) call allocqt(ia,ic,is,itl,n2)

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

		if (master) then
			it1 = mapt%map(1); ! inition/current hilbert space location
			allocate(Hte(eig(it1)%n2,n2))
			Hte = matmul(eig(it1)%evec,HtUf);
			! Hte is Ht in eigenbasis.
			!	split rhovt into degenerate sectors
			!	and do: Hte^T.rhovtsec.Hte OR Hte.rhovtsec.Hte^T ?
			! for all t-points on coarse grid.
			! get amp2. Multiply penalty function and get rates.
			call merates(ih,ic,is,Hte,eig(it1)%n2,n2,it1,itl)
			deallocate(Hte)
		else
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
		endif
		deallocate(HtUf)  
	endif

	
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



!--------------------------------------------------------
	!	split rhovt into degenerate sectors
	!	and do: Hte^T.rhovtsec.Hte OR Hte.rhovtsec.Hte^T ?
	! for all t-points on coarse grid.
	! get amp2. Multiply penalty function and get rates.
	subroutine merates(ih,ic,is,Hte,n1,n2,it1,it2)
	use modmain, only: rhovt, lrhov,ntcoarse, mrate, eig, maprho,ts
	!use lindblad, only: ? or rhovt in modmain?
	implicit none
	integer, intent(in):: ih,ic,is,n1,n2,it1,it2
	double precision, dimension(n1,n2), intent(in):: Hte ! Ht in eigenbasis
	! local
	double precision:: r
	double precision, allocatable, dimension(:):: amp2
	double precision, allocatable, dimension(:,:):: rho, aux
	integer :: irow,icol,l1,l2,it,j,dl,i
	
	allocate(amp2(eig(it2)%nsec))

	do it=1,ntcoarse
	do i=1,eig(it1)%nsec ! initial spectrum degen sectors
		l1 = eig(it1)%ind(i); l2 = eig(it1)%ind(i+1)-1;
		dl = l2 - l1 + 1;
		!---------------- degen block: rho vector to matrix; 
		allocate(rho(dl,dl)); ! dense block for this degenerate sector
		! fold part of rhov of a degen block to make upper triangular of rho
		j=maprho(l1);
		do irow=1,dl;
			do icol = irow,dl;
				rho(irow,icol) = rhovt(j,it);
				j = j + 1;
			enddo
		enddo
		! fill lower triangular
		do irow=1,dl;
			do icol = 1,irow-1;
				rho(irow,icol) = rho(icol,irow)
			enddo
		enddo
		!--------------------------------
		allocate(aux(dl,n2))
		aux = matmul(rho,Hte(l1:l2,:));
		aux = Hte(l1:l2,:) * aux; ! elem by elem *
		aux(1,:) = sum(aux,dim=1);! amp
		! amp2 for degen sectors:
		call GetAmp2(aux(1,:), n2, amp2, eig(it2)%nsec, eig(it2)%ind)
		deallocate(aux, rho)
		! now get rates using penalty function
		call ratehcst(ih,ic,is,eig(it1)%esec(i),
     .              eig(it2)%nsec,eig(it2)%esec,amp2,r)
		! set global variable mrate:
		! accumulate for all init sectors
		mrate(ih)%rcst(ic,is,it) = mrate(ih)%rcst(ic,is,it)
     .                         + ts(ic,is)*r 
	enddo ! i; initial spec sec
	enddo ! it; time

	deallocate(amp2)
	
	return
	end 	subroutine merates
!--------------------------------------------------------



	subroutine ratehcst(ih,ic,is,Einit,nsec,esec,amp2,r)
	use modmain, only: nx,ways,dqc,dEQs
	implicit none
	integer, intent(in) :: ih,ic,is,nsec
	double precision, intent(in) :: Einit
	double precision, dimension(nsec),intent(in):: esec, amp2
	double precision, intent(out):: r
	
	! local
	double precision, dimension(nsec):: de
	integer :: itl,i
	double precision:: x,y
	
	! if no available hops, set rate = 0
	if (ih <= 8 .and. ways(ih)%ns == 0) then
		r = 0.d0
		return
	endif

	! conditions on nx for ih=1-4
	if(nx==0) then
						! Dhops Homo-Homo channel (ic=1) blocked, 
						! Lumo-homo (ic=3) not available because
						! all spin down: homo filled on all active sites
						! similarly, D,Phi creation only ic=4 possible
						! (conditions in ratehcs()
						! amplitudes for these ic's and ih's not calculated
		if(((ih==1 .or. ih==2) .and. (ic==1 .or. ic==3)) .or. 
     . ((ih==3 .or. ih==4) .and. (ic==2 .or. ic==3))) then
			r = 0.0d0 
			return
		endif
	endif

	! conditions on nx for ih=7,8
	if (ih==7 .or. ih==8) then
		if ( (ic < 3 .and. nx .lt. 1) .or. 
     .   (ic == 3 .and. nx .lt. 2) ) then
			r = 0.d0;
			!write(*,*)"rate: nx cond; ih,ic,is=",ih,ic,is
			return 
		endif
	endif

	! change in energy for all transitions
	de = esec + dqc(ih,ic); ! changes due to efield, barries, etc.
	de = de - Einit; ! total change in energy

	! charging energy contact hops
	if(ih .ge. 9 .and. ih .le. 24) then
		de = de + dEQs(ih-8)
	endif

	r = sum(Penalty(de,nsec) * amp2)
	return
	end subroutine ratehcst
!--------------------------------------------
	function Penalty(de,ne)
	use modmain, only: beta
	implicit none
		integer,intent(in) :: ne
		double precision, dimension(ne),intent(in):: de
		double precision, dimension(ne):: tmp
		double precision, dimension(ne):: Penalty
		! local
		integer:: i
		double precision:: x, fac=30.0d0
		
		Penalty = 1.0d0;
		tmp = de*beta;
		do i=1,ne
			if(tmp(i) > fac) then
				Penalty(i)=0.0d0
			elseif (tmp(i) > 0.0d0) then
					Penalty(i)= dexp(-tmp(i))
			endif
		end do

	return
	end function Penalty
!------------------------------------------
	end 	module amplitudes
