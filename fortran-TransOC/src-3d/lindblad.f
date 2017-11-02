
	module lindblad
	use modmain, only: dt, ntmax, ntcoarse, beta, wcut, J0,
     .  rhovt,maprho,lrhov,na,commonbath,pauli,ne, nsec
	implicit none
	! to use various subroutines in this module
	!integer :: ne, nsec
	integer, allocatable, dimension(:):: ind !start index of sectors
	double precision, allocatable, dimension(:):: esec, amp
	double precision, allocatable, dimension(:):: rhov ! rho in vectorised form at t=0
	!integer, allocatable, dimension(:):: maprho ! location of diagonal elements of rho in rhov
	!integer:: lrhov
	double precision, allocatable, dimension(:,:):: SD
	double precision, allocatable, dimension(:,:):: Onx !exciton number operator in eigenbasis
	double precision, allocatable, dimension(:,:,:):: Onxi !exciton number operator for individual sites in eigenbasis

	!double precision :: gphi
	! 1/kbT, Ohmic spectral density cut-off and prefix with pi etc absorbed
	!double precision :: beta, wcut, J0
	double precision :: dth,dt6
	!integer, parameter :: ntcoarse = 50 !
	! rhovt(:,:): rhov for a coarse grained time grid

! Lindblad dynamics
! 1. integrate Lindblad equation and find density matrix on a coarse grid
!		solve in eigenbasis: include dephasing and decay, 
!		at t=0: assume zero coherence between states at different energies
!		i- population transfer to lower energy eigenstates
!		ii- decoherence between state in a degenerate sector
!		Use some typical J(w) for environment for dephasing and decay rates
! 2. find transition amplitude squared and rates
!	3. use bisection method to find time t=1/R(t), R(t)= total rate
!	4. make transition at t with rates {rij(t)}
!	..... etc.

	contains
!---------------------------------------------
	function PenaltySqrt(de,ne)
	implicit none
		integer,intent(in) :: ne
		double precision, dimension(ne),intent(in):: de
		double precision, dimension(ne):: tmp
		double precision, dimension(ne):: PenaltySqrt
		! local
		integer:: i
		double precision:: x, fac=40.0d0
		
		PenaltySqrt = 1.0d0;
		tmp = de*beta;
		do i=1,ne
			if(tmp(i) > fac) then
				PenaltySqrt(i)=0.0d0
			elseif (tmp(i) > 0.0d0) then
					PenaltySqrt(i)= dexp(-tmp(i)/2.0d0)
			endif
		end do

	return
	end function PenaltySqrt
!---------------------------------------------
	subroutine initme(ih,ic,is)
	! sets some gloabal variables of this module
	use modmain, only: eig,qt,itypes,mapt,maph,mapc,PermSym,
     .               dqc,dEQs,Einit,Eqtot
	implicit none
	integer, intent(in) :: ih,ic,is
	!local
	integer :: itl,icl,ia,nsecocc,neocc,i,i1,i2
	double precision, allocatable :: de(:)

	! locations of data
	itl = mapt%map(itypes(ih,ic)) ! final hilber space
	icl = mapc(ih,ic); ! channel; reverse of mapc = mapc
	ia = maph(ih,icl); ! amplitudes

	ne = eig(itl)%n2;
	nsec = eig(itl)%nsec;

	allocate(de(nsec))
	de = eig(itl)%esec(1:nsec) + dqc(ih,ic); ! efield, barries, etc.
	de = de - Einit-Eqtot; ! total change
	! charging energy contact hops
	if(ih .ge. 9 .and. ih .le. 24) then
		de = de + dEQs(ih-8)
	endif

	de = PenaltySqrt(de,nsec);

	! leave the high energy sectors with zero populations
	nsecocc=-1;neocc=-1;
	do i=1,nsec
		if(de(i) .lt. 1.0d-10) then ! sqrt=1.d-5 ==> pf=1.d-10
			nsecocc = i-1;
			neocc = eig(itl)%ind(i) - 1;
			exit
		endif
	enddo
	if(nsecocc .ne. -1) then
		! set to occupied sectors only
		nsec = nsecocc
		ne = neocc
	endif
	
	! occupied sectors/states only, pop >= 1.d-10
	if(allocated(esec))deallocate(esec)
	if(allocated(ind))deallocate(ind)
	if(allocated(amp))deallocate(amp)
	allocate(esec(nsec))
	allocate(ind(nsec+1))
	allocate(amp(ne))

	esec = eig(itl)%esec(1:nsec);
	ind = eig(itl)%ind(1:nsec+1);

	if(PermSym) then
		amp = qt(ia)%cs(icl,1)%amp(1:ne)
	else
		amp = qt(ia)%cs(icl,is)%amp(1:ne)
	endif

	! multiply with penalty function sqrt
	do i=1,nsec
		i1 = ind(i); i2=ind(i+1)-1;
		amp(i1:i2) = amp(i1:i2) * de(i);	
	enddo
	! normalise
	amp = amp/sum(amp**2);
	!write(*,*)" norm = ",norm2(amp)
	!write(*,*)" amp = ",amp

	
	return
	end 	subroutine initme	
!----------------------------------------	
	subroutine makerho() !in: amp,ne,nsec,ind
	implicit none
	! out: sets global arrays rhov, maprho, and size of rhov = lrhov
	! local
	integer:: i,i1,i2,j,k,indx,nst
	double precision :: x(ne,ne),y

	! calculate size of vectorised rho with upper triangular only
	lrhov=0;
	do i=1,nsec
		nst = ind(i+1) - ind(i);
		lrhov = lrhov + nst*(nst+1)/2;
	end do

	if(allocated(rhov))deallocate(rhov)
	allocate(rhov(lrhov))

	if(allocated(maprho))deallocate(maprho)
	allocate(maprho(ne))

	! calculate rhov and map for its diagonal elements
	rhov = 0.0d0; indx=0; 
	do i=1,nsec
		i1 = ind(i); i2=ind(i+1)-1;
		do j=i1,i2
			indx = indx + 1;
			! diagonal element
			rhov(indx) = amp(j) * amp(j);
			maprho(j) = indx;
			!off-diagonal elements
			do k=j+1,i2
				indx = indx + 1;
				rhov(indx) = amp(j) * amp(k);
			enddo
		enddo		
	enddo

	return
	end subroutine makerho
!----------------------------------------	
!	diagonal elements of rho only
	subroutine makedrho() !in: amp,ne,nsec,ind
	implicit none
	! out: sets global arrays rhov, and size of rhov = lrhov

	lrhov=ne;	!size
	if(allocated(rhov))deallocate(rhov)
	allocate(rhov(lrhov))
	rhov = amp**2;

	!write(*,*)" lrhov, trace of rho = ", lrhov, sum(rhov)
	!write(*,*)"shape(amp) = ",shape(amp)
	return
	end subroutine makedrho
!-----------------------------------------
	subroutine makeOnx()
	use modmain, only: eig, mapt,na,nx,basis,mapb
	implicit none
	! local
	double precision, allocatable, dimension(:,:):: aux
	integer :: n1,n2,k,i1,i2,itl,ib
	
	itl = mapt%map(1) ! location of eig data
	n1 = eig(itl)%n1; ! dim of hilbert space
	n2 = eig(itl)%n2; ! num of evec we have
	!--------------------------------------
	! calculate Onxc
	allocate(aux(n2,n1))
	aux = transpose(eig(itl)%evec); ! U^T

	ib = mapb%map(3); ! location of basis data

	! U^T.Onxc
	do k=0,min(na,nx);
		! range of basis with exciton number k
		i1=basis(ib)%pntr(k+1)+1;
		i2=basis(ib)%pntr(k+2);
		! relevant terms of U^T.Onxc
		aux(:,i1:i2) = aux(:,i1:i2) * k
		write(*,*)" Onx: k, i1,i2 =",k, i1,i2
	enddo

	if(allocated(Onx))deallocate(Onx)
	allocate(Onx(n2,n2))
	Onx = matmul(aux,eig(itl)%evec) ! (U^T.Onxc).U
	deallocate(aux)


	! debug...
	if (1==1) then
	do k=1,min(10,n2)
		write(*,'(10f10.3)') Onx(k,1:min(10,n2))
	enddo
	endif


	
	!allocate(aux(n2,n2))
	!aux = matmul(transpose(eig(itl)%evec),eig(itl)%evec) ! U^T.U
	open(105,file='UTNU.out',action='write')
	do k=1,n2
		write(105,'('//repeat("E15.3,3x,", n2)//'2x)') Onx(k,:)
	enddo
	close(105)
	!deallocate(aux)

	
	return
	end 	subroutine makeOnx
!---------------------------------------
!---------------------------------------
	subroutine makeOnxi()
	use modmain, only: eig, mapt,na,nx,basis,mapb,isk
	implicit none
	! local
	double precision, allocatable, dimension(:,:,:):: aux
	double precision, allocatable, dimension(:,:):: evect
	integer :: n1,n2,k,i1,i2,itl,ib,ik,i,j
	integer(kind=isk):: l
	
	itl = mapt%map(1) ! location of eig data
	n1 = eig(itl)%n1; ! dim of hilbert space
	n2 = eig(itl)%n2; ! num of evec we have
	!--------------------------------------
	! calculate Onxc

	allocate(evect(n2,n1))
	allocate(aux(n2,n1,na))
	aux = 0.0d0;
	evect = transpose(eig(itl)%evec); ! U^T

	ib = mapb%map(3); ! location of basis data

	!basis(ib)%sec(j)%sets(ntot,j)

	! U^T.Onxc
	do k=1,min(na,nx);
		! range of basis with exciton number k
		i1=basis(ib)%pntr(k+1)+1;
		i2=basis(ib)%pntr(k+2);
		do j=1,i2-i1+1; ! basis in this k-sector
			i = i1 + j - 1 ! absolute index of basis
			do ik=1,k
				! location of jth up spin (in kth sector)
				l = basis(ib)%sec(k)%sets(j,ik);
				! relevant terms of U^T.Onxc
				aux(:,i,l) = evect(:,i);
			enddo
		enddo
	enddo

	if(allocated(Onxi))deallocate(Onxi)
	allocate(Onxi(n2,n2,na))
	do i=1,na
		Onxi(:,:,i) = matmul(aux(:,:,i),eig(itl)%evec) ! (U^T.Onxc).U
	enddo
	deallocate(aux,evect)
	
	return
	end 	subroutine makeOnxi
!---------------------------------------
	double precision function Sw(w)
	double precision, intent(in) :: w
	double precision :: bw
	!Sw = Jw(w) * [ nw(w) + 1 ]; w >= 0
	bw = beta*w;
	if(abs(w) < 1.0d-6) then ! avoid 0/0 
		Sw = J0*dexp(-(w/wcut)**2) *(bw + 1.d0)/(beta*(1.d0 + bw))
		! only w->0 case, beta->0 means T-> inf ==> inf S
	elseif(beta > 0.0d0) then
		Sw = J0 * w * dexp(-(w/wcut)**2) * (1.d0+1.0d0/(dexp(bw)-1.d0))
	else
		Sw	 = 1.0d5 ! a large rate
		write(*,*)"Warning(lindblad): beta=0 ==> S(w) = 10^5 for all w"
	endif
	return
	end function
!---------------------------------------
	subroutine makeSD()
	implicit none
	! local
	integer:: s,r
	double precision:: ws
	
	if(allocated(SD))deallocate(SD)
	allocate(SD(nsec,nsec))
	SD = 0.0d0
	do s=1,nsec
		ws=esec(s)
		do r=1,s
			SD(s,r) = Sw(ws - esec(r))
		enddo
	enddo	

	! debug...
	!sd = J0;
	!write(*,*)"SD = J0 for testing... , J0 = ",J0

	
	return
	end 	subroutine makeSD
!---------------------------------------
	subroutine mesolve()
	implicit none
	if(pauli) then
		call pmesolve()
	else
		call rmesolve()
	endif
	return
	end subroutine mesolve
!---------------------------------------
	subroutine rmesolve()
	!use modmain, only: rhovt	, ntcoarse,J0,wcut
	implicit none

	! aux
	!integer :: ntmax !, ntcoarse
	integer :: i, it1,j,niter1,ii,p
	double precision:: norm
	double precision, allocatable, dimension(:) :: rhovx

	!	calculate rhov, maprho, lrhov
	call makerho()

	! debug...
	allocate(rhovx(lrhov))

	
	! calculate SD matrix, power spectral density; S(w_{sr})
	call makeSD()

	! calculate Onx, exciton number operator in eigenbasis
	if(commonbath) then
		call makeOnx()
	else
		call makeOnxi()
	endif
	
	! set time increments for rk4 routine
	!dt = 0.005; ! ?? choose a suitable value
							!	Ttot = dt*ntmax should be 
							!	large enough > 1/max(th,tl,thl,....)
	dth = dt/2.0d0;
	dt6 = dt/6.0d0;

	if(allocated(rhovt))deallocate(rhovt)
	allocate(rhovt(lrhov,ntcoarse))
	rhovt = 0.0d0;

	!write(*,*)"rhov = ",rhov


	! debug... set off diag to 0
	if (1==0)then
		write(*,*) " testing... diag rho"
		rhovt(:,1) = rhov;	
		rhov = 0.0d0;
		do i=1,ne
			j = maprho(i)
			rhov(j) = rhovt(j,1) 
		enddo
	endif




	niter1 = int(ntmax/ntcoarse)
	do i=1,ntcoarse
		rhovt(:,i) = rhov;	
		call writepop(rhov,lrhov,dt*(i-1)*niter1, 0)
		do j=1,niter1
 			call rk4step(rhov);
 		enddo 		

	! debug...
	if(mod(i-1,100)==0 .and. 1==0) then
		call yprime(rhov,rhovx)
		norm = 0.0d0
		do ii=1,nsec
			j = maprho(ii)
			norm = norm + rhovx(j)
			write(*,*) "it: i, L(i,i) = ",i, ii, rhovx(j)
		enddo
		write(*,*) "it: total change in pop = ",norm
	endif


 		
	end do

	call writepop(rhov,lrhov,0.0d0, 1)

	! debug...
	if (1==0) then
	call yprime(rhovt(:,ntcoarse),rhov)
	norm = 0.0d0
	do i=1,nsec
		j = maprho(i)
		norm = norm + rhov(j)
		write(*,*) "tmax: i, L(i,i) = ",i, rhov(j)
	enddo
	write(*,*) "tmax: total change in pop = ",norm
	endif
	
	return
	end subroutine rmesolve
!---------------------------------------



!---------------------------------------
! Pauli's master equation; just the populations
	subroutine pmesolve()
	implicit none

	! aux
	!integer :: ntmax !, ntcoarse
	integer :: i, it1,j,niter1,ii,p
	double precision:: norm
	double precision, allocatable, dimension(:) :: rhovx

	!	calculate rhov, maprho, lrhov
	call makedrho()

	! calculate SD matrix, power spectral density; S(w_{sr})
	call makeSD()

	! calculate Onx, exciton number operator in eigenbasis
	if(commonbath) then
		call makeOnx()
	else
		call makeOnxi()
	endif
	
	! set time increments for rk4 routine
	!dt = 0.005; ! ?? choose a suitable value
							!	Ttot = dt*ntmax should be 
							!	large enough > 1/max(th,tl,thl,....)
	dth = dt/2.0d0;
	dt6 = dt/6.0d0;

	if(allocated(rhovt))deallocate(rhovt)
	allocate(rhovt(lrhov,ntcoarse))
	rhovt = 0.0d0;

	write(*,*)"shape(rhovt) = ",shape(rhovt)

	niter1 = int(ntmax/ntcoarse)
	do i=1,ntcoarse
		rhovt(:,i) = rhov;	
		call writeppop(rhov,lrhov,dt*(i-1)*niter1, 0)
		do j=1,niter1
 			call rk4step(rhov);
 		enddo 		
	end do

	call writeppop(rhov,lrhov,0.0d0, 1)

	return
	end subroutine pmesolve
!---------------------------------------

	subroutine writepop(rhovx,lrhov,t,gnu)
	!use modmain, only: lrhov
	implicit none
	integer, intent(in):: lrhov, gnu
	double precision, dimension(lrhov), intent(in) :: rhovx
	double precision, intent(in):: t


	integer :: i,ii
	double precision, dimension(nsec):: pop

	!write(*,*) "rhov = ",rhov
	!write(*,*) "lrhov, nsec, ind ",lrhov, nsec , ind
	!write(*,*) "2. maprho", maprho(1:min(100,ne))
	open(120,file='pop.out',action='write',position='append')

	if(gnu==1)then
		write(120,*) 
		write(120,*) 
	else
		pop = 0.0d0;
		do i=1,nsec
			do ii=ind(i),ind(i+1)-1
				pop(i) = pop(i) + rhovx(maprho(ii))
				!write(*,*)" i, pop(i) = ", ii, rhov(maprho(ii))
			enddo
			!sum((/(rhov(maprho(ii)),ii=ind(i),ind(i+1)-1)/))
			!write(*,*)" i, pop(i) = ", i, pop(i)
		enddo
		!write(*,*) "t, pop = ",t, pop
		write(120,*) t, sum(pop), pop
	endif

	close(120)
	
	return
	end 	subroutine writepop
!--------------------------------------	


	subroutine writeppop(rhovx,lrhov,t,gnu)
	!use modmain, only: lrhov
	implicit none
	integer, intent(in):: lrhov, gnu
	double precision, dimension(lrhov), intent(in) :: rhovx
	double precision, intent(in):: t


	integer :: i,ii
	double precision, dimension(nsec):: pop

	!write(*,*) "rhov = ",rhov
	!write(*,*) "lrhov, nsec, ind ",lrhov, nsec , ind
	!write(*,*) "2. maprho", maprho(1:min(100,ne))
	open(120,file='pop.out',action='write',position='append')

	if(gnu==1)then
		write(120,*) 
		write(120,*) 
	else
		pop = 0.0d0;
		do i=1,nsec
			do ii=ind(i),ind(i+1)-1
				pop(i) = sum(rhovx(ind(i):ind(i+1)-1))
			enddo
		enddo
		write(120,*) t, sum(pop), pop
	endif
	close(120)
	
	return
	end 	subroutine writeppop
!--------------------------------------	


	subroutine rk4step(y)
	implicit none
	double precision, intent(inout) :: y(lrhov)
	! local
	double precision, dimension(lrhov) :: k1,k2,k3   
	call yprime(y,k1)
	call yprime(y+k1*dth,k2)
	k1 = k1 + 2*k2
	call yprime(y+k2*dth,k3)
	k1 = k1 + 2*k3
	call yprime(y+k3*dt,k2)
	k1 = k1 + k2
	y = y + k1*dt6
	return
	end subroutine rk4step
!----------------------------------------
	subroutine yprime(rhov,Ldiss)
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	if(pauli) then
		if(commonbath) then
			call pdissipatorcb(rhov,Ldiss)
		else
			call pdissipatorib(rhov,Ldiss)
		endif
	else
		if(commonbath) then
			call dissipatorcb(rhov,Ldiss)
		else
			call dissipatorib(rhov,Ldiss)
		endif
	endif
	return
	end subroutine yprime
!---------------------------------------
! Redfied equation
!	with degeneracies, populations and coherence
!	do not decouple in the secular approximation.
!	For only upper triangular blocks of degen sectors in rho
!	uses exciton number operator Onxs=sum_{i=1,na}(c_i^daggar.c_i)
! Onxs is diagonal in computational basis.
!	In eigenbasis, Onx=U^T.Onxs.U, U eigenstates (columns) in comp basis.
	subroutine dissipatorcb(rhov,Ldiss)
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	integer:: is,js,i1,i2,j1,j2,k1,k2
	integer:: p,q,r,s
	double precision:: rhops, Aps, Nxrs, Swsr
	integer:: i,ir,ld,i0,ld2,j0,lq,ls,i11,j11, indx
	integer, allocatable, dimension(:) :: rcoo1, rcoo2, rcoo3
	
	! before calling this routine:
		! Onx: exciton number operator in eigenbasis
		! declare global in this module and set when mesolve starts
		! Calculate SD(is,js) matrix....

	Ldiss = 0.0d0
	do is=1,nsec
		i1=ind(is); i2=ind(is+1)-1;
		! rcoo: vectorised rho/L coordinates:
		ld = i2-i1 + 1; ! size of this degen block
		i0 = maprho(i1); ! coord of (1,1) elem of this degen block
		allocate(rcoo1(ld)); allocate(rcoo2(ld))
		i11 = i1 - 1;
		do p=i1,i2
			call rowcoor(p-i11,ld,i0,rcoo1); ! coor of rho(p,:)/L(p,:)
			do s=i1,i2
				if(s==p) then
					rcoo2 = rcoo1
				else
					call rowcoor(s-i11,ld,i0,rcoo2); ! coor of rho(s,:)/L(s,:)
				endif
				ls = rcoo1(s-i11); ! location of rho(p,s)
				rhops = rhov(ls)
				Aps = 0.0d0;
				do js=1,is ! js <= is : wsr >=0; no excitations, just decay/dephasing
					j1=ind(js); j2=ind(js+1)-1;
					Swsr = SD(is,js); ! wsr= w(is)-w(js); s in is, r in js sector
					ld2 = j2-j1 + 1;
					j0 = maprho(j1);
					allocate(rcoo3(ld2));
					j11 = j1 - 1;
					do r=j1,j2
						Nxrs = Swsr*Onx(r,s);
						call rowcoor(r-j11,ld2,j0,rcoo3); ! coor of rho(r,:)/L(r,:)
						!-------------------------------------------------------
						! only q=r: Sw(s,r)*N(p,r)*N(r,s)* [ rho(s,:) at (p,:)]
						! How to do this?
						! Aps = Aps + Sw(s,r)*N(p,r)*N(r,s) here, and 
						!	multiply with rho(s,:) after r loop completes
						Aps = Aps + Onx(p,r)*Nxrs;
						!-------------------------------------------------------
						do q=j1,j2; ! lower triangular part of row for h.c
							!	SD(s,r)*N(p,q)*N(r,s)* [ rho(p,s) at (r,q)]
							lq = rcoo3(q-j11); ! location of L(r,q)
							Ldiss(lq) = Ldiss(lq) + Onx(p,q)*Nxrs*rhops
						enddo
					enddo ! r
					deallocate(rcoo3)
				enddo ! js
				!-------------------------------
				! set Ldiss(p,:) = Aps * rho(s,:) here, after sum over r & js
				do i=1,ld; ! lower triangular part of row for h.c
					! rcoo1 = locations of L(p,:)
					! rcoo2 = locations of rho(s,:)
					ir = rcoo1(i);
					Ldiss(ir)	= Ldiss(ir)	- Aps * rhov(rcoo2(i))	
				enddo
			enddo ! is
		enddo ! p
		deallocate(rcoo1,rcoo2)
	enddo ! is

	! h.c: diagonal doubling is missing in above loops
	do p=1,ne
		indx = maprho(p);
		Ldiss(indx) = Ldiss(indx)*2.0d0
	enddo

	!write(*,*)"lindblad: L = ", sum(abs(Ldiss))

	return
	end subroutine dissipatorcb
!---------------------------------------
! for pauli's master eq; diagonal rho case
	subroutine pdissipatorcb(rhov,Ldiss)
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	integer:: is,js,i1,i2,j1,j2,k1,k2
	integer:: p,q,r,s
	double precision:: rhops, Nxrs
	integer:: i,ir,ld,i0,ld2,j0,i11,j11, indx
	integer, allocatable, dimension(:) :: rcoo1, rcoo2, rcoo3
	
	Ldiss = 0.0d0
	do is=1,nsec
		i1=ind(is); i2=ind(is+1)-1;
		do s=i1,i2
				rhops = rhov(s)
				do js=1,is ! js <= is : wsr >=0; no excitations, just decay/dephasing
					j1=ind(js); j2=ind(js+1)-1;
					! wsr= w(is)-w(js); s in is, r in js sector
					rhops = rhops * SD(is,js);
					do r=j1,j2
						Nxrs = rhops * Onx(r,s)**2
						Ldiss(r) = Ldiss(r) + Nxrs
						Ldiss(s)	= Ldiss(s)	- Nxrs
					enddo ! r
				enddo ! js
			enddo ! is
	enddo ! is

	! h.c: diagonal doubling is missing in above loops
	Ldiss = Ldiss*2.0d0

	return
	end subroutine pdissipatorcb
!---------------------------------------
! Redfied equation
!	with degeneracies, populations and coherence
!	do not decouple in the secular approximation.
!	For only upper triangular blocks of degen sectors in rho
!	uses exciton number operator Onxsi=c_i^daggar.c_i
! All Onxi's are diagonal in computational basis.
!	In eigenbasis, Onxi=U^T.Onxsi.U, U eigenstates (columns) in comp basis.
	subroutine dissipatorib(rhov,Ldiss)
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	integer:: is,js,i1,i2,j1,j2,k1,k2
	integer:: p,q,r,s
	double precision:: rhops, Aps, Swsr
	integer:: i,ir,ld,i0,ld2,j0,lq,ls,i11,j11, indx
	integer, allocatable, dimension(:) :: rcoo1, rcoo2, rcoo3
	double precision, dimension(na) :: Nxrs

	!write(*,*)"Individual baths: shape(Nxrs) = ",shape(Nxrs)
	!write(*,*)"shape(Onxi) = ",shape(Onxi)

	
	Ldiss = 0.0d0
	do is=1,nsec
		i1=ind(is); i2=ind(is+1)-1;
		! rcoo: vectorised rho/L coordinates:
		ld = i2-i1 + 1; ! size of this degen block
		i0 = maprho(i1); ! coord of (1,1) elem of this degen block
		allocate(rcoo1(ld)); allocate(rcoo2(ld))
		i11 = i1 - 1;
		do p=i1,i2
			call rowcoor(p-i11,ld,i0,rcoo1); ! coor of rho(p,:)/L(p,:)
			do s=i1,i2
				if(s==p) then
					rcoo2 = rcoo1
				else
					call rowcoor(s-i11,ld,i0,rcoo2); ! coor of rho(s,:)/L(s,:)
				endif
				ls = rcoo1(s-i11); ! location of rho(p,s)
				rhops = rhov(ls)
				Aps = 0.0d0;
				do js=1,is ! js <= is : wsr >=0; no excitations, just decay/dephasing
					j1=ind(js); j2=ind(js+1)-1;
					Swsr = SD(is,js); ! wsr= w(is)-w(js); s in is, r in js sector
					ld2 = j2-j1 + 1;
					j0 = maprho(j1);
					allocate(rcoo3(ld2));
					j11 = j1 - 1;
					do r=j1,j2
						Nxrs = Swsr*Onxi(r,s,:);
						call rowcoor(r-j11,ld2,j0,rcoo3); ! coor of rho(r,:)/L(r,:)
						!-------------------------------------------------------
						! only q=r: Sw(s,r)*N(p,r)*N(r,s)* [ rho(s,:) at (p,:)]
						! How to do this?
						! Aps = Aps + Sw(s,r)*N(p,r)*N(r,s) here, and 
						!	multiply with rho(s,:) after r loop completes
						Aps = Aps + sum(Onxi(p,r,:)*Nxrs); ! elem by elem multiply
						!-------------------------------------------------------
						do q=j1,j2; ! lower triangular part of row for h.c
							!	SD(s,r)*N(p,q)*N(r,s)* [ rho(p,s) at (r,q)]
							lq = rcoo3(q-j11); ! location of L(r,q)
							Ldiss(lq) = Ldiss(lq) + sum(Onxi(p,q,:)*Nxrs)*rhops
						enddo
					enddo ! r
					deallocate(rcoo3)
				enddo ! js
				!-------------------------------
				! set Ldiss(p,:) = Aps * rho(s,:) here, after sum over r & js
				do i=1,ld; ! lower triangular part of row for h.c
					! rcoo1 = locations of L(p,:)
					! rcoo2 = locations of rho(s,:)
					ir = rcoo1(i);
					Ldiss(ir)	= Ldiss(ir)	- Aps * rhov(rcoo2(i))	
				enddo
			enddo ! is
		enddo ! p
		deallocate(rcoo1,rcoo2)
	enddo ! is

	! h.c: diagonal doubling is missing in above loops
	do p=1,ne
		indx = maprho(p);
		Ldiss(indx) = Ldiss(indx)*2.0d0
	enddo

	!write(*,*)"lindblad: L = ", sum(abs(Ldiss))

	return
	end subroutine dissipatorib
!---------------------------------------
! for pauli's master eq; diagonal rho case
	subroutine pdissipatorib(rhov,Ldiss)
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	integer:: is,js,i1,i2,j1,j2,k1,k2
	integer:: p,q,r,s
	double precision:: rhops, Nxrs
	integer:: i,ir,ld,i0,ld2,j0,i11,j11, indx
	integer, allocatable, dimension(:) :: rcoo1, rcoo2, rcoo3
	
	Ldiss = 0.0d0
	do is=1,nsec
		i1=ind(is); i2=ind(is+1)-1;
		do s=i1,i2
				rhops = rhov(s)
				do js=1,is ! js <= is : wsr >=0; no excitations, just decay/dephasing
					j1=ind(js); j2=ind(js+1)-1;
					! wsr= w(is)-w(js); s in is, r in js sector
					rhops = rhops * SD(is,js);
					do r=j1,j2
						Nxrs = rhops * sum(Onxi(r,s,:)**2)
						Ldiss(r) = Ldiss(r) + Nxrs
						Ldiss(s)	= Ldiss(s)	- Nxrs
					enddo ! r
				enddo ! js
			enddo ! is
	enddo ! is

	! h.c: diagonal doubling is missing in above loops
	Ldiss = Ldiss*2.0d0

	return
	end subroutine pdissipatorib
!--------------------------------------------
!	finds coordinates of rows: rho(s,:)
	! since rho and L have only upper triangular of degen blocks:
	! 	bend any row of rho or L to 90deg up and use columns instead.
	subroutine rowcoor(s,ld,i0,rcoo)
	implicit none
	integer, intent(in) :: s ! row number that needs coor in vectzd format
	integer, intent(in) :: ld ! size of the block
	integer, intent(in) :: i0 ! i0 = coordinate of block's (1,1) element
	integer, dimension(ld),intent(out) :: rcoo ! output coordinates
	! local
	integer :: i,ld1

	!write(*,*)"rowcoor: s,ld,i0 = ", s,ld,i0

	! part of sth row bended to sth col
	rcoo(1) = (i0-1) + s;
	ld1 = ld+1
	do i=2,s
		rcoo(i) = rcoo(i-1) + ld1 - i;
	enddo
	! the rest of sth row, not bended
	do i=s+1,ld
		rcoo(i) = rcoo(i-1) + 1
	enddo

	return
	end 	subroutine rowcoor
!---------------------------------------
	end module lindblad
