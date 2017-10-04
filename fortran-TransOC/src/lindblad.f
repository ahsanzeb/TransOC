
	module lindblad
	use modmain, only: dt, ntmax, ntcoarse, beta, wcut, J0,
     .              rhovt,maprho,lrhov
	implicit none
	! to use various subroutines in this module
	integer :: ne, nsec
	integer, allocatable, dimension(:):: ind !start index of sectors
	double precision, allocatable, dimension(:):: esec, amp
	double precision, allocatable, dimension(:):: rhov ! rho in vectorised form at t=0
	!integer, allocatable, dimension(:):: maprho ! location of diagonal elements of rho in rhov
	!integer:: lrhov
	double precision, allocatable, dimension(:,:):: SD
	double precision, allocatable, dimension(:,:):: Onx !exciton number operator in eigenbasis
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
	subroutine initme(ih,ic,is)
	! sets gloabal variables of this module
	use modmain, only: eig,qt,itypes,mapt,maph,mapc,PermSym
	implicit none
	integer, intent(in) :: ih,ic,is
	!local
	integer :: itl,icl,ia

	! locations of data
	itl = mapt%map(itypes(ih,ic)) ! final hilber space
	icl = mapc(ih,ic); ! channel; reverse of mapc = mapc
	ia = maph(ih,icl); ! amplitudes

	ne = eig(itl)%n2;
	nsec = eig(itl)%nsec;

	if(allocated(esec))deallocate(esec)
	if(allocated(ind))deallocate(ind)
	if(allocated(amp))deallocate(amp)
	allocate(esec(nsec))
	allocate(ind(nsec+1))
	allocate(amp(ne))
	
	esec = eig(itl)%esec(1:nsec);
	ind = eig(itl)%ind;

	write(*,*)"lindblad: esec(1:5) = ", esec(1:min(nsec,5))
	write(*,*)"lindblad: ind(1:5) = ", ind(1:min(nsec,5))
		

	
	if(PermSym) then
		amp = qt(ia)%cs(icl,1)%amp
	else
		amp = qt(ia)%cs(icl,is)%amp
	endif

	amp = amp/norm2(amp);
	!if(dabs(sum(amp**2)-1.0d0) > 1.0d-3) then
	!	write(*,*)"initme: sum(amp**2) < 1 ???"
	!endif
	write(*,*)"lindblad: amp(1:5) = ", amp(1:min(ne,20))

	
	return
	end 	subroutine initme	
!---------------------------------------
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

	ib = mapb%map(3) ! location of basis data

	! U^T.Onxc
	do k=0,min(na,nx);
		! range of basis with exciton number k
		i1=basis(ib)%pntr(k+1)+1;
		i2=basis(ib)%pntr(k+2);
		! relevant terms of U^T.Onxc
		aux(:,i1:i2) = aux(:,i1:i2) * k
	enddo

	if(allocated(Onx))deallocate(Onx)
	allocate(Onx(n2,n2))
	Onx = matmul(aux,eig(itl)%evec) ! (U^T.Onxc).U
	deallocate(aux)
	
	return
	end 	subroutine makeOnx
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
	integer:: is,js
	double precision:: wi
	
	if(allocated(SD))deallocate(SD)
	allocate(SD(nsec,nsec))
	SD = 0.0d0
	do is=1,nsec
		wi=esec(is)
		do js=is,nsec
			SD(js,is) = Sw(esec(js) - wi)
		enddo
	enddo	

	! debug...
	sd = 1.0d0

	
	return
	end 	subroutine makeSD
!---------------------------------------
	subroutine mesolve()
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
	call makeOnx()

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
	rhovt(:,1) = rhov;	
	rhov = 0.0d0;
	do i=1,nsec
		j = maprho(i)
		rhov(j) = rhovt(j,1) 
	enddo




	niter1 = int(ntmax/ntcoarse)
	do i=1,ntcoarse
		rhovt(:,i) = rhov;		
		call writepop(rhov,lrhov,dt*(i-1)*niter1, 0)
		
		do j=1,niter1
 			call rk4step(rhov);
 		enddo 		

	if(mod(i-1,10)==0 .and. i < 100) then
		! debug...
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
	call yprime(rhovt(:,ntcoarse),rhov)
	norm = 0.0d0
	do i=1,nsec
		j = maprho(i)
		norm = norm + rhov(j)
		write(*,*) "tmax: i, L(i,i) = ",i, rhov(j)
	enddo
	write(*,*) "tmax: total change in pop = ",norm

	
	return
	end subroutine mesolve
!---------------------------------------

	subroutine writepop(rhov,lrhov,t,gnu)
	!use modmain, only: lrhov
	implicit none
	integer, intent(in):: lrhov, gnu
	double precision, dimension(lrhov), intent(in) :: rhov
	double precision, intent(in):: t


	integer :: i,ii
	double precision, dimension(nsec):: pop


	!write(*,*) "rhov = ",rhov
	!write(*,*) "lrhov, t, gnu, nsec ",lrhov, t, gnu, nsec 
	!write(*,*) "maprho", maprho(1:min(100,ne))
	open(120,file='pop.out',action='write',position='append')

	if(gnu==1)then
		write(120,*) 
		write(120,*) 
	else
		pop = 0.0d0;
		do i=1,nsec
			do ii=ind(i),ind(i+1)-1
				pop(i) = pop(i) + rhov(maprho(ii))
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
	subroutine makerho() !in: amp,ne,nsec,ind
	implicit none
	! out: sets global arrays rhov, maprho, and size of rhov = lrhov
	! local
	integer:: i,i1,i2,j,k,indx,nst
	double precision :: x(ne,ne)

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

	!write(*,*)"norm |amp><amp| = ",
  !   .matmul(reshape(amp,shape=(/1,ne/)),reshape(amp,shape=(/ne,1/)))
	
	x=matmul(reshape(amp,shape=(/ne,1/)),reshape(amp,shape=(/1,ne/)))
	!write(*,*)"makerho: res rho-rho^2 = ",
  !   . sum(abs(x - matmul(x,x)))
	write(*,*)"============================="
	do i=1,ne
		write(*,*)i, x(i,:)
	enddo
	write(*,*)"============================="
	
	return
	end subroutine makerho

!---------------------------------------
! Redfied equation
!	with degeneracies, populations and coherence
!	do not decouple in the secular approximation.
!	For only upper triangular blocks of degen sectors in rho
!	uses exciton number operator Onxs=sum_{i=1,na}(c_i^daggar.c_i)
! Onxs is diagonal in computational basis.
!	In eigenbasis, Onx=U^T.Onxs.U, U eigenstates (columns) in comp basis.
	subroutine yprime(rhov,Ldiss)
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	integer:: is,js,i1,i2,j1,j2,k1,k2
	integer:: p,q,r,s
	integer:: indx,indx2,indx3
	double precision:: rhops, Aps, Nxrs, Swsr
	integer:: i,ir,ld,i0
	integer, allocatable, dimension(:) :: rcoo1, rcoo2
	
	! before calling this routine:
		! Onx: exciton number operator in eigenbasis
		! declare global in this module and set when mesolve starts
		! Calculate SD(is,js) matrix....

	Ldiss = 0.0d0
	do is=1,nsec
	i1=ind(is); i2=ind(is+1)-1;
	do p=i1,i2
	indx = maprho(p); ! location of rho_pp
	! rcoo: vectorised rho/L coordinates:
	ld = i2-i1 + 1; ! size of this degen block
	i0 = maprho(i1); ! coord of (1,1) elem of this degen block
	allocate(rcoo1(ld)); allocate(rcoo2(ld))
	call rowcoor(p-i1+1,ld,i0,rcoo1); ! coor of L(p,:)
	!=============== s >= p ===================
	do s=p,i2 ! s >= p
		rhops=rhov(indx); ! use: indx+s-p OR indx += 1
		indx = indx + 1;
		call rowcoor(s-i1+1,ld,i0,rcoo2); ! coor of rho(s,:)
		Aps = 0.0d0;
		do js=1,is ! js <= is
		j1=ind(js); j2=ind(js+1)-1;
		Swsr = SD(is,js); ! wsr= w(is)-w(js); s in is, r in js sector
		do r=j1,j2
		indx2 = maprho(r); ! location of rho_rr
		Nxrs = Swsr*Onx(r,s);
		!-------------------------------------------------------
		! only q=r: Sw(s,r)*N(p,r)*N(r,s)* [ rho(s,:) at (p,:)]
		! How to do this?
		! Aps = Aps + Sw(s,r)*N(p,r)*N(r,s) here, and 
		!	multiply with rho(s,:) after r loop completes
		Aps = Aps + Onx(p,r)*Nxrs;
		!-------------------------------------------------------
		do q=r,j2; ! q >= r
			!	SD(s,r)*N(p,q)*N(r,s)* [ rho(p,s) at (r,q)]
			Ldiss(indx2) = Ldiss(indx2) + Onx(p,q)*Nxrs*rhops
			indx2 = indx2 + 1;
		enddo
		enddo
		enddo
		!-------------------------------
		! set Ldiss(p,:) = Aps * rho(s,:) here, after sum over r; all jsec
		! since rho and L have only upper triangular of degen blocks:
		! 		bend the s and p rows of rho and L to 90deg up
		!		and use columns instead.
		do i=1,ld
			ir = rcoo1(i);
			Ldiss(ir)	= Ldiss(ir)	- Aps * rhov(rcoo2(i))	
		enddo
		! also do s<p case. before or after this s loop, with inner r loop
		!-------------------------------
	enddo ! is >= p
	!=============== s < p ===================
	do s=i1,p-1
		call rowcoor(s-i1+1,ld,i0,rcoo2); ! coor of rho(s,:)
		Aps = 0.0d0;
		do js=1,is ! js <= is
			Swsr = SD(is,js); ! wsr= w(is)-w(js); s in is, r in js sector
			do r=ind(js),ind(js+1)-1
				Aps = Aps + Swsr*Onx(p,r)*Onx(r,s)
			enddo
		enddo
		do i=1,ld
			ir = rcoo1(i);
			Ldiss(ir)	= Ldiss(ir)	- Aps * rhov(rcoo2(i))	
		enddo	
	enddo
	!==================================
	deallocate(rcoo1,rcoo2)	
	enddo ! p
	enddo

	! + h.c ===> diagonal = diagonal*2
	do p=1,ne
		indx = maprho(p);
		Ldiss(indx) = Ldiss(indx)*2.0d0
	enddo

	!write(*,*)"lindblad: L = ", sum(abs(Ldiss))



	return
	end subroutine yprime
!---------------------------------------
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
