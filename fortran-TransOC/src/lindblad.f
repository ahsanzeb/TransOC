
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
	double precision, allocatable, dimension(:,:):: gdecay
	double precision :: gphi
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
	if(PermSym) then
		amp = qt(ia)%cs(icl,1)%amp
	else
		amp = qt(ia)%cs(icl,is)%amp
	endif

	amp = amp/norm2(amp);
	!if(dabs(sum(amp**2)-1.0d0) > 1.0d-3) then
	!	write(*,*)"initme: sum(amp**2) < 1 ???"
	!endif
	
	return
	end 	subroutine initme	
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
	subroutine makegammas()
	implicit none
	! local
	integer:: is,js
	double precision:: w
	
	if(allocated(gdecay))deallocate(gdecay)
	allocate(gdecay(nsec,nsec))
	gdecay = 0.0d0
	! decay rates
	do js=1,nsec
		do is=1,js-1
			w = esec(js) - esec(is);
			gdecay(js,is) = 2*Sw(w)
		enddo
	enddo	
	gphi = 2*Sw(0.0d0); ! dephasing rate

	!write(*,*)"gphi, gdecay="
	!write(*,*)gphi
	!write(*,*)gdecay
	
	return
	end 	subroutine makegammas
!---------------------------------------
	subroutine mesolve()
	!use modmain, only: rhovt	, ntcoarse,J0,wcut
	implicit none

	! aux
	!integer :: ntmax !, ntcoarse
	integer :: i, it1,j,niter1,ii
	double precision:: norm

	!	calculate rhov, maprho, lrhov
	call makerho()
	
	! calculate gdecay, gphi
	call makegammas()

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

	niter1 = int(ntmax/ntcoarse)
	do i=1,ntcoarse
		rhovt(:,i) = rhov;		
		call writepop(rhov,lrhov,dt*(i-1)*niter1+1, 0)
		
		do j=1,niter1
 			call rk4step(rhov);
 		enddo
	end do

	call writepop(rhov,lrhov,0.0d0, 1)

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
		write(120,*) t, pop
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
!---------------------------------------
	subroutine yprime(rhov,Ldiss) ! dissipators
	! calculates Ldiss using input rhov
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	double precision, dimension(ne):: Ldecay
	integer:: i,i1,i2,j,k,indx,lrhov,g,is,js

	! NOTE:
	! calculate gdecay, gphi before calling this.
	
	! decay; diagonal elements only
	Ldecay = 0.0d0;
	do js=1,nsec
		do j=ind(js),ind(js+1)-1 ! eigenstates in js sector

			indx = maprho(j); ! location of jth diagonal element of rho in rhov

			! lower sectors reduce populations
			do is=1,js-1
				g = ind(is+1) - ind(is); ! degeneracy of ith sector
				Ldecay(j) = Ldecay(j) - gdecay(js,is)*g*rhov(indx);
			enddo

			! higher sectors increase populations
			do is=js+1,nsec
				do i=ind(is),ind(is+1)-1
					indx = maprho(i);
					Ldecay(j) = Ldecay(j) + gdecay(is,js)*rhov(indx);
				enddo
			enddo

		enddo !j
	enddo ! js

	! dephasing ==> off-diagonal of rho
	!	Lphi = gphi * rhov;
	! fill diagonal elements with Ldecay
	Ldiss = gphi * rhov; 
	do i=1,ne
		indx = maprho(i);
		Ldiss(indx) = Ldecay(i);
	enddo

	if((.not. isnan(sum(rhov))) .and. isnan(sum(Ldiss))) then
		write(*,*)"Ldiss=",Ldiss
	endif

	return
	end subroutine yprime
!---------------------------------------

	end module lindblad
