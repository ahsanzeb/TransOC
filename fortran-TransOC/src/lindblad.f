
	module lindblad
	implicit none
	! to use various subroutines in this module
	integer :: ne, nsec
	integer, allocatable, dimension(:):: ind !start index of sectors
	double precision, allocatable, dimension(:):: esec, amp
	double precision, allocatable, dimension(:):: rhov
	integer, allocatable, dimension(:):: maprho ! location of diagonal elements of rho in rhov
	integer:: lrhov


	
	! rhovt(:,:): rhov for a coarse grained time grid

	! conversion of deg sectors in rhov back to comput basis
	! and calc of rates etc

	! we need to have data for the prev hop so that in case
	! we are trapped, we pretend to go back and solve master eq etc...


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
	use modmain, only: eig,qt,itypes,mapt,maph,mapc
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
	amp = qt(ia)%cs(icl,is)%amp

	return
	end 	subroutine initme	
!---------------------------------------
	subroutine mesolve(dt,ntmax)
	implicit none
	double precision :: dt
	! aux
	double precision, dimension(lrhov) :: psit,v
	integer :: i
	double precision :: dth,dt6

	dth = dt/2.0d0;
	dt6 = dt/6.0d0;

	!	calculate rhov, maprho, lrhov
	call makerho()

	! calculate gammadec, gammaphi




	do i=2,ntmax
		if (mod(i,prntstep) == 0)
     .   write(6,'(a20,2x,i10,a10,i10)') ' time evolution step ',i,
     .   ' out of ',ntmax
		call rk4step(psit)
	end do

	end subroutine mesolve
!---------------------------------------
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
	integer:: i,i1,i2,j,k,indx

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
	subroutine makerho
!---------------------------------------
	subroutine yprime(rhov,Ldiss) ! dissipators
	! calculates Ldiss using input rhov
	use modmain, only: gammadec, gammaphi ! define these in this module?
	implicit none
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	double precision, dimension(ne):: Ldecay
	integer:: i,i1,i2,j,k,indx,lrhov,g,is,js

	! NOTE:
	! calculate gammadec, gammaphi before calling this.
	
	! decay; diagonal elements only
	Ldecay = 0.0d0;
	do js=1,nsec
		do j=ind(js),ind(js+1)-1 ! eigenstates in js sector

			indx = maprho(j); ! location of jth diagonal element of rho in rhov

			! lower sectors reduce populations
			do is=1,js-1
				g = ind(is+1) - ind(is); ! degeneracy of ith sector
				Ldecay(j) = Ldecay(j) - gammadec(js,is)*g*rhov(indx);
			enddo

			! higher sectors increase populations
			do is=js+1,nsec
				do i=ind(is),ind(is+1)-1
					indx = maprho(i);
					Ldecay(j) = Ldecay(j) + gammadec(is,js)*rhov(indx);
				enddo
			enddo

		enddo !j
	enddo ! js

	! dephasing ==> off-diagonal of rho
	!	Lphi = gammaphi * rhov;
	! fill diagonal elements with Ldecay
	Ldiss = gammaphi * rhov; 
	do i=1,ne
		indx = maprho(i);
		Ldiss(indx) = Ldecay(i);
	enddo

	return
	subroutine yprime
!---------------------------------------



	end module lindblad
