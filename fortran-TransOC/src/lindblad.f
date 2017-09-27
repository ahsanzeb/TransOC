
	module lindblad
	implicit none
	! to use various subroutines in this module
	integer :: ne, nsec
	! make it allocatable or parameter and set using modmain's qt.
	integer, dimension(nsec+1):: ind !start index of sectors

	! rhovt(:,:): rhov for a coarse grained time grid

	! conversion of deg sectors in rhov back to comput basis
	! and calc of rates etc

	! we need to have data for the prev hop so that in case
	! we are trapped, we pretend to go back and solve master eq etc...
	itl = mapt%map(itypes(ih,ic));
	ia = maph(ih,ic); ! location of amplitudes

	!amp, ne,nsec,ind = 
	!	qt(ia)%cs(ic,is)%amp, eig(itl)%n2, eig(itl)%nsec, eig(itl)%ind
  

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
!---------------------------------------
	subroutine mesolve(ntot,dt,ntmax)
	implicit none
	double precision :: dt
	! aux
	double precision, dimension(ntot) :: psit,v
	integer :: i
	double precision :: dth,dt6

	dth = dt/2.0d0;
	dt6 = dt/6.0d0;

	!	calculate rhov, maprho, lrhov
	! calculate gammadec, gammaphi
	! ne,nsec,ind

	do i=2,ntmax
		if (mod(i,prntstep) == 0)
     .   write(6,'(a20,2x,i10,a10,i10)') ' time evolution step ',i,
     .   ' out of ',ntmax
		call rk4step(psit)
	end do

	contains
      !----------------------------------------
      subroutine yprime(vin,vout)
      implicit none
      double precision, intent(in) :: vin(ntot)
      double precision, intent(out) :: vout(ntot)
			call dissipators(ne,nsec,ind,vin,vout)
			return
      end subroutine yprime
      !----------------------------------------
      subroutine rk4step(y)
      implicit none
      double precision, intent(inout) :: y(ntot)
      ! local
			double precision, dimension(ntot) :: k1,k2,k3    
      call yprime(y,k1)
      call yprime(y+k1*dth,k2)
      k1 = k1 + 2*k2
      call yprime(y+k2*dth,k3)
      k1 = k1 + 2*k3
      call yprime(y+k3*dt,k2)
      k1 = k1 + k2
      y = y + k1*dt6
      end subroutine rk4step
      !----------------------------------------
	end subroutine mesolve
!---------------------------------------
	subroutine makerho(amp,ne,nsec,ind)
	use modmain, only: rhov, maprho, lrhov
	implicit none
	! out: sets global arrays wij, rhov, maprho,
	!											 and size of rhov = lrhov
	integer, intent(in) :: ne, nsec
	double precision,dimension(ne),intent(in):: amp
	integer,dimension(nsec+1),intent(in):: ind !start index of sectors
	!integer, intent(out) :: ntot ! size of rhov array

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
	subroutine dissipators(ne,nsec,ind,rhov,Ldiss) ! yprime
	! calculates Ldiss using input rhov
	use modmain, only: maprho, lrhov,
     .               Ldiss, gammadec, gammaphi
	implicit none
	integer, intent(in) :: ne, nsec
	integer,dimension(nsec+1),intent(in):: ind !start index of sectors
	double precision, dimension(lrhov),intent(in):: rhov
	double precision, dimension(lrhov),intent(out):: Ldiss
	
	! local
	!double precision, dimension(lrhov):: Lphi
	double precision, dimension(ne):: Ldecay
	integer:: i,i1,i2,j,k,indx,ntot,g,is,js

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
	subroutine dissipators ! yprime
!---------------------------------------



	end module lindblad
