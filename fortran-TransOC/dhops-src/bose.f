
! calculates bosonic states of microcavity. That is, any number of LP,UP,Dark states...
!	from the vacuum.

	module bose
	use basisstates, only: nCr
	use hamiltonian, only: HgMap, getDownSpins
	use modmain, only: 

	implicit none



	! to set kind of integers in ksub
	! kind 1 (8bits) can go upto 127, which should be enough 
	! for us if nsites <= 127.
	!	kind 2 (16 bits) max is 32,767
	! selected_int_kind(R) ==> smallest kind, to 10^R exclusive 
	integer, parameter :: isk=selected_int_kind(2);

	
	type :: statesec
		integer :: ntot,nx,np
		integer(kind=isk), allocatable :: dn(:,:) ! sets of down spins (:,1:n-k) ! same order as in map(:,1:n-k)
		double complex, allocatable :: c(:) ! coefficient of a state for this sector
		integer, allocatable :: map(:,:) ! dim = ntot, n-nx
	end type statesec

	type :: bostate
		integer :: ntot
		!logical :: xst
		integer :: n
		integer :: maxk ! max k of k-sub it has.
		integer, allocatable :: pntr(:)
		type(statesec), allocatable :: sec(:)
	end type bostate

	type(statesec), allocatable :: fmap(:) ! full map for all sectors
	type(bostate) :: vac, vec1, vec2




	contains

	! routines... 
!......................................................................
	subroutine getfullmap(n)
	implicit none
	integer, intent(in) :: n
	! local
	integer :: k,ntot, ibl

	ibl = 1 ! basis for N sites, determined by mapb%map
	
	allocate(fmap(0:n-1))
	do k=0,n-1
		ntot = 	nCr(n,k)
		fmap(k)%ntot = ntot
		allocate(fmap(k)%map(ntot,n-k))
		call HgMap(ibl,n,k,ntot,fmap(k)%map)
		allocate(fmap(k)%dn(ntot,n-k))
		call getDownSpins(ibl,n,k,ntot,fmap(k)%dn)
	end do
	
	return
	end 	subroutine getfullmap
!......................................................................
	subroutine addLP(n,nex,alpha,beta)
	implicit none
	integer, intent(in) :: n,nex
	double precision, intent(in):: alpha, beta
	integer :: k,i,j,l

	! restrict k loop over existing k-sectors, dont processes zero-coeff sectors

	! alpha * photon creation 
	do k=0,n ! might be extending non-zero coeff sectors, let it go to max n, not costly.
		v2%sec(k)%np = v2%sec(k)%np + 1;
		v2%sec(k)%c(:) = v2%sec(k)%c(:) + alpha*v1%sec(k)%c(:)
	end do

	! beta * bright exciton creation
	do k=0,min(nex,n-1)
		do i=1,ntot ! basis state of k-th sector
			do l=1,n-k
				j = fmap(k)%map(i,l)
				v2%sec(k+1)%c(j) = v2%sec(k+1)%c(j) + beta*v1%sec(k)%c(i)
			end do
		end do
	end do

	return
	end 	subroutine addLP
!......................................................................
	subroutine addDark(n,nex,kk)
	implicit none
	integer, intent(in) :: n,nex,kk
	integer :: k,i,j,l,is
	double complex :: tpn, phase

	tpn = (0.0d0, 2*3.14159265359d0/n); ! iota * 2 *pi/n

	! restrict k loop over existing k-sectors, dont processes zero-coeff sectors
	do k=0,min(nex,n-1)
		kktpn = kk * tpn
		! dark exciton creation with index kk
		do i=1,ntot ! basis state of k-th sector
			do l=1,n-k
				j = fmap(k)%map(i,l)
				is = fmap(k)%dn(i,l)
				phase = dexp(is * kktpn)
				v2%sec(k+1)%c(j) = v2%sec(k+1)%c(j) + 
     .                     beta * v1%sec(k)%c(i) * phase
			end do
	end do

	return
	end 	subroutine addDark
!......................................................................






	

	end 	module bose
