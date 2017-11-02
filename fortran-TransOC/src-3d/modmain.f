
! Copyright (C) 2017 M. Ahsan Zeb
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

	module modmain
	implicit none

	! to set kind of integers in ksub
	! kind 1 (8bits) can go upto 127, which should be enough 
	! for us if nsites <= 127.
	!	kind 2 (16 bits) max is 32,767
	! selected_int_kind(R) ==> smallest kind, to 10^R exclusive 
	integer, parameter :: isk=selected_int_kind(2);

	logical :: debug

	! total number of sites in the system
	integer:: nsites

	! number of electrons 
	! DopingRange
	integer :: nel,nelmin,nelmax,dnelec ! nel will stay fixed if no contacts

	! net charge on the system
	integer :: Qnet
	! Charging energy prop constant
	double precision:: Eq
	! charging energy for contact hops
	double precision:: dEQs(16), Eqtot
	
	! undoped case can get into traps,
	! and until no clear release mechanism is decided,
	!	e.g., full dissipative dynamics between hops,
	logical :: onlydoped

	! excitations:
	integer :: mexmin,mexmax,dmex

	! DetuningRange
	double precision :: dwmin, dwmax, ddw
	integer :: ndw
	! active sites
	integer:: na
	! number of excitatons
	integer:: nx

	! total number of iterations/hops in the trajectory
	integer:: niter

	! number of trajectories
	integer :: ntraj
	
	! L-H and H-L cross hops allowed? 
	logical :: crosshops 

	! output options:
	logical::  ztout, ntrapout, ratesout

	! no matter-light coupling?
	logical :: nog
	integer :: ipsi ! to store index of state when nog=T


	! minimum D, Phi carriers at start of a trajectory?
	logical :: mincarriers
	
	! all molecules-light couplings same?
	logical :: sameg
	
	! reuse basis/hamiltonians (fixmap=F) or not (fixmap=T)
	logical :: fixmap
	! list of active sites in order they appear in the basis 
	integer, allocatable :: ASites(:)

	! PermSym: if all sites have same w0 and g, and AlwaysLP=True, i.e,
	!	 	a quick relaxation to Lower polariton state after every hop,
	!		then we are always in the LP state when any type of hopping
	!		occurs, all sites will be similar and we can use amplitudes
	!		calculated for one site for others.
	logical:: PermSym

	!	na, nx lists for 13 Hilbert spaces
	integer, dimension(13):: nalist,nxlist !CHECK  not used????

	! index for jump, channel, site
	integer:: itype

	! hopping parameters for Homo-Homo, etc
	!	set defauls after reading input file
	double precision, dimension(8):: tpar
	
	! cavity and exciton loss rates
	double precision :: kappa, gamma	
	logical :: nokappa, nogamma
	! hopping parameters
	double precision, dimension(26,4):: ts
	! energy changes due to contact barriers, applied field, etc
	double precision, dimension(26,4):: dqc
	! hopping parameters
	double precision:: th, tl, tlh, thl
	double precision:: JhR, JlR, JhL, JlL
	! bare exciton energy, right and left contact barriers
	!	Electric field energy Er assumed constant for all hops
	double precision:: w0, Ebr,Ebl, Er
	! beta = 1/KbT for penalty function
	double precision:: beta	
	! block injection of electron/holes?
	!logical, dimension(2):: BlockInjection
	logical:: EBlock, HBlock
	!	periodic boundary conditions?
	logical:: periodic ! XXXXXXX make it about geometry of the system
	! only bulk processes? with or without periodic boundary conditions
	logical :: onlybulk ! no contacts?
	logical :: nolosses ! kappa=0=gamma ?
	logical :: AlwaysLP ! quickly relax to Lower Polariton after a hop
	double precision:: Einit! intial energy in every iteration
	double precision, allocatable:: Einit2(:)
	!	--------------- maps ---------------
	integer, dimension(13)::
     . dna	= (/ 0,0,0,2,2,2,-2,-2,-2,1,1,-1,-1 /);
	integer, dimension(13)::
     . dnx = (/ 0,-1,1,1,0,2,-1,-2,0,1,0,-1,0 /);
	integer, dimension(13)::
     . ibs = (/ 3,3,3,5,5,5,1,1,1,4,4,2,2 /);
	integer, dimension(5)::
     . dns = (/ -2,-1,0,1,2 /); 
	integer, dimension(26,4)::
     . itypes = reshape( ( / 1, 1, 2, 3, 1, 1, 2, 3, 1, 1,
     . 2, 3, 1, 1, 2, 3, 4, 4, 5, 6,
     . 4, 4, 5, 6, 7, 7, 8, 9, 7, 7,
     . 8, 9, 11, 11, 11, 11, 10, 10, 
     . 10, 10, 11, 11, 11, 11, 10, 10,
     . 10, 10, 11, 11, 11, 11, 10,
     . 10, 10, 10, 11, 11, 11, 11,
     . 10, 10, 10, 10, 13, 13, 13,
     . 13, 13, 13, 13, 13, 12, 12,
     . 12, 12, 12, 12, 12, 12, 13,
     . 13, 13, 13, 13, 13, 13, 13,
     . 12, 12, 12, 12, 12, 12, 12,
     . 12, 2, 2, 2, 2, 2, 2, 2, 2 /),
     . (/ 26,4 /), order=(/2,1/) ); 


	type :: BMaps
		integer :: nnu
		logical, dimension(5):: req ! required in an iteration or not?
		integer, dimension(5) :: map	 ! map for ib
		integer, dimension(5) :: cal ! which ib
	end type BMaps

	type :: TMaps
		integer :: nnu
		integer, dimension(13) :: map ! map for ib
		integer, dimension(13) :: cal ! which ib
		! 13 types: which ones are required 
		! based on availability of relevant hops
		logical, dimension(13):: req ! ReqType
		! which itypes for a given ib; ntb,grouptb only for mapt
		integer, dimension(5) :: ntb 
		integer, dimension(5,13) :: grouptb 
	end type TMaps

	type(BMaps) :: mapb
	type(TMaps) :: mapt
	!	------------------------------


	!---------------------------------------	
	! for 5 values of N; smaller m can use larger m's data
	! basis sectors for a given na,k-up spins
	!---------------------------------------
	type :: BSectors
		integer :: ntot ! will see later if these ntot are needed?
		integer(kind=isk), allocatable :: sets(:,:)
	end type BSectors

	type :: BasisSet
		integer :: ntot
		logical :: xst
		! number of active sites
		integer :: n
		integer :: maxk ! max k of k-sub it has.
		integer, allocatable :: pntr(:)
		type(BSectors), allocatable :: sec(:)
	end type BasisSet

	!---------------------------------------	
	! for 13 different (N,m)
	!---------------------------------------	
	type :: Eigensystems
		!integer :: n,m ! same as Ham n,m
		integer :: ntot, n1,n2 ! size of hilbert space, =nrows, ncols=nevecs
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
		integer :: nsec ! number of degenerate sectors
		integer, allocatable:: ind(:) ! start index of all sectors
		double precision, allocatable :: esec(:)	 ! energy of deg sec
	end type Eigensystems
	!---------------------------------------	
	! hamiltonian for 13 diff (N,m)
	!---------------------------------------	
	type :: Ham
		logical :: xst,dense
		integer :: n,m,m1
		integer :: nev, ncv
		integer :: ntot ! HilbertSpace dimension
		integer :: nnz  ! no of non-zero elements
		integer :: srptr ! size of row pointers, not eq to ntot
		integer, allocatable :: col(:)
		integer, allocatable :: rowpntr(:)
		double precision, allocatable :: dat(:)
		!integer, allocatable :: row(:)
		 ! pointers to ksub sectors in rowpntr
		integer, allocatable :: spntr(:) 
	end type Ham

		! Number of Hilbert spaces and basis (same N)
	!	use extra to save repeated computations but 
	!	dont use too many to keep smaller memory usage
	! NBasisSets = 5, NHilbertSpaces=13
	integer :: NBasisSets, NHilbertSpaces

	! 5 BasisSet
	type(BasisSet), allocatable:: basis(:)
	! change name of HilbertSpace to something like eigensystem ??
	!	13 hamiltonians and eigensystems
	type(Eigensystems), allocatable:: eig(:)
	type(Ham), allocatable:: Hg(:)


	!---------------------------------------
	! hoppings: transition matrices
	!---------------------------------------	
	!---------------------------------------	
	type :: TransitionMatrix
		! transition matrix data
		integer :: nnz
		integer, allocatable :: row(:)
		integer, allocatable :: col(:)
	end type TransitionMatrix
	!---------------------------------------	
	type :: HoppingProcesses
		character(len=3) :: frmt ! format of sparse matrix; coo,  csr, diag, etc.
		integer :: nc, ns
		integer, allocatable, dimension(:) :: sites
		type(TransitionMatrix), allocatable :: ht(:,:) ! for channel, sites
	end type HoppingProcesses
	!---------------------------------------	


	!---------------------------------------
	! hoppings: transition amplitudes
	!---------------------------------------	
	type :: TransitionAmplitudes
		integer :: namp ! size of amp array
		double precision, allocatable :: amp(:) ! amplitudes
		! amplitudes squared of degen sectors
		integer :: nsec ! number of degenerate sectors
		!integer, allocatable:: ind(:) ! start index of all sectors
		double precision, allocatable :: amp2(:) 
	end type TransitionAmplitudes

	type :: QuantumTransitions
		integer :: nc,ns ! number of channels, sites
		type(TransitionAmplitudes), allocatable:: cs(:,:) ! channels and sites
	end type QuantumTransitions
	!---------------------------------------	
	! map from ih,ic to ia,icl
	!---------------------------------------	
	integer, dimension(26,4) :: maph ! ih  to ia
	integer, dimension(26,4) :: mapc ! ic to icl
	!	mapc only needed for channel 8,
	!	channel 8: swap channel 1 and 2 for amplitudes
	!---------------------------------------	


	!---------------------------------------
	! hoppings: transition rates
	!---------------------------------------	
	type :: TransitionRates
		integer :: nc,ns ! number of channels, sites		
		double precision :: r ! total
		double precision, allocatable :: rcs(:,:) !channel,site resolved
	end type TransitionRates

	! Master equation case, transition rates
	type :: MasterRates
		!integer :: nc,ns,nt ! number of channels, sites, timesteps		
		double precision, allocatable :: r(:) ! total
		double precision, allocatable :: rcst(:,:,:) !channel,site,timestep resolved
	end type MasterRates

	logical :: master ! alter CalAmp0/CalAmp behaviour when master eq is solved
	integer :: ntcoarse ! size of coarse time grid for rhovt and mrate etc
	integer :: ntmax ! mesolve integration grid size
	double precision :: dt ! incremental time, for integration of master eq in mesolve
	double precision :: wcut, J0 ! bath Ohmic spectral density parameters
	double precision, allocatable, dimension(:,:) :: rhovt
	integer, allocatable, dimension(:):: maprho ! location of diagonal elements 
	integer :: lrhov
	!integer :: lrhov
	!---------------------------------------	
	! define system
	!---------------------------------------	
	type :: System
		integer:: nsites
		integer:: n0,n1,n2 ! number of phi, acive, D
		integer, allocatable :: occ(:) ! occupancy of the site
	end type System
	type(System) :: sys
	!---------------------------------------	
	! define ways; sites for various hopping processes
	!---------------------------------------	
	type :: HoppingWays
		integer:: ns ! number of ways/sites for the given hop type
		! sites: 	d or phi if d/phi hop, 
		!				 	phi if d,phi pair annihilation
		!					right active if other is also active
		integer, allocatable :: sites(:) 
		!	active: active site or D if d,phi annihilation
		integer, allocatable :: active(:) ! active site(s) involved
	end type HoppingWays
	type(HoppingWays), dimension(26) :: ways ! 8 bulk, others contact


	!	26 hopping processes
	!type(HoppingProcesses), dimension(26) :: hop

	! not 26: in some cases, multiple hopping processes share amplitudes.
	! maph would give location for a given hopping process, ih ---> ihl
	type(QuantumTransitions), dimension(14) :: qt 
	double precision, allocatable :: psi(:,:) ! to store quantum state
	double precision, allocatable :: psi2(:,:) ! to store quantum state

	type(TransitionRates), dimension(26):: rate
	! set the format for sparse matrix
	!hop(1:4)%spfrmt = 'diagonal';

	type(MasterRates), dimension(26):: mrate

	logical :: leads ! include contact or not?

	double precision :: dw,g
	logical :: detuning

	!integer :: EvecKind

	! above this size, iterative arpack solver will be used
	!	to diagonalise Hamiltonians
	integer :: smalln ! key to set it in input: DirectSolverSize
	! max iterations for iterative solver
	integer :: diagmaxitr ! key to set it in input: DirectSolverSize

	integer :: ndsec ! number of degenerate sectors in the spectrum to calculate
	! ~ 3 seems fine, energetic penalty is not going to
	! allow much higher transitions anyway.

	logical :: commonbath
	logical :: pauli 
	integer :: ne, nsec
	! total hop rates of 26 hopping processes vs time (iterations)
	! write to file?!
	!double precision, allocatable :: HopRates(:,:)


	contains

!--------------- timer -----------------
	double precision function clock()
		real*4:: etime, tm(2)
		clock = etime( tm )
	return
	end function clock
!---------------------------------------
!	write rates in file 'rates.out'
!	order: hop type bulk hops(1-8), kappa(25),gamma(26),contacts(9:24)
	subroutine writeout(trajend,zt,ntrap)
	implicit none
	logical,intent(in) :: trajend
	integer,intent(in) :: zt,ntrap
	! local
	integer :: i
	logical, save :: first = .true.;
	integer :: iu

	if (first) then
	first = .false.
	i = 1+int((mexmax-mexmin)/dmex);

	iu=101
	if(ratesout)then
		open(iu,file='rates.out',action='write')
		if (onlydoped) then
			write(iu,*) i,ndw, 1+int((nelmax-nelmin-1)/dnelec),ntraj
		else
			write(iu,*) i,ndw, 1+int((nelmax-nelmin)/dnelec),ntraj
		endif
		close(iu)
	endif

	iu=102
	if(ztout)then
		open(iu,file='zt.out',action='write')
		if (onlydoped) then
			write(iu,*) i,ndw, 1+int((nelmax-nelmin-1)/dnelec),ntraj
		else
			write(iu,*) i,ndw, 1+int((nelmax-nelmin)/dnelec),ntraj
		endif
		close(iu)
	endif

	iu=103
	if(ntrapout)then
		open(iu,file='traps.out',action='write')
		if (onlydoped) then
			write(iu,*) i,ndw, 1+int((nelmax-nelmin-1)/dnelec),ntraj
		else
			write(iu,*) i,ndw, 1+int((nelmax-nelmin)/dnelec),ntraj
		endif
		close(iu)
	endif
	
	endif
	!------------------------------------------------------

	iu=104
	if(ratesout) then
		open(iu,file='rates.out',action='write',position='append')
		if(periodic .or. onlybulk) then
			write(iu,'(f15.8,5x,10f15.8)')
     . sum(rate(:)%r),rate(1:8)%r,rate(25:26)%r
		else
			write(iu,'(f15.8,5x,26f15.8)')
     . sum(rate(:)%r),rate(1:8)%r,rate(25:26)%r,rate(9:24)%r
		endif
		close(iu)
	endif

	iu=105
	if(trajend .and. ztout)then
		open(iu,file='zt.out',action='write',position='append')
		write(iu,*) zt
		close(iu)
	endif

	iu=106
	if(trajend .and. ntrapout)then
		open(iu,file='traps.out',action='write',position='append')
		write(iu,*) ntrap
		close(iu)
	endif

	
	return
	end 	subroutine writeout
!---------------------------------------
	! for uncoupled case:
	! pick one of ksub with max possible excited sites
	! for nog case, we never need N_ex > 0 unless
	! we want to have some up spins
	subroutine chooseipsi(ipsi,Einit);
	implicit none
	integer, intent(out):: ipsi
	double precision, intent(out):: Einit
	! local
	integer :: it

	it = mapt%map(1);
	! here, the last one, just because
	! it is easy to get it's index
	ipsi = eig(it)%ntot
	Einit = eig(it)%eval(ipsi)
	return
	end subroutine chooseipsi

	end module
