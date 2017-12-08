
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

	integer :: nproc
	logical :: EqualDistr

	logical :: impurity ! dopant or trap level: 1LS with occ=0,1 only
	double precision :: Eimp ! enegy of impuity level (w.r.t HOMO levels at 0)
	integer :: impocc ! occ of impurity level

	logical :: vqout ! output Vq, q at all iterations?
	logical:: coulomb ! include coulomb interaction?
	double precision :: epsr != 4.0d0; !dielectric constant of organic material
	double precision :: Kq != 0.3556d0; ! Coulomb law constant K for the medium for rij in nm, q's in |e|. Units: eV.nm/e
										! Kq = (1/epsr=4.0) * (8.89d9*1.6d-19/1d-9)
	
	double precision :: Etotq ! total Coulomb energy
	type :: CoulombEnergy
		double precision, allocatable, dimension(:) :: dEq ! size= ways(ih)%ns
	end type CoulombEnergy
	type(CoulombEnergy) :: Ecoul(42) ! size = types of hops

	logical :: vrh ! positional disorder? variable range hopping
	double precision, allocatable, dimension(:,:):: bondlengths
	integer, dimension(34):: signEr =(/
     .   1, -1, -1, 1, 1, -1, -1, 1,
     .   1, 1, -1, -1, -1, -1, 1, 1,
     .   -1, 1, -1, 1, 1, -1, 1, -1,
     .   0, 0,
     .   0, 0, 0, 0, 0, 0, 0, 0 /);
	!double precision, dimension(34):: Efieldh
	! average of r_nns, std dev in r_nns gaussian distribution
	double precision :: a0, sigma0, nsigma, dinvl

	logical :: debug

	! total number of sites in the system
	integer:: nsites

	double precision :: wcut, J0 ! bath Ohmic spectral density parameters
	logical :: simplepf 
	
	! number of electrons 
	! DopingRange
	integer :: nel,nelmin,nelmax,dnelec ! nel will stay fixed if no contacts

	! net charge on the system
	integer :: Qnet
	! Charging energy prop constant
	double precision:: Eq
	! charging energy for contact hops
	double precision:: dEQs(16)
	
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
	double precision, dimension(42,4):: ts
	! energy changes due to contact barriers, applied field, etc
	double precision, dimension(42,4):: dqc
	! hopping parameters
	double precision:: th, tl, tlh, thl
	double precision:: JhR, JlR, JhL, JlL
	! bare exciton energy, right and left contact barriers
	!	Electric field energy Er assumed constant for all hops
	double precision:: w0, Ebr,Ebl, Er, Exb
	double precision:: Ermin, Ermax, dEr
	integer :: ner

	double precision:: Eblmin, Eblmax,Ebrmin, Ebrmax, dEbr, dEbl
	double precision:: gmin, gmax, Tmin, Tmax,dtemp,dg, kb=8.33d-5; !(eV/K)
	integer :: nEbl,	 nEbr, nTemp, ng

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
	integer, dimension(42,4)::
     . itypes = reshape( ( / !bulk: left, right
     . 1, 1, 2, 3, 1, 1, 2, 3, 1, 1,
     . 2, 3, 1, 1, 2, 3, 4, 4, 5, 6,
     . 4, 4, 5, 6, 7, 7, 8, 9, 7, 7,
     . 8, 9, ! contact:
     . 11, 11, 11, 11, 10, 10, 
     . 10, 10, 11, 11, 11, 11, 10, 10,
     . 10, 10, 11, 11, 11, 11, 10,
     . 10, 10, 10, 11, 11, 11, 11,
     . 10, 10, 10, 10, 13, 13, 13,
     . 13, 13, 13, 13, 13, 12, 12,
     . 12, 12, 12, 12, 12, 12, 13,
     . 13, 13, 13, 13, 13, 13, 13,
     . 12, 12, 12, 12, 12, 12, 12,
     . 12, ! kappa, gamma:
     . 2, 2, 2, 2, 2, 2, 2, 2, ! bulk, up/down:
     . 1, 1, 2, 3, 1, 1, 2, 3, 1, 1,
     . 2, 3, 1, 1, 2, 3, 4, 4, 5, 6,
     . 4, 4, 5, 6, 7, 7, 8, 9, 7, 7,
     . 8, 9, ! impurity: ih=35:42 <===> left contact hops 11,12,15,16,21:24
     . 11, 11, 11, 11, 10, 10, 10, 10,
     . 11, 11, 11, 11, 10, 10, 10, 10,     
     . 13, 13, 13, 13, 13, 13, 13, 13,
     . 12, 12, 12, 12, 12, 12, 12, 12 /),
     . (/ 42,4 /), order=(/2,1/) ); ! ih=1-8 <==> 26-34 for in-plane hops
 
	integer, dimension(8) ::
     .                  ihdphops=(/1,2,3,4,27,28,29,30/)

	integer, dimension(4) ::
     .                  ihcreat = (/7,8,33,34/)
	integer, dimension(4) ::
     .                  ihannih = (/5,6,31,32/)

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
	! to manage cases of e/h injection barriers, make devices types etc.
	type :: Device1
		double precision, allocatable, dimension(:,:) :: Eb
	end type Device1
	type :: Device0
		integer :: ntyps, ntot
		integer, allocatable, dimension(:) :: typ, ndev
		type(Device1), allocatable, dimension(:) :: X
		integer, allocatable, dimension(:,:) :: idv
	end type Device0
	type(Device0)	:: dev
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
	! max dim 1 = 42 if no permsym and impurity;
	! can be made allocatable, and allocated in calmaphc()
	integer, dimension(42,4) :: maph ! ih  to ia 
	integer, dimension(42,4) :: mapc ! ic to icl
	!integer, allocatable, dimension(:,:) :: maph, mapc ! ih  to ia; ic to icl
	
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
	!---------------------------------------	
	! define system
	!---------------------------------------	
	type :: System
		integer:: nsites
		integer:: n0,n1,n2,nimp ! number of phi, acive, D
		integer, allocatable :: occ(:) ! occupancy of the site
		double precision, allocatable, dimension(:,:) :: r
		double precision, allocatable, dimension(:) :: q0, q ! integer?
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
		double precision, allocatable :: rij(:) ! distance for the hop
	end type HoppingWays
	type(HoppingWays), dimension(42) :: ways ! 8 bulk, others contact, + 8 impurity


	!	26 hopping processes
	!type(HoppingProcesses), dimension(26) :: hop

	! not 26: in some cases, multiple hopping processes share amplitudes.
	! maph would give location for a given hopping process, ih ---> ihl
	type(QuantumTransitions), dimension(21) :: qt ! max possible 21 if no permsym and impurity
	double precision, allocatable :: psi(:,:) ! to store quantum state
	double precision, allocatable :: psi2(:,:) ! to store quantum state

	type(TransitionRates), dimension(42):: rate ! max possible 42=34+8 if impurity
	! set the format for sparse matrix
	!hop(1:4)%spfrmt = 'diagonal';

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

!========================================
	function statistics(x,l)
	implicit none
	integer,intent(in) :: l
	double precision, dimension(l),intent(in):: x
	double precision, dimension(2):: statistics
	
	! local
	integer:: i
	double precision:: av, var

	statistics(1) = sum(x)/l;
	statistics(2) = sum((x - statistics(1))**2)/l
	return
	end function statistics
!========================================

!==============================
	subroutine getnphoton(ipsi,itype,nphoton)
	implicit none
	integer, intent(in):: ipsi,itype
	integer, intent(out) :: nphoton
	integer:: n,m,m1,ib,ibi,i
	logical :: found

	!write(*,*)'nphoton: ipsi,itype = ',ipsi,itype
	n = na + dna(itype); ! no of active sites
	m = nx + dnx(itype); ! no of excitations
	m1 = min(m,n); ! max no of up spins possible
	ibi = ibs(itype);
	ib = mapb%map(ibi);	
	found=.false.;
	do i=1,m1+2
	!write(*,*)'ipsi, m1+2, i = ',ipsi, m1+2, i
	!write(*,*)basis(ib)%pntr(i), basis(ib)%pntr(i+1)
	
	if(ipsi > basis(ib)%pntr(i) .and. 
     .  ipsi <= basis(ib)%pntr(i+1) ) then   
		nphoton = m-i+1;
		found=.true.;
		exit
	endif
	enddo
	if(.not. found) then
		write(*,*)"main: nphoton not found!!!"
		stop
	endif

	return
	end subroutine getnphoton
!==============================

	end module
