
! Copyright (C) 2017 M. Ahsan Zeb
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

	module modmain
	implicit none
	
	! total number of sites in the system
	integer(kind=1):: nsites
	! active sites
	integer(kind=1):: na
	! number of excitatons
	integer(kind=1):: nx

	! list of active sites in order they appear in the basis 
	integer, allocatable :: ASites(:)

	! PermSym: if all sites have same w0 and g, and AlwaysLP=True, i.e,
	!	 	a quick relaxation to Lower polariton state after every hop,
	!		then we are always in the LP state when any type of hopping
	!		occurs, all sites will be similar and we can use amplitudes
	!		calculated for one site for others.
	logical:: PermSym

	!	na, nx lists for 13 Hilbert spaces
	integer(kind=1), dimension(13):: nalist,nxlist !CHECK  not used????

	! index for jump, channel, site
	integer:: itype

	! hopping parameters for Homo-Homo, etc
	!	set defauls after reading input file
	double precision, dimension(8):: tpar
	
	! cavity and exciton loss rates
	double precision :: kappa, gamma	

	! hopping parameters
	double precision, dimension(26,4):: hpar
	! energy changes due to contact barriers, applied field, etc
	double precision, dimension(26,4):: dqc
	! bare exciton energy, right and left contact barriers
	!	Electric field energy Er assumed constant for all hops
	double precision:: w0, Ebr,Ebl, Er
	! beta = 1/KbT for penalty function
	double precision:: beta	
	! block injection of electron/holes?
	logical, dimension(2):: BlockInjection
	!	periodic boundary conditions?
	logical:: periodic ! XXXXXXX make it about geometry of the system
	! only bulk processes? with or without periodic boundary conditions
	logical :: onlybulk ! no contacts?
	logical :: nolosses ! kappa=0=gamma ?
	double precision:: Einit ! intial energy in every iteration
	!	--------------- maps ---------------
	integer(kind=1), dimension(13)::
     . dna	= (/ 0,0,0,2,2,2,-2,-2,-2,1,1,-1,-1 /);
	integer(kind=1), dimension(13)::
     . dnx = (/ 0,-1,1,1,0,2,-1,-2,0,1,0,-1,0 /);
	integer(kind=1), dimension(13)::
     . ibs = (/ 3,3,3,5,5,5,1,1,1,4,4,2,2 /);
	integer(kind=1), dimension(5)::
     . dns = (/ -2,-1,0,1,2 /); 
	integer(kind=1), dimension(26,4)::
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
	!integer(kind=1), dimension(5):: mapb 	! map for ib
	!integer(kind=1), dimension(13):: mapt	! map for itype

	type :: BTMaps
		integer :: nnu
		!integer :: ntot
		integer(kind=1), allocatable :: map(:)	 ! map for ib or it
		integer(kind=1), allocatable :: cal(:) ! which ib to be calculated?
		! which itypes for a given ib; ntb,grouptb only for mapt
		integer(kind=1), dimension(5) :: ntb 
		integer(kind=1), dimension(5,13) :: grouptb 
	end type BTMaps

	type(BTMaps) :: mapb, mapt
	!	------------------------------


	logical crosshops ! L-H and H-L cross hops allowed? 



	!---------------------------------------	
	! for 5 values of N; smaller m can use larger m's data
	! basis sectors for a given na,k-up spins
	!---------------------------------------
	type :: BSectors
		integer :: ntot ! will see later if these ntot are needed?
		integer(kind=1), allocatable :: sets(:,:)
	end type BSectors

	type :: BasisSet
		integer :: ntot
		integer :: maxk ! max k of k-sub it has.
		integer(kind=4), allocatable :: pntr(:)
		type(BSectors), allocatable :: sec(:)
	end type BasisSet

	!---------------------------------------	
	! for 13 different (N,m)
	!---------------------------------------	
	type :: Eigensystems
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
		integer(kind=4) :: ntot
		integer(kind=4), allocatable :: row(:)
		integer(kind=4), allocatable :: col(:)
		double precision, allocatable :: dat(:)
	end type Ham
	!---------------------------------------
	! hoppings: transition matrices
	!---------------------------------------	
	!---------------------------------------	
	type :: TransitionMatrix
		! transition matrix data
		integer(kind=4) :: nnz
		integer(kind=4), allocatable :: row(:)
		integer(kind=4), allocatable :: col(:)
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
		integer(kind=4) :: namp ! size of amp array
		double precision, allocatable :: amp(:) ! amplitudes
		! amplitudes squared of degen sectors
		integer :: nsec ! number of degenerate sectors
		!integer, allocatable:: ind(:) ! start index of all sectors
		double precision, allocatable :: amp2(:) 
	end type TransitionAmplitudes

	type :: QuantumTransitions
		integer(kind=4) :: nc,ns ! number of channels, sites
		type(TransitionAmplitudes), allocatable:: cs(:,:) ! channels and sites
	end type QuantumTransitions
	!---------------------------------------	
	! map from ih,ic to ia,icl
	!---------------------------------------	
	integer(kind=4), dimension(26,4) :: maph ! ih  to ia
	integer(kind=4), dimension(26,4) :: mapc ! ic to icl
	!	mapc only needed for channel 8,
	!	channel 8: swap channel 1 and 2 for amplitudes
	!---------------------------------------	


	!---------------------------------------
	! hoppings: transition rates
	!---------------------------------------	
	type :: TransitionRates
		integer(kind=4) :: nc,ns ! number of channels, sites		
		double precision :: r ! total
		double precision, allocatable :: rcs(:,:) !channel,site resolved
	end type TransitionRates

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
		integer:: ns ! number of sites
		!XXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!XXXXXXXXXXXXXXXXXXX
		!XXXXXXXXXXXXXXXXXXXXXXXXX
		integer, allocatable :: sites(:) ! change to sensible....
		integer, allocatable :: active(:) ! active site(s) involved
	end type HoppingWays
	type(HoppingWays), dimension(26) :: ways ! 8 bulk, others contact

	! 5 BasisSet
	type(BasisSet), dimension(5) :: basis
	! change name of HilbertSpace to something like eigensystem ??
	!	13 hamiltonians and eigensystems
	type(Eigensystems), dimension(13) :: eig
	type(Ham), dimension(13) :: Hg

	!	26 hopping processes
	!type(HoppingProcesses), dimension(26) :: hop

	! not 26: in some cases, multiple hopping processes share amplitudes.
	! maph would give location for a given hopping process, ih ---> ihl
	type(QuantumTransitions), dimension(14) :: qt 
	double precision, allocatable :: psi(:) ! to store quantum state

	type(TransitionRates), dimension(26):: rate
	! set the format for sparse matrix
	!hop(1:4)%spfrmt = 'diagonal';



	double precision :: dw,g
	logical :: detuning





	!integer :: EvecKind


	end module
