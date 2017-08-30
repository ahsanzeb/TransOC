
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

	!	na, nx lists for 13 Hilbert spaces
	integer(kind=1), dimension(13):: nalist,nxlist !CHECK  not used????

	! index for jump, channel, site
	integer itype

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
	type :: HilbertSpace
		integer :: ntot, n1,n2
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
	end type HilbertSpace
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
		!double precision, allocatable :: dat(:)
		!	dat: t's will be multiplied latter
		!	dat ====>  1 always; so no need to store it.
		!integer(kind=4), allocatable :: map(:,:)
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
		integer, allocatable :: sites(:)
	end type HoppingWays
	type(HoppingWays) :: ways


	! 5 BasisSet
	type(BasisSet), dimension(5) :: basis
	! change name of HilbertSpace to something like eigensystem ??
	!	13 hamiltonians and eigensystems
	type(HilbertSpace), dimension(13) :: hspace
	type(Ham), dimension(13) :: Hg
	!	26 hopping processes
	type(HoppingProcesses), dimension(26) :: hop
	! not 26: in some cases, multiple hopping processes share amplitudes.
	! maph would give location for a given hopping process, ih ---> ihl
	type(QuantumTransitions), dimension(14) :: qt 
	double precision, allocatable :: psi(:) ! to store quantum state

	! set the format for sparse matrix
	!hop(1:4)%spfrmt = 'diagonal';


	double precision :: dw,g
	logical :: detuning





	!integer :: EvecKind


	end module
