
! Copyright (C) 2017 M. Ahsan Zeb
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

	module modmain
	implicit none
	
	! total number of sites in the system
	integer nmax
	! active sites
	integer na
	! number of excitatons
	integer nx
	!	na, nx lists for 13 Hilbert spaces
	integer, dimension(13):: nalist,nxlist
	integer, dimension(13):: dna=[0,0,0,2,1,2,-2,-2,-2,1,1,-1,-1];
	integer, dimension(13):: dnx=[0,-1,1,1,0,2,-1,-2,0,1,0,-1,0];
	integer, dimension(5):: nalist5
	integer, dimension(5):: dnalist5=[-2,-1,0,1,2];
	integer, dimension(13):: nainds=[3,3,3,5,4,5,1,1,1,4,4,2,2];
	! index for jump, channel, site
	integer ij, ic, is
	! 
	!integer*4, allocatable :: pntr(:)
	!integer :: ind													! use as work array
	!integer(kind=1), allocatable :: comb(:) 				! use as work array
	!integer(kind=1), allocatable :: combset(:,:) 	! use as work array

	logical crosshops ! L-H and H-L cross hops allowed? 
	
	! for 5 values of N; smaller m can use larger m's data
	! basis sectors for a given na,k-up spins
	type :: BSectors
		integer :: ntot ! will see later if these ntot are needed?
		integer(kind=1), allocatable :: sets(:,:)
	end type BSectors

	type :: BasisSet
		integer :: ntot
		integer(kind=4), allocatable :: pntr(:)
		type(BSectors), allocatable :: sec(:)
	end type BasisSet


	! for 13 different (N,m)
	type :: HilbertSpace
		integer :: ntot
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
	end type HilbertSpace

	type(BasisSet), dimension(5) basis
	type(HilbertSpace), dimension(13) hspace

	integer(kind=1), allocatable :: sites(:)






	end module
