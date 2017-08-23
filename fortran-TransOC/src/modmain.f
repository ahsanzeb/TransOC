
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
	integer m
	! index for jump, channel, site
	integer ij, ic, is
	! 
	integer*4, allocatable :: pntr(:)
	integer :: ind													! use as work array
	integer, allocatable :: comb(:) 				! use as work array
	integer, allocatable :: combset(:,:) 	! use as work array



	
	! for 13 different (N,m)
	type :: HilbertSpace
		integer :: ntot
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
	end type HilbertSpace

	! for 5 values of N; smaller m can use larger m's data
	type :: BasisSet
		integer, allocatable :: pntr(:)
		integer, allocatable :: sets(:,:)
	end type BasisSet






	end module
