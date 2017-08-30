

	subroutine init()
	use modmain
	implicit none

	!-----------------------------------------
	! initialise mapb, mapt
	!-----------------------------------------
	allocate(mapb%map(5))
	allocate(mapb%cal(5))
	allocate(mapt%map(13))
	allocate(mapt%cal(13))
	mapb%nnu = 5; mapt%nnu = 13 ! calc all 
	mapb%map = (/ 1,2,3,4,5 /)
	mapb%cal = (/ 1,2,3,4,5 /)
	mapt%map = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	mapt%cal = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	!-----------------------------------------


	!-----------------------------------------
	! calculate maphc, the map from hopping index to amplitude index,
	!	and channel to channel location, needed for channel 8 only,
	! but does not cost much to make map for all 26 x 4.
	call calmaphc()
	!-----------------------------------------















	end subroutine init
	
