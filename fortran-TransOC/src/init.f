

	subroutine init()
	use modmain
	implicit none


	! initialise mapb, mapt
	allocate(mapb%map(5))
	allocate(mapb%cal(5))
	allocate(mapt%map(13))
	allocate(mapt%cal(13))
	mapb%nnu = 5; mapt%nnu = 13 ! calc all 
	mapb%map = (/ 1,2,3,4,5 /)
	mapb%cal = (/ 1,2,3,4,5 /)
	mapt%map = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	mapt%cal = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)

















	end subroutine init
	
