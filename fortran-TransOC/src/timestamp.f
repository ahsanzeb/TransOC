


	Subroutine timestamp
!--------------------------------------------------------------------
! 	prints welcome message with date and time
!--------------------------------------------------------------------

	implicit none

	integer	:: dt(8)
	character*10 :: bb(3)
	double precision :: dummy

!--------------------------------------------------------------------
	call date_and_time(bb(1), bb(2), bb(3), dt)

	write(6,*)
	write(6,*)"**************************"//
     .		"******************************"
	write(6,'(/,a)')"			WELCOME ! "
	write(6,'(/,a)')"		   TransOC Started !"
        write(6,'(/,a)')
     .        '	        . . . Running v 0.0 . . . '

	write(6,'(/,a,i2,a,i2,2a,i2,a,i2,a,i4,a)')
     . 	"	   ",
     .                   dt(5)," Baj kar ",dt(6), " min",
     .               " DATE :  ",dt(3),"/",dt(2),"/",dt(1)
	write(6,*)"**************************"//
     .		"******************************"




	! start random sequence with seed = milliseconds in current time
	call srand(dt(8)) !rand(dt(8))

	END Subroutine timestamp

