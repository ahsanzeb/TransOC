


	Subroutine timestamp
!--------------------------------------------------------------------
! 	prints welcome message with date and time
!--------------------------------------------------------------------

	implicit none

	integer	:: dt(8)
	character*10 :: bb(3)
	double precision :: dummy
	logical, save :: start = .true.;
	integer, save	:: dti(8)
!--------------------------------------------------------------------
	call date_and_time(bb(1), bb(2), bb(3), dt)

	if (start) then
		dti = dt;
		start = .false.;
		! start random sequence with seed = milliseconds in current time
		call srand(dt(8)) !rand(dt(8))
		write(6,*)
		!write(6,*)"**************************"//
    ! .		"******************************"
		write(6,'(/,a)')"			TransOC Started! "
		!write(6,'(/,a)')" TransOC: Transport in Organic Cavities"
		!write(6,'(/,a)')
    ! .        '	        . . . Running v 0.0 . . . '
		write(6,'(/,a,i2,a,i2,2a,i2,a,i2.2,a,i4,a)')
     . 	"	   ",
     .                   dt(5)," Baj kar ",dt(6), " min",
     .               " aor tarekh hai ",dt(3),"/",dt(2),"/",dt(1)
		write(6,*)"         -----------------"//
     .		"----------------------------"
	else
		dt= dt-dti;
		write(6,
     . '("Time taken: ",i2,":",i2.2,":",i2.2,":",i2.2,"(d:h:m:s)")')
     .   dt(3),dt(5:7)
	endif




	END Subroutine timestamp

