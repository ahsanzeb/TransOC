


	Subroutine timestamp()
!--------------------------------------------------------------------
! 	prints welcome message with date and time
!--------------------------------------------------------------------

	implicit none

	integer	:: dt(8), x
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
		! assume prog will run on the same day
		x= dt(3)*86400+dt(5)*3600 + dt(6)*60 + dt(7) ! sec
		x = x - (dti(3)*86400+dti(5)*3600 + dti(6)*60 + dti(7));
		dt(3) = x/86400; ! days
		x = x - dt(3)*86400;
		dt(5) = x/3600; ! hours
		x = x - dt(5)*3600;
		dt(6) = x/60; ! min
		x = x - dt(6)*60;
		dt(7) = x; ! sec

		write(6,
     . '("Time taken: ",i2,":",i2.2,":",i2.2,":",i2.2,"(d:h:m:s)")')
     .   dt(3),dt(5:7)
	endif

	return
	END Subroutine timestamp
		

