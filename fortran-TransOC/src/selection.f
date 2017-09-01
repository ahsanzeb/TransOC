

	module selection
	use modmain, only: rate
	implicit none


	contains

!***************************************************
	integer function ihSelect()
	use modmain, only: rate
	implicit none
	! local
	integer:: ih
	double precision, dimension(27):: rlist ! accumulated rates

	! select ih stochastically
	rlist(1) = 0.0d0;
	do ih=1,26
		rlist(ih+1) = 	rlist(ih) + rate(ih)%r
	end do
	! eta, random real in range 0-rlist(27)
	eta = random_number() * rlist(27);
	! location of eta on accumulated rates
	ihSelect = -1;
	do ih=1,26
		if (eta > rlist(ih) .and. eta .le. rlist(ih) ) then
			ihSelect = ih;
			exit
		endif
	end do	
	! error?
	if (ihSelect == -1) then
		write(*,*) "ihSelect: something wrong.... "
		write(*,*) "ihSelect: rlist = ",	rlist	
		stop
	endif

	return
	end function ihSelect
!***************************************************









	end module selection
