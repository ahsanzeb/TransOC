

	module selection
	implicit none

	public ::ihSelect, icsSelect 

	contains
!***************************************************
	integer function ihSelect()
	use modmain, only: rate
	implicit none
	! local
	integer:: ih
	double precision, dimension(27):: rlist ! accumulated rates
	double precision:: eta


	!write(*,*)"rates= ",rate(1:8)%r
	! select ih stochastically
	rlist(1) = 0.0d0;
	do ih=1,26
		rlist(ih+1) = 	rlist(ih) + rate(ih)%r
	end do
	! eta, random real in range 0-rlist(27)
	eta = rand(0) * rlist(27);
	! location of eta on accumulated rates
	ihSelect = -1;
	do ih=1,26
		if (eta > rlist(ih) .and. eta .le. rlist(ih+1) ) then
			ihSelect = ih;
			exit
		endif
	end do	
	! error?
	if (ihSelect == -1) then
		write(*,*) "ihSelect: something wrong.... "
		write(*,*) "ihSelect: eta = ",	eta		
		write(*,*) "ihSelect: rlist = ",	rlist	
		stop
	endif

	return
	end function ihSelect
!***************************************************
	subroutine icsSelect(ih,ic,is)
	use modmain, only: rate, ways,PermSym,crosshops
	implicit none
	integer, intent(in) :: ih
	integer, intent(out) :: ic,is
	! local
	double precision, allocatable:: rlist(:) ! accumulated rates
	integer :: ntot,nc,ns,which,i
	double precision:: eta

	nc = 1;
	if (ih .le. 8) then
		nc = 2
		if (crosshops) nc = 4
	endif
	
	ns = 1; ! only in=1-4, 7,8
	if(.not. PermSym) then	
		if ( (ih .le. 4) .or. 
     . ih==7 .or. ih==8 .or. ih==26) ns = ways(ih)%ns
	endif
	
	ntot = ns * nc + 1;
	allocate(rlist(ntot))

	i = 1;
	rlist(1) = 0.0d0;
	do is=1,ns
		do ic =1,nc
			rlist(i+1) = rlist(i) + rate(ih)%rcs(ic,is)
			i = i + 1;
		enddo
	enddo
			
	! eta, random real in range 0-rlist(27)
	eta = rand(0)* rlist(ntot);
	! location of eta on accumulated rates
	which = -1;
	do i=1,ntot-1
		if (eta > rlist(i) .and. eta .le. rlist(i+1) ) then
			which = i;
			exit
		endif
	end do	
	! error?
	if (which == -1) then
		write(*,*) "icsSelect: something wrong.... "
		write(*,*) "icsSelect: rlist = ",	rlist	
		stop
	endif


	if(ns == 1) then
		ic = which;
		ns = ways(ih)%ns;
		if(ns > 1) then
			is = int(rand(0)*ns) + 1; ! when PermSym
		else
			is = 1;
		endif
	elseif(nc == 1) then
			ic = 1;
			is = which;
	else
		is = 1 + mod((which-1),nc); ! shifted => Mathematica Mod[which,nc,1]
		ic = which - (is-1)*nc;
	endif


	write(*,*)"ih,ic,is = ",ih,ic,is



	return
	end subroutine icsSelect
!***************************************************

	end module selection
