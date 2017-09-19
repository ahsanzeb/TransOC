

	module selection
	
	implicit none

	public ::ihSelect, icsSelect,getpsi2

	contains

!***************************************************
	integer function ihSelect()
	use modmain, only: rate
	implicit none
	! local
	integer:: ih
	double precision, dimension(27):: rlist ! accumulated rates
	double precision:: eta

	!write(*,*)"rates= ",rate(:)%r
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
		if (eta .ge. rlist(ih) .and. eta .lt. rlist(ih+1) ) then
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
	use modmain, only: rate, ways,PermSym,crosshops,nog,sys
	implicit none
	integer, intent(in) :: ih
	integer, intent(out) :: ic,is
	! local
	double precision, allocatable:: rlist(:) ! accumulated rates
	integer :: ntot,nc,ns,which,i
	double precision:: eta

	
	!write(*,*)"rate%rcs= ",rate(ih)%rcs(:,:)

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
	
	if(nog)then
		nc = 4; ns = sys%nsites;
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


	!write(*,*) "icsSel: rlist=",rlist
			
	! eta, random real in range 0-rlist(27)
	eta = rand(0)* rlist(ntot);
	! location of eta on accumulated rates
	which = -1;
	do i=1,ntot-1
		if (eta .ge. rlist(i) .and. eta .lt. rlist(i+1) ) then
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
			ic = 1 + mod((which-1),nc); ! shifted => Mathematica Mod[which,nc,1]
			is = (nc-1+which)/nc; !floor((nc-1+which)/nc)
	endif
	!write(*,*)"ih,ic,is = ",ih,ic,is



	return
	end subroutine icsSelect
!***************************************************
!	superposition of states in the
! first excited degen sector above the lowest state
	subroutine getpsi2(ih,ic,is)
	use modmain, only: itypes,mapt,maph,eig,qt,
     .       psi2,Einit2,PermSym,debug
	implicit none
	integer, intent(in):: ih,ic,is
	! local
	integer :: i,ia,it,itp,nsec,s,iss

	if(PermSym) then
		iss=1;
	else
		iss=is
	endif

	itp = itypes(ih,ic)
	it = mapt%map(itp);
	ia = maph(ih,ic);
	nsec = eig(it)%nsec;

	if (allocated(Einit2)) deallocate(Einit2)
	allocate(Einit2(nsec))
	Einit2 = eig(it)%esec(1:nsec); ! second deg sector

	if (allocated(psi2)) deallocate(psi2)
	allocate(psi2(nsec,eig(it)%n1))
	psi2 = 0.0d0

	!write(*,*)"=====> sel: itp,it,ia,nsec=", itp,it,ia,nsec
	!write(*,*)shape(eig(it)%ind)
	!write(*,*)shape(eig(it)%evec)
	!write(*,*)shape(qt(ia)%cs(ic,is)%amp)
	!write(*,*)eig(it)%ind
	!write(*,*)eig(it)%eval

	! debugging, remove later
	if(debug) then
		write(*,*)eig(it)%ind
		write(*,*)eig(it)%evec
		write(*,*)qt(ia)%cs(ic,iss)%amp
		write(*,*) "selection: it, eig(it)%nsec=", it, nsec
		write(*,*)"sel: eig(it)%evec ",shape(eig(it)%evec)
		write(*,*)"sel: qt(ia)%cs ",shape(qt(ia)%cs)
		write(*,*)"sel: qt(ia)%cs(ic,is)%amp ",
     .				shape(qt(ia)%cs(ic,iss)%amp)
		write(*,*)"sel: is = is"
	endif

	do s=1,nsec
		!write(*,*)"i1, i2 = ",eig(it)%ind(s),eig(it)%ind(s+1)-1
		!write(*,*)"............  s = ",s
		do i=eig(it)%ind(s),eig(it)%ind(s+1)-1 ! second degenerate sector
			psi2(s,:) = psi2(s,:)  +
     .  eig(it)%evec(:,i) * qt(ia)%cs(ic,iss)%amp(i)
		enddo
	enddo

	return
	end subroutine getpsi2
!***************************************************

	end module selection
