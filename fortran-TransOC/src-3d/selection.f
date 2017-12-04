

	module selection
	
	implicit none

	public ::ihSelect, icsSelect,getpsi2

	contains

!***************************************************
	integer function ihSelect()
	use modmain, only: rate, nproc
	implicit none
	! local
	integer:: ih
	double precision, dimension(43):: rlist ! accumulated rates
	double precision:: eta

	!write(*,*)"rates= ",rate(:)%r
	! select ih stochastically
	rlist(1) = 0.0d0;
	do ih=1,nproc
		rlist(ih+1) = 	rlist(ih) + rate(ih)%r
	end do
	! eta, random real in range 0-rlist(27)
	eta = rand(0) * rlist(nproc+1);
	! location of eta on accumulated rates
	ihSelect = -1;
	do ih=1,nproc
		if (eta .ge. rlist(ih) .and. eta .lt. rlist(ih+1) ) then
			ihSelect = ih;
			exit
		endif
	end do	
	! error?
	if (ihSelect == -1 .or. ihSelect>nproc) then
		write(*,*) "ihSelect: something wrong.... ihSelect=",ihSelect
		write(*,*) "ihSelect: eta = ",	eta		
		write(*,'(43G18.5)') "ihSelect: rlist = ",	rlist	
		!stop
	endif

	return
	end function ihSelect
!***************************************************
	subroutine icsSelect(ih,ic,is)
	use modmain, only: rate, ways,PermSym,crosshops,nog,sys,
     .               vrh,coulomb
	implicit none
	integer, intent(in) :: ih
	integer, intent(out) :: ic,is
	! local
	double precision, allocatable:: rlist(:) ! accumulated rates
	integer :: ntot,nc,ns,which,i
	double precision:: eta

	!write(*,*)"rate%rcs= ",rate(ih)%rcs(:,:)
	if(ih==25) then ! output not significant, not used anywhere
		ic=1; is=1;
		return
	endif

	nc = 1;
	if (ih .le. 8 .or. (ih .ge. 27 .and. ih <=34 )) then ! ih=1-8, 27-34
		nc = 2
		if (crosshops) nc = 4
	endif

	ns = 1; ! just for selection; not actual ns 
	if(.not. PermSym .or. vrh .or. coulomb) then	 ! also applies to nog case
		ns = ways(ih)%ns 
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

	!write(*,*) "ih, ways(ih)%ns = ",ih,ways(ih)%ns
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
		write(*,*) "icsSelect: something wrong.... ih,ic,is =",ih,ic,is
		write(*,*) "icsSelect: rlist = ",	rlist	
		stop
	endif


	if(ns == 1) then
		ic = which;
		ns = ways(ih)%ns; ! resume
		is = int(rand(0)*ns) + 1;
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
		write(*,*)"sel: is =",is
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

!--------------------------------------
	logical function MemberQ4(a,la,x)
	implicit none
	integer, intent(in):: la,x
	integer, dimension(la), intent(in):: a
	integer:: i
	
	MemberQ4 = .false.
	do i=1,la
		if (a(i) == x) then
			MemberQ4 = .true.
			exit
		endif
	end do
	return
	end function
!--------------------------------------


	end module selection
