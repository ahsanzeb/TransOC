
	module modways
	use modmain, only: sys, ways, periodic
	implicit none

	public :: UpdateWays, UpdateOcc
	private:: WhichBulkHop

	contains
!**********************************************
	subroutine UpdateWays()
	implicit none
	! local
	integer(kind=1):: ntot, p,la,lo,is,is1
	integer, dimension(26) :: nps
	integer, allocatable, dimension(:,:) :: act,oth

	! initialise
	nps = 0	
	! copy current sys
	ntot = sys%nsites
	allocate(act(26,ntot))
	allocate(oth(26,ntot))
	act = -2
	oth = -2

	!write(*,*)'sys%nsites====>',  sys%nsites
	! bulk processes
	do is=1,ntot ! loop over left site

		! choose the right site
		if (is < ntot) then
			is1 = is+1;
		elseif (periodic) then
			is1 = 1;
		else
			exit ! if not periodic then is=1,ntot-1
		endif
	
		call WhichBulkHop(is,is1,p,la,lo)
		!write(*,*) "=======>>> is, p, la,lo =",is, p, la,lo 
		if (p .gt. 0 ) then
			nps(p) = nps(p) + 1;
			act(p,nps(p)) = la
			oth(p,nps(p)) = lo
			if (p .eq. 7) then ! count both 7,8
				! la,lo = left,right sites for ih=7,8
				nps(8) = nps(8) + 1; 	!7: D,Phi; 8: Phi,D
				act(8,nps(8)) = lo 		! put d in active?
				oth(8,nps(8)) = la 		
			endif
		endif
		
	end do

	! contact processes: DO LATER

	!write(*,*) 'nps(:) = ',nps(:)
	!write(*,*) 'sys%occ = ',sys%occ

	! set gloable variable ways
	do p=1,8 ! only 8 at the moment, contact processes later
		ways(p)%ns = nps(p)
		if(allocated(ways(p)%active)) deallocate(ways(p)%active)
		if(allocated(ways(p)%sites)) deallocate(ways(p)%sites)
		if (nps(p) .gt. 0) then
			!write(*,*) "nps(p) > 0 ?????   p, nps(p) = ",p, nps(p)		
			allocate(ways(p)%active(nps(p)))
			allocate(ways(p)%sites(nps(p)))
			ways(p)%active = act(p,1:nps(p))
			ways(p)%sites = oth(p,1:nps(p))
			!write(*,*) 'ways(p)%active = ',ways(p)%active
			!write(*,*) 'ways(p)%sites = ',ways(p)%sites
		endif
	end do
	
	deallocate(act,oth)

	!-----------------------------------
	! allocate qt etc
	!-----------------------------------



	

	end subroutine UpdateWays
!**********************************************
!-----------------------	-----------------------	
! 	WhichBulkHop(l,r):
!	 find which hop is possible
!	 for two nearest neighbour sites
!	 with occupancies l,r for left and right sites
!-----------------------	-----------------------	
! list of hops:
! 1, D hops Right
! 2, D hops Left
! 3, Phi hops right
! 4, Phi hops left
! 5, D, Phi annihilated
! 6, Phi,D annihilated
! 7, D, Phi created
! 8, Phi,D created
	subroutine WhichBulkHop(is,is1,p,la,lo)
	implicit none
	integer(kind=1), intent(in) :: is,is1
	integer(kind=1), intent(out) :: p,la,lo 
	! process ind, active site, other site (or phi if D,Phi annihilation)
	! local
	integer:: l,r ! occupation on left, right sites

	l = sys%occ(is);
	r = sys%occ(is1)
	! both D or both Phi
	if ((l==0 .and. r==0).or.(l==2 .and. r==2))then
		p=-1; la=-1;lo=-1; ! no way
	return
	endif

	! both active
	if ((l==1 .and. r==1))then
		p = 7; la=is;lo=is1 ! 7,8 both possible: phi,d | d,phi created
	return
	endif

	! left active
	if (l==1) then
		la=is;lo=is1;
		if (r==0) then ! phi
			p = 4 ! phi hops left
		else ! D
			p = 2 ! D hops left
		endif
	return
	endif

	! right active
	if (r==1) then
		la=is1;lo=is;
		if (l==0) then ! phi
			p = 3 ! phi hops right
		else ! D
			p = 1 ! D hops right
		endif
	return
	endif

	! left D
	if (l==2) then
		! right must be a phi, since all other options crossed above.
		la=is;lo=is1; ! la=D
		p = 5 ! d,phi annihilated
	elseif (l==0) then 	! left Phi
		! right must be a D, since all other options crossed above.
		la=is1;lo=is; ! la=D
		p = 6 ! phi,d annihilated
		!write(*,*) "=============== la,lo = ",la,lo		

	else
		p = -1; la=-1;lo=-1;
		write(*,*) "modways: something wrong .... "
		stop
	endif


	return
	end subroutine WhichBulkHop
!**********************************************

	subroutine UpdateOcc(ih,is)
	! only for bulk hoppings at the moment
	use modmain, only: sys,Asites,ways
	use lists, only: Drop, Join, Substitute
	implicit none
	integer, intent(in):: ih,is
	! local
	integer(kind=1) :: la,lo,n0,n1,n2,two=2

	if (ih > 8 .and. ih < 25 ) then
		write(*,*) "UpdateOcc: ERROR!!!"
		write(*,*) "contact processes not done yet"
		stop
	endif

	la = ways(ih)%active(is)
	lo = ways(ih)%sites(is)

	n0 = sys%n0; ! no of phi site
	n1 = sys%n1; ! no of active site
	n2 = sys%n2; ! no of d site

	!write(*,*)"la,lo=",la,lo
	!write(*,*)"n0,n1,n2=",n0,n1,n2
	!write(*,*) "sys%occ 1=",sys%occ
		
	select case(ih)	
		case(1,2)
			if(sys%occ(la)==1 .and. sys%occ(lo)==2)then
				sys%occ(la) = 2
				sys%occ(lo) = 1
			else
				write(*,*)"Error(UpdateOcc): case 1,2"
				stop
			endif
			call Substitute(Asites,n1,la,lo)
		case(3,4)
			if(sys%occ(la)==1 .and. sys%occ(lo)==0)then
				sys%occ(la) = 0
				sys%occ(lo) = 1
			else
				write(*,*)"Error(UpdateOcc): case 3,4"
				stop
			endif
			call Substitute(Asites,n1,la,lo)
		case(5,6)
			if(sys%occ(la)==2 .and. sys%occ(lo)==0)then
				sys%occ(la) = 1
				sys%occ(lo) = 1
			else
				write(*,*)"Error(UpdateOcc): case 5,6"
				stop
			endif
			! append first D then Phi; la=D for ih=5,6
			call Join(Asites,n1,(/la,lo/),two,Asites)
			sys%n0 = n0-1
			sys%n1 = n1+2
			sys%n2 = n2-1
		case(7,8)
			sys%occ(la) = 2 ! la=D,lo=phi for ih=7,8
			sys%occ(lo) = 0
			call Drop(Asites,n1,la,Asites)
			call Drop(Asites,n1,lo,Asites)
			sys%n0 = n0+1
			sys%n1 = n1-2
			sys%n2 = n2+1
	end select		

	if(sys%n0+sys%n1+sys%n2 .ne. sys%nsites)then
		write(*,*)"modways: nsites != n0+n1+n2 "
		stop
	endif

	return
	end subroutine UpdateOcc
!**********************************************
	end module modways
