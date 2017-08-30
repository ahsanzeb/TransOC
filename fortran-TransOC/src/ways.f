	subroutine UpdateWays()
	use modmain, only: sys, ways

	! local
	integer:: ntot, n0,n1,n2,i,j
	integer, allocatable, dimension(:) :: occ
	integer, dimension(26) :: nps
	integer, dimension, dimension(:,:) :: waysl

	! initialise
	nps = 0
	waysl = 0
	
	! copy current sys
	ntot = sys%nsites
	allocate(occ(ntot))
	occ = sys%occ
	!n0 = sys%n0
	!n1 = sys%n1
	!n2 = sys%n2

	allocate(waysl(26,ntot))


	! bulk processes
	do is=1,ntot-1
		p = WhichBulkHop(occ(is),occ(is+1))
		if (p .ne. 0 ) then
			nps(p) = nps(p) + 1;
			waysl(p,nps(p)) = is
			if (p .eq. 7) then ! count both 7,8
				nps(p+1) = nps(p+1) + 1; ! p=8 
				waysl(p+1,nps(p+1)) = is
			endif
		endif
	end do

	! contact processes: DO LATER

	! set gloable variable ways
	do p=1,8 ! only 8 at the moment, contact processes later
		ways(p)%ns = nps(p)
		if(allocated(ways(p)%sites)) deallocate(ways(p)%sites)
		if (nps(p) > 0) then
			allocate(ways(p)%sites(nps(p)))
			ways(p)%sites = waysl(p,1:nps(p+1))
		endif
	end do
	
	
	deallocate(waysl)


	contains
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
	integer function WhichBulkHop(l,r)
	implicit none
	integer, intent(in) :: l,r ! occupation on left, right sites
	integer :: p
	
	if (l==0) then
		select case (r)
			case (0)
      		p = 0 ! nothing
			case(1)
				p = 3 ! phi hps right
			case(2)
				p = 6 ! phi,d annihilated
		end select
	elseif(l==1) then
		select case (r)
			case (0)
      		p = 4 ! phi hops left
			case(1)
				p = 7 ! 7,8 both possible: phi,d | d,phi created
			case(2)
				p = 2 ! d hops left
		end select
	elseif(l==2)
		select case (r)
			case (0)
      		p = 5 ! d,phi annihilated
			case(1)
				p = 1 ! d hops right
			case(2)
				p = 0 ! nothing
		end select
	endif
	WhichHop = p
	return
	end function WhichHop
!-----------------------	

	! DO IT LATER.... 
	integer function WhichLCHop(occ)
	implicit none
	integer :: occ,p ! occupation on left, right sites

		select case (occ)
			case (0)
      		p = 0 ! 
			case(1)
				p = 0 ! 
			case(2)
				p = 0 ! 
		end select
	endif
	end function WhichLCHop

	
	end subroutine UpdateWays
