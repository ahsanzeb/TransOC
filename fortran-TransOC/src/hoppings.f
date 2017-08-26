
	module hoppings
	! calculates tranistion matrices for all allowed hops
	implicit none

	contains
!------------------------------------------

	! loops over various types of hops
	!	avoid calc matrices calculated in prev iteration
	!	that can be used in this iteration.... ?????

	subroutine dhops()
	use modmain, only: hop
	implicit none
	integer(kind=1):: ns,nc,ih,i,is

	
	!use ways, only: 
	! NOTE: 	hop(ih)%chan(:)%site(:) dimensions of chan/sites allocated elsewhere 


		ih = 1;
		hop(ih)%sites = (/ 1,2,3 /)
		ns = size(hop(ih)%sites);	
		
		nc=4; 
		hop(ih)%nc = nc
		hop(ih)%ns = ns
		allocate(hop(ih)%ht(nc,ns))
 
		!do i=1,ns
		!	is = hop(ih)%sites(i);
		!	call dhops1(ih,is)
		!end do

		!write(*,*) " dhops1 done.... "

		do i=1,ns
			is = hop(ih)%sites(i);
			call dhops2(ih,is)
		end do

		write(*,*) " dhops2 done.... "



		
	end subroutine dhops




















	end module hoppings

