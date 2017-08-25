
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
	integer:: ns

	
	!use ways, only: 
	! NOTE: 	hop(ih)%chan(:)%site(:) dimensions of chan/sites allocated elsewhere 
		ns = hop(1)%chan(:)%ns
		allocate(hop(1)%chan(1)%site(ns))
		allocate(hop(1)%chan(2)%site(ns))


	end subroutine dhops




















	end module hoppings

