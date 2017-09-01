
	!	calls routines that calculate amplitudes for all hops
	module Hoppings
	use modmain, only: ways,PermSym,na,nx,crosshops
	use dhops, only: dhops1, dhops2
	use creation !, only: DPhiCreat1,DPhiCreat3,DPhiCreat4
	use annihilation, only: DPhiAn1,DPhiAn2
	use losses, only: losskappa, lossgamma

	implicit none

	public :: AllHops

	contains
	!-------------------------------------------------------
	! perform all available hoppings. 
	!	the routines called set the gloabl variables related to 
	!	respective quantum amplitudes
	!-------------------------------------------------------	
	subroutine AllHops()

	implicit none
	! local
	integer(kind=1) :: ih,is,l

	!p=1,8 bulk processes
	!ways(p)%ns
	!ways(p)%sites ! ???? correct this array to store something meaningful
	!ways(ih)%active(is); ! active site 
			

	!-------------------------------
	! ih:1-4 DHopsR/L, PhiHopsR/L
	! ?? L/R also same if PermSym. D/phi also same, swaped 1,2 channels?
	do ih=1,4
		do is=1,ways(ih)%ns,1
			l = ways(ih)%active(is); ! active site 
			if (.not. crosshops) then
				call dhops1(ih,l) ! chan 1-2
			else
				call dhops2(ih,l) ! chan 1-4
			endif
			if(PermSym) exit; ! only a single site/case for each hop type
		end do
	end do
	!-------------------------------
	! ih:5-6 (D,Phi), (Phi,D) annihilation
	if (ways(5)%ns + ways(6)%ns > 0) then
		if (.not. crosshops) then
			call DPhiAn1() ! 1-2
		else
			call DPhiAn2() ! 1-4
 		endif
 	endif
	!-------------------------------
	! ih=7,8 (D,Phi), (Phi,D) created
	ih=7;	! ih=8 has the same amplitudes,
				! chan 1,2 swapped, handled via maph, mapc
	do is=1,ways(ih)%ns,1
			l = ways(ih)%active(is); ! active site, on the left in the lattice
			if (nx .ge. 1) then
				call DPhiCreat1(is,l) ! 1-2
			endif
			if (crosshops) then
				if (nx .ge. 2) then
				call DPhiCreat3(is,l) ! 3
				endif
				call DPhiCreat4(is,l) ! 4
			endif
			if(PermSym) exit! only a single site/case for each hop type
	end do
	!-------------------------------

	!-------------------------------
	! cavity and exciton losses
	!-------------------------------
	if (nx .gt. 0) then
		! Exciton non-radiative decay
		do is=1,na
			call LossGamma(is);
			if(PermSym) exit; ! only a single site/case for each hop type
		end do
		! Cavity photon losses
		call LossKappa() ! kappa
	endif
	!-------------------------------

	return
	end subroutine AllHops
!----------------------------------------------------------


	end module Hoppings
