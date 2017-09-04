
	!	calls routines that calculate amplitudes for all hops
	module Hoppings
	use modmain, only: ways,PermSym,na,nx,crosshops,
     .              nokappa,nogamma
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

	write(*,*) "hops: ---------- 1"
	!-------------------------------
	! ih:1-4 DHopsR/L, PhiHopsR/L
	! ?? L/R also same if PermSym. D/phi also same, swaped 1,2 channels?
	do ih=1,4
		do is=1,ways(ih)%ns,1
			if (.not. crosshops) then
				call dhops1(ih,is) ! chan 1-2
			else
				call dhops2(ih,is) ! chan 1-4
			endif
			if(PermSym) exit; ! only a single site/case for each hop type
		end do
	end do
		write(*,*) "hops: ---------- 2"

	!-------------------------------
	! ih:5-6 (D,Phi), (Phi,D) annihilation
	if (ways(5)%ns + ways(6)%ns > 0) then
		if (.not. crosshops) then
			call DPhiAn1() ! 1-2
		else
			call DPhiAn2() ! 1-4
 		endif
 	endif
 		write(*,*) "hops: ---------- 3"

	!-------------------------------
	! ih=7,8 (D,Phi), (Phi,D) created
	ih=7;	! ih=8 has the same amplitudes,
				! chan 1,2 swapped, handled via maph, mapc

	do is=1,ways(ih)%ns,1
			if (nx .ge. 1) then
				call DPhiCreat1(is) ! 1-2
				write(*,*) "hops: ---------creat 12"

			endif
			if (crosshops) then
				if (nx .ge. 2) then
				call DPhiCreat3(is) ! 3
				write(*,*) "hops: ---------creat 3"
				
				endif

				call DPhiCreat4(is) ! 4
				write(*,*) "hops: ---------creat 4"

			endif
			if(PermSym) exit! only a single site/case for each hop type
	end do
		write(*,*) "hops: ---------- 4"

	!-------------------------------
	!-------------------------------
	! cavity and exciton losses
	!-------------------------------
	if (nx .gt. 0) then
		! Exciton non-radiative decay
		if( .not. nogamma) then
			do is=1,na
				call LossGamma(is);
				if(PermSym) exit; ! only a single site/case for each hop type
			end do
		endif
		! Cavity photon losses
		if( .not. nokappa) then
			call LossKappa() ! kappa
		endif
	endif
	!-------------------------------

	write(*,*) "hops: ---------- 5"

	return
	end subroutine AllHops
!----------------------------------------------------------


	end module Hoppings
