
	!	calls routines that calculate amplitudes for all hops
	module Hoppings
	use modmain, only: ways,PermSym,na,nx,crosshops,
     . nokappa,nogamma,leads,nsites
	use dhops, only: dhops1, dhops2
	use creation !, only: DPhiCreat1,DPhiCreat3,DPhiCreat4
	use annihilation, only: DPhiAn1,DPhiAn2
	use losses, only: losskappa, lossgamma
	use Contacts, only: CDAnnihil, CDCreat
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
	integer:: ih,ih1
	integer:: is,l
	integer, dimension(8) ::
     .                  ihdphops=(/1,2,3,4,27,28,29,30/)

	!-------------------------------
	! ih:1-4 DHopsR/L, PhiHopsR/L
	! ih=27-30  DHopsU/D, PhiHopsU/D (Up/Down 3d cases)
	! ?? L/R also same if PermSym. D/phi also same, swaped 1,2 channels?
	do ih1=1,8
		ih = ihdphops(ih1)
		do is=1,ways(ih)%ns
			if (.not. crosshops) then
				call dhops1(ih,is) ! chan 1-2
			else
				call dhops2(ih,is) ! chan 1-4
			endif
			if(PermSym) exit; ! only a single site/case for each hop type
		end do
	end do
		!write(*,*) "hops: ---------- 2"

	!-------------------------------
	! ih:5-6,31-32 (D,Phi), (Phi,D) annihilation
	if(ways(5)%ns+ways(6)%ns+ways(31)%ns+ways(32)%ns >0)then
		if (.not. crosshops) then
			call DPhiAn1() ! 1-2
		else
			call DPhiAn2() ! 1-4
 		endif
 	endif
 	!write(*,*) "hops: ---------- 3"

	!-------------------------------
	! ih=7,8 (D,Phi), (Phi,D) created
	ih=7;	! ih=8 has the same amplitudes,
				! chan 1,2 swapped, handled via maph, mapc
	do is=1,ways(ih)%ns,1
		if (nx .ge. 1) then
			call DPhiCreat1(is) ! 1-2
		endif
		if (crosshops) then
			if (nx .ge. 2) then
				call DPhiCreat3(is) ! 3                                        
				endif
				call DPhiCreat4(is) ! 4
		endif
		if(PermSym) exit! only a single site/case for each hop type
	end do
	!write(*,*) "hops: ---------- 4"
	!-------------------------------	
	! ih=33,34 for in-plane hops  <====> ih=7,8
	if(.not. PermSym) then
		ih = 33; ! ih=34 will be treated like ih=8 with ih=7
		do is=1,ways(ih)%ns,1
			if (nx .ge. 1) then
				call DPhiCreat1(is) ! 1-2
			endif
			if (crosshops) then
				if (nx .ge. 2) then
				call DPhiCreat3(is) ! 3
				endif
				call DPhiCreat4(is) ! 4
			endif
			!if(PermSym) exit! only a single site/case for each hop type
		end do
	endif
	!-------------------------------
	!-------------------------------
	! cavity and exciton losses
	!-------------------------------
	if (nx .gt. 0) then
		! Exciton non-radiative decay
		if( .not. nogamma) then
			do is=1,na
				!write(*,*)"hop: gamma; is=",is
				call LossGamma(is);
				if(PermSym) exit; ! only a single site/case for each hop type
			end do
		endif

		!write(*,*) "hops: ---------- kappa"
	
		! Cavity photon losses
		if( .not. nokappa) then
			call LossKappa() ! kappa
		endif
	endif
	!-------------------------------

	! Contacts:
	if (leads) then
		! annihilation at contacts
		if (sum(ways(15:16)%ns+ways(11:12)%ns+
     .   ways(13:14)%ns+ways(9:10)%ns) > 0) then
			call CDAnnihil()
		endif

		! creation of D/Phi at contacts
		if (sum(ways(21:24)%ns)>0) call CDCreat('l')
		if (sum(ways(17:20)%ns)>0) call CDCreat('r')
	endif

	return
	end subroutine AllHops
!----------------------------------------------------------


	end module Hoppings
