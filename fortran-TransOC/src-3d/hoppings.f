
	!	calls routines that calculate amplitudes for all hops
	module Hoppings
	use modmain, only: ways,PermSym,na,nx,crosshops,
     . nokappa,nogamma,leads,nsites,ihdphops, qt,maph
	use dhops, only: dhops1, dhops2
	use creation !, only: DPhiCreat1,DPhiCreat3,DPhiCreat4
	use annihilation, only: DPhiAn1,DPhiAn2
	use losses, only: losskappa, lossgamma
	use Contacts, only: CDAnnihil, CDCreat
	implicit none
	integer, dimension(4,2) ::
     . dpih = reshape( (/1,2,27,28, 3,4,29,30/),
     .            (/ 4,2 /), order=(/1,2/) );
	integer, dimension(4) :: ihs
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
	integer:: is,l,ii, ia1,ia2
	logical :: found
	!-------------------------------
	! ih:1-4 DHopsR/L, PhiHopsR/L
	! ih=27-30  DHopsU/D, PhiHopsU/D (Up/Down 3d cases)
	! ?? L/R also same if PermSym. D/phi also same, swaped 1,2 channels?

	do ii=1,2 ! DHops, PhiHops
			ihs = dpih(:,ii)	
		do ih1=1,4
			ih = ihs(ih1)
			do is=1,ways(ih)%ns
				found = .false.
				if(PermSym) then
					call CheckAmp(ihs,ih,found,ia1,ia2)		
					if(found) then ! copy qt instead of repeat calc
						call copyqt(ia1,ia2); ! takes care of all channels 1-2/1-4	
					endif
				endif
				if(.not. found) then ! due to PermSym=F or from CheckAmp()
					if (.not. crosshops) then
						call dhops1(ih,is) ! chan 1-2
					else
						call dhops2(ih,is) ! chan 1-4
					endif
				endif
				if(PermSym) exit; ! only a single site/case for each hop type
			enddo ! is
		enddo ! ih
	end do ! ii

	!-------------------------------
	! ih:5-6,31-32 (D,Phi), (Phi,D) annihilation
	if(ways(5)%ns+ways(6)%ns+ways(31)%ns+ways(32)%ns >0)then
		! use ih=5 in DPhiAn1()/DPhiAn2() to find ia in Amp/Amp0
		! ih=5,6,31,32 all have same ampliudes/ia
		if (.not. crosshops) then
			call DPhiAn1() ! 1-2
		else
			call DPhiAn2() ! 1-4
 		endif
 	endif
 	!write(*,*) "hops: ---------- 3"

	!-------------------------------
	! ih=7,8 (D,Phi), (Phi,D) created
	!ih=7;	! ih=8 has the same amplitudes,
				! chan 1,2 swapped, handled via maph, mapc
	ihs(1:2) = (/7,33/);
	do ih1=1,2
	ih = ihs(ih);
	do is=1,ways(ih)%ns,1
		found = .false.
		if(PermSym) then
			if(ih==33 .and. ways(7)%ns > 0) then ! copy qt of ih=7 to ih=33
				ia1 = maph(7,1); ia2 = maph(33,1);
				call copyqt(ia1,ia2);
				found = .true.
			endif
		endif
		if(.not. found) then ! calculate 
			if (nx .ge. 1) then
				call DPhiCreat1(ih,is) ! 1-2
			endif
			if (crosshops) then
				if (nx .ge. 2) then
					call DPhiCreat3(ih,is) ! 3                                        
				endif
				call DPhiCreat4(ih,is) ! 4
			endif
		endif
		if(PermSym) exit! only a single site/case for each hop type
	enddo
	enddo

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

		! creation of D/Phi at contacts;
		do is=1,ways(21)%ns ! ih=21:24 same amplitudes
				call CDCreat('l',is)
				if(PermSym) exit;
		end do
			
		do is=1,ways(17)%ns ! ih=17:20 same amplitudes
			call CDCreat('r',is)
			if(PermSym) exit;
		end do
	endif
	!-------------------------------

	return
	end subroutine AllHops
!---------------------------------

	subroutine CheckAmp(a,ih,found,ia1,ia2)
	implicit none
	integer, dimension(4), intent(in):: a
	integer, intent(in):: ih
	integer, intent(out):: ia1, ia2
	logical, intent(out):: found
	integer:: i
	
	ia1=0; ia2=0;
	found=.false.;
	do i=1,4
		if (a(i) < ih .and. ways(a(i))%ns > 0) then
			found=.true.;
			ia1 = maph(a(i),1);
			ia2 = maph(ih,1); ! icl=ic: only for D,Phi hops, R,L,U,D;
			exit
		endif
	end do

	return
	end 	subroutine CheckAmp
!---------------------------------
	subroutine copyqt(ia1,ia2)
	implicit none
	integer, intent(in) :: ia1, ia2
	integer:: ic, nc, n2
	integer :: is=1

		! some conditions on nx etc might stop calc of ia1 amp
		! should skip such cases; but since rates() does not use qt in such cases,
		! we can just copy incorrect data from prev iterations to avoid checking all
		! the conditions here.

	if(crosshops) then
		nc=4
	else
		nc =2
	endif
	
	do ic=1,nc
		! --------- allocate amp-----------
		! find the size of ia1 array amp
		n2 = qt(ia1)%cs(ic,is)%namp
		if(n2>0) then
			! allocate ia2 array amp
			if(allocated(qt(ia2)%cs(ic,is)%amp))
     .						deallocate(qt(ia2)%cs(ic,is)%amp)
			allocate(qt(ia2)%cs(ic,is)%amp(n2))
			qt(ia2)%cs(ic,is)%namp = n2
			! --------- allocate amp2 -----------
			! find the size of ia1 array amp2
			n2 = qt(ia1)%cs(ic,is)%nsec
			! 	allocate ia2 array amp2
			if(allocated(qt(ia2)%cs(ic,is)%amp2))
     .					deallocate(qt(ia2)%cs(ic,is)%amp2)
			allocate(qt(ia2)%cs(ic,is)%amp2(n2))
			qt(ia2)%cs(ic,is)%nsec = n2
			! --------------assign --------------
			qt(ia2)%cs(ic,is)%amp = qt(ia1)%cs(ic,is)%amp
			qt(ia2)%cs(ic,is)%amp2 = qt(ia1)%cs(ic,is)%amp2
			!------------------------------------
		endif
	enddo

	return
	end subroutine copyqt
!---------------------------------


	logical function AmpExist(a,x)
	implicit none
	integer, dimension(4), intent(in):: a
	integer, intent(in):: x
	integer:: i, ns

	ns = 0;
	do i=1,4
		if (a(i) < x) then
			ns = ns + ways(a(i))%ns
		else
			exit
		endif
	end do
	AmpExist = (ns > 0);
	
	return
	end function
!-------------------------------------------
	end module Hoppings
