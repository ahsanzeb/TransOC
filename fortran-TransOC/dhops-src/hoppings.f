
	!	calls routines that calculate amplitudes for all hops
	module Hoppings
	use modmain, only: ways,PermSym,na,nx,crosshops,debug,
     . nokappa,nogamma,leads,nsites,ihdphops, qt,maph,impurity
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
	integer, dimension(8) :: dpih1 = (/1,2,3,4,27,28,29,30/);

	public :: AllHops0, AllHops1

	contains
	!-------------------------------------------------------
	! perform all available hoppings. 
	!	the routines called set the gloabl variables related to 
	!	respective quantum amplitudes
	!-------------------------------------------------------
	! call when PermSym	
	subroutine AllHops0()
	implicit none
	! local
	integer:: ih,ih1
	integer:: l,ii, ia1,ia2, is=1
	logical :: contannih, contcreat, notdone

	!-------------------------------
	! if PermSym then amp same for Dhops/Phihops/left/right/up/down
	! see maph for sets of hops types sharing amplitudes
	!notdone = .true.
	!if(sum(ways(1:4)%ns + ways(27:30)%ns) >0) then
	do ih1=1,8
		ih=dpih1(ih1);
		if(ways(ih)%ns > 0) then
			if (.not. crosshops) then
				call dhops1(ih,is) ! chan 1-2
			else
				call dhops2(ih,is) ! chan 1-4
			endif
			exit
		endif
	enddo
	!endif

	return
	end subroutine AllHops0
!======================================================================
	! call when no PermSym	
	subroutine AllHops1()
	implicit none
	! local
	integer:: ih,ih1
	integer:: is,l,ii, ia1,ia2
	logical :: contannih
	!-------------------------------
	! ih:1-4 DHopsR/L, PhiHopsR/L
	! ih=27-30  DHopsU/D, PhiHopsU/D (Up/Down 3d cases)

	do ii=1,2 ! DHops, PhiHops
			ihs = dpih(:,ii)	
		do ih1=1,4
			ih = ihs(ih1)
			do is=1,ways(ih)%ns
				if (.not. crosshops) then
					call dhops1(ih,is) ! chan 1-2
				else
					call dhops2(ih,is) ! chan 1-4
				endif
			enddo ! is
		enddo ! ih
	end do ! ii
	!-------------------------------

	return
	end subroutine AllHops1
!======================================================================

!---------------------------------
! check if amp for ih are already calcuated?
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
!-------------------------------------------
	end module Hoppings
