

	subroutine UpdateMapB(wj,wc)
	!	map for ib: ==> 5 basis
	use modmain, only : dna,dnx,ibs,dns,itypes,mapb
	integer, intent(in):: wj,wc
	! local
	integer(kind=1), dimension(5):: temp,map,notused
	integer :: itype, ib,d,jb,jb2,i,inu
	logical :: used

		map = mapb%map;
		
		itype = itypes(wj,wc);
		ib = ibs(itype);
		d = dns(ib);
		temp = map;
		map = -1 ;! set for error checking
		notused = -1;
		inu=1;
		do jb=1,5
			used=.false.
			do jb2=1,5
				if (jb2==jb-d) then
					map(jb2) = temp(jb);
					used=.true.
					exit
				endif
			enddo			
			if(.not. used) then
				notused(inu) = jb
				inu = inu + 1;
			endif			
		end do
		inu = inu - 1
		
		write(*,*) map
		! where to put to be calculated results?

		do i = 1,inu,1		
			do jb=1,5
				if (map(jb)==-1) then
						map(jb) = temp(notused(i))
						exit ! exit jb loop
				endif
			end do
		end do

	! set global variables
	mapb%map = map;
	mapb%cal = notused;
	mapb%nnu = inu;
	
	return
	end subroutine UpdateMapB

	subroutine UpdateMapT(wj,wc)
	!	map for itype: ==> ib, 13 Hilbert spaces
	use modmain, only : dna,dnx,ibs,dns,itypes,mapt
	integer, intent(in):: wj,wc
	! local
	integer(kind=1), dimension(13):: temp,map,notused
	integer :: itype, it1,it2,dn,dm,d,dnsel,dmsel
	integer :: i,inu
	logical :: used

		map = mapt%map;	
				
		itype = itypes(wj,wc);
		dnsel = dna(itype)
		dmsel = dnx(itype)

		! copy current map
		temp = map;
		! set to -1 etc for error check
		map = -1; notused = -1;
		inu = 1
		! calc new map
		do it1=1,13
			! next iteration n,m
			dn = dna(it1) - dnsel ! dnselected
			dm = dnx(it1) - dmsel;
			used=.false.
			do it2=1,13
				if ( dn==dna(it2) .and. dm==dnx(it2)) then
					map(it2) = temp(it1);
					used=.true.
					exit ! exit it2 loop
				endif
			end do
			if(.not. used) then
				notused(inu) = it1
				inu = inu + 1;
			endif
		end do
		inu = inu-1
		
		write(*,*) map
		! where to put to be calculated results?
		do i = 1,inu,1		
			do it1=1,13
				if (map(it1)==-1) then
						map(it1) =  temp(notused(i));
						write(*,*) "Calc new Hg etc for itype ",it1
						exit
				endif
			end do
		end do

!		do it1=1,13
!			if (map(it1)==-1) then
!				write(*,*) "Calc new Hg etc for itype ",it1
!			end if
!		end do
		
	! set global variables
	mapt%map = map;
	mapt%cal = notused;
	mapt%nnu = inu;

	return
	end subroutine UpdateMapT




	subroutine UpdateGroupTB()
	!	group itypes with the same ib
	use modmain, only : mapt,ibs
	! local
	integer(kind=1), dimension(13):: cal
	integer(kind=1), dimension(5,13):: grouptb
	integer(kind=1), dimension(5):: ntb
	
	integer :: ib,i,nnu,itype


	cal = mapt%cal
	nnu = mapt%nnu

	ntb=0; grouptb=0
	do i=1,nnu
		itype = cal(i);
		ib = ibs(itype);
		ntb(ib) = ntb(ib) + 1; ! number of cases for ib
		grouptb(ib,ntb(ib)) = itype; ! itype for ib
	end do

	! set global arrays
	mapt%ntb = ntb
	mapt%grouptb = grouptb
	
	return
	end subroutine UpdateGroupTB

