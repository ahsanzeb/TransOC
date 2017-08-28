

	!module maps
	program maps
	implicit none

	integer(kind=1), dimension(26,4):: itypes
	integer(kind=1), dimension(13):: dna,dnx,ibs
	integer(kind=1), dimension(5):: dns
	integer :: wj,wc,i
	integer(kind=1), dimension(5):: mapb
	integer(kind=1), dimension(13):: mapt,num
	
	dna	= (/ 0,0,0,2,2,2,-2,-2,-2,1,1,-1,-1 /);
	dnx = (/ 0,-1,1,1,0,2,-1,-2,0,1,0,-1,0 /);
	ibs = (/ 3,3,3,5,5,5,1,1,1,4,4,2,2 /);
	dns = (/ -2,-1,0,1,2 /); 
	itypes = reshape( ( / 1, 1, 2, 3, 1, 1, 2, 3, 1, 1,
     . 2, 3, 1, 1, 2, 3, 4, 4, 5, 6,
     . 4, 4, 5, 6, 7, 7, 8, 9, 7, 7,
     . 8, 9, 11, 11, 11, 11, 10, 10, 
     . 10, 10, 11, 11, 11, 11, 10, 10,
     . 10, 10, 11, 11, 11, 11, 10,
     . 10, 10, 10, 11, 11, 11, 11,
     . 10, 10, 10, 10, 13, 13, 13,
     . 13, 13, 13, 13, 13, 12, 12,
     . 12, 12, 12, 12, 12, 12, 13,
     . 13, 13, 13, 13, 13, 13, 13,
     . 12, 12, 12, 12, 12, 12, 12,
     . 12, 2, 2, 2, 2, 2, 2, 2, 2 /),
     . (/ 26,4 /), order=(/2,1/) );     


	num = (/ (i,i=1,13,1)/)

	!do wj =1,26
	!	write(*,*) "wj = ",wj," itype = ", itypes(wj,:)
	!end do



	! start with some simple values
	mapb = (/ (i,i=1,5,1) /)
	mapt = (/ (i,i=1,13,1) /)

	write(*,*) "mapb; mapt: "
	write(*,*) mapb
	write(*,*) mapt


	do i=1,10
		! select random wj,wc
		wj=int(rand(0)*(26))+1
		wc=int(rand(0)*(4))+1
		write(*,*) " ----- iteration ------- ",i
		write(*,*) "wj,wc = ",wj,wc
		! change maps accordingly
		call GetMapIb(wj,wc,mapb)
		call GetMapItype(wj,wc,mapt)
		write(*,*) mapb
		write(*,*) num
		write(*,*) mapt
	end do



	
	contains

	subroutine GetMapIb(wj,wc,map)
	!	map for ib: ==> 5 basis
	integer, intent(in):: wj,wc
	integer(kind=1), dimension(5),intent(inout):: map
	integer(kind=1), dimension(5):: temp,notused
	integer :: itype, ib,d,jb,jb2,i,inu
	logical :: used
	
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

		write(*,*) map
		! where to put to be calculated results?
		do i = 1,inu-1,1		
			do jb=1,5
				if (map(jb)==-1) then
						map(jb) = temp(notused(i))
						exit ! exit jb loop
				endif
			end do
		end do
		
	return
	end subroutine GetMapIb

	subroutine GetMapItype(wj,wc,map)
	!	map for itype: ==> ib, 13 Hilbert spaces
	integer, intent(in):: wj,wc
	integer(kind=1), dimension(13),intent(inout):: map
	integer(kind=1), dimension(13):: temp,notused
	integer :: itype, it1,it2,dn,dm,d,dnsel,dmsel
	integer :: i,iu,inu
	logical :: used

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

		write(*,*) map
		! where to put to be calculated results?
		do i = 1,inu-1,1		
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
		

	return
	end subroutine GetMapItype




	end program
	!end module maps
