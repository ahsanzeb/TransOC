

! call in the following order
!-1. update ways
! 0. update na, nx
!	1. call UpdateReqType ! ReqType
! 2. call UpdateMapT ! mapt%map, mapt%cal
! 3. call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
!		grouptb has itypes that would be calculated in mkHamilt()
! 4. call UpdateMapB

	module maps
	implicit none

	contains
!-------------------------------------------
! REQUIRED ITYPES
! updated na,nx needed. so update na,nx before calling this.
	subroutine UpdateReqType()
	use modmain, only: nx,mapt,crosshops,nogamma,nokappa,
     .    ways !dna,dnx,ibs,dns,itypes
	implicit none
	integer :: it
	logical :: ch,exst,logNmm1,logkg,ldpa,lap,lcda,lcAp

	ch = crosshops;
	logNmm1= (sum(ways(1:4)%ns) > 0 .and. ch);
	logkg= ((.not. nokappa) .or. (.not. nogamma));
	ldpa = sum(ways(5:6)%ns) > 0;
	lap = ways(7)%ns > 0;
	lcda=sum(ways(9:16)%ns) > 0;
	lcAp=sum(ways(17:24)%ns) > 0;

	!write(*,*)"ch",ch
	!write(*,*)"logNmm1",logNmm1,ways(1:4)%ns
	!write(*,*)"logkg",logkg
	!write(*,*)"ldpa",ldpa,ways(5:6)%ns
	!write(*,*)"lap",lap
	!write(*,*)"lcda",lcda
	!write(*,*)"lcAp",lcap
	!write(*,*)

	do it=1,13
		exst = .false.
		select case(it)
		case(1)
			exst = .true. ! always
		case(2)
			if((logNmm1 .or. logkg) .and. nx>0)exst=.true.
		case(3)
			if(logNmm1) exst = .true.
		case(4)
			if(ldpa)exst = .true.
		case(5,6)
			if(ldpa .and. ch)exst = .true.
		case(7)
			if(lap .and. nx>0)exst = .true.
		case(8)
			if(lap .and. ch .and. nx>1)exst = .true.
		case(9)
			if(lap .and. ch)exst = .true.
		case(10,11)
			if(lcda) exst = .true.
		case(12)
			if(lcAp .and. nx>0) exst = .true.
		case(13)
			if(lcAp) exst = .true.
		end select

		mapt%req(it) = exst
	
	end do

	return
	end subroutine UpdateReqType
!-------------------------------------------
	subroutine UpdateMapT()
	!	map for itype: ==> ib, 13 Hilbert spaces
	use modmain, only : dna,dnx,ibs,dns,mapt,
     .  Hg,NHilbertSpaces, na,nx
	use lists, only: FreeQ
	! local
	integer(kind=1), dimension(13):: map,notexist
	integer :: it,n,m,i,inxt
	integer(kind=1):: j,thrt=13
	logical :: found,exst


	!write(*,*)"ReqType=",mapt%req
	map = -1;
	notexist = -1;
	inxt = 0;
	do it=1,13
			if (.not. mapt%req(it)) cycle
			n = na + dna(it); ! n,m values for itype=it
			m = nx + dnx(it);
			! search for these n,m locations among NHilbertSpaces
			exst = .false.
			do i=1,NHilbertSpaces
				if (Hg(i)%xst) then
					if (Hg(i)%n == n .and. Hg(i)%m == m)then
						!il = i;
						map(it) = i;
						exst = .true.
						exit
					endif
				endif
			enddo
			if (.not. exst) then
				inxt=inxt+1;
				notexist(inxt) = it
			endif	
	enddo ! it

	!write(*,*)"mapt  1 = ",map
	!write(*,*)"notfound",notexist


	!---------------------------------------
	! where to put to be calculated results?
	do i = 1,inxt
		it = notexist(i); ! finding slot/space for this itype
		found = .false.
		do j=1,NHilbertSpaces ! search for free slot
			if(.not. Hg(j)%xst .and. FreeQ(map,thrt,j)) then
				!Hg(j)%xst = .true.
				map(it) = j;
				found = .true.
				exit
			endif
		end do

		!write(*,*)"mapt 2 = ",map

		
		! if not found, search for existing
		!			 slots that are not required
		if (.not. found) then
			do j=1,NHilbertSpaces ! search for free slot
				if(FreeQ(map,thrt,j)) then ! slot j not used for another itype
					map(it) = j;	
					found = .true.
					exit
				endif
			end do
		endif
		! still not found? Error!
		if (.not. found) then
			write(*,*)"Error(maps): slot not found for itype=",it
			stop
		endif
	end do ! it

	!write(*,*)"mapt 3 = ",map

	!---------------------------------------
	! set global variables
	mapt%map = map;
	mapt%cal = notexist;
	mapt%nnu = inxt;

	!write(*,*)"mapt%map",mapt%map
	!write(*,*)"mapt%req",mapt%req
	return
	end subroutine UpdateMapT
!------------------------------------
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
!-------------------------------------------











	subroutine UpdateMapB()
	!	map for ib: ==> 5 basis
	use modmain, only : dns,mapb,NBasisSets,basis,na
	use lists, only: FreeQ
	! local
	integer(kind=1), dimension(5):: map,notused
	integer :: jb,jb2,i,inu
	integer(kind=1):: j,n,five=5
	logical :: used, found

	map = -1 ;! set for error checking
	notused = -1; inu=0;		
	do jb=1,5	
		n = na + dns(jb);
		! search in exisiting record
		used=.false.
		do jb2=1,NBasisSets 
			if (n == basis(jb2)%n) then
				map(jb) = jb2;
				used=.true.
				exit
			endif
		enddo	
		! add to the list that is to be calcualted
		if(.not. used) then
			inu = inu + 1;
			notused(inu) = jb
		endif			
	end do

	!write(*,*)"mapb 1= ",map

	! where to put to be calculated results?
	do i=1,inu
		jb = notused(i);
		! try to find a free slot
		found = .false.
		do j=1,NBasisSets
			if(.not. basis(j)%xst.and.FreeQ(map,five,j)) then
				map(jb) = j;
				found =.true.
				exit
			endif
		end do

	!write(*,*)"mapb 2= ",map

		! overwrite basis in an exisiting slot
		if (.not. found) then
			do j=1,NBasisSets
				if(FreeQ(map,five,j)) then ! slot j not used for another itype
					map(jb) = j;	
					found = .true.
					exit
				endif
			end do
		endif
		! still not found? Error!
		if (.not. found) then
			write(*,*)"Error(maps): slot not found for ib=",jb
			stop
		endif
	end do

	!write(*,*)"mapb 3 = ",map

	! set global variables
	mapb%map = map;
	mapb%cal = notused;
	mapb%nnu = inu;

	!write(*,*)"mapb%map",mapb%map

	
	return
	end subroutine UpdateMapB
!--------------------------------------------

	subroutine calmaphc()
	use modmain, only: maph,mapc
	implicit none

	!integer(kind=4), dimension(26,4,2) :: maphc
	integer(kind=4) :: ic,ih

	! input: ih,ic
	!	output: ia,icl
	! index: ih,ic; 1 ==> ia, 2 ==> icl

	! ia: 1 - 14 in Mathematica notaitons:
	! RDAR,RDAL,RPhiAR,RPhiAL
	! RDPhiA, RDPhiAinv, [ih=8: ic swap 1,2]
	! RCDAl, RCDAh
	! RCAlR,RCAhR
	! RCAlL,RCAhL
	! Rkappa, Rgamma
	
	do ih=1,26
		do ic=1,4
			! icl=ic except for ih=8
			mapc(ih,ic) = ic
			if(ih<5) then
				maph(ih,ic) = ih
			elseif(ih<7) then
				maph(ih,ic) = 5   
			elseif(ih==7) then
				maph(ih,ic) = 6 	
			elseif(ih==8) then
				maph(ih,ic) = 6
				if(ic==1) then
					mapc(ih,ic) = 2
				elseif(ic==2) then
					mapc(ih,ic) = 1
				endif	
			elseif(ih==9 .or. ih==11 .or. ih==13 .or. ih==15) then
				maph(ih,ic) = 7
			elseif(ih==10 .or. ih==12 .or. ih==14 .or. ih==16) then
				maph(ih,ic) = 8
			elseif(ih==17 .or. ih==18) then
				maph(ih,ic) = 9
			elseif(ih==19 .or. ih==20) then
				maph(ih,ic) = 10
			elseif(ih==21 .or. ih==22) then
				maph(ih,ic) = 11
			elseif(ih==23 .or. ih==24) then
				maph(ih,ic) = 12
			elseif(ih==25) then
				maph(ih,ic) = 13
			elseif(ih==26) then
				maph(ih,ic) = 14
			else
				write(*,*) "ih > 26 ?????"
				stop
			endif
		enddo !ic
	enddo ! ih

	end subroutine calmaphc
!-------------------------------------------
	subroutine DONTUSEmaps()
	! this would also ask to calculated basis and Hg/eig that are not even required! very bad! just use for testing.
	use modmain, only: mapt,mapb
	implicit none
	
 	mapb%nnu = 5; mapt%nnu = 13 ! calc all 
	mapb%map = (/ 1,2,3,4,5 /)
	mapb%cal = (/ 1,2,3,4,5 /)
	mapt%map = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	mapt%cal = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	!------------------------------------------
	!	initialise mapt% ntb/GroupTB
	!------------------------------------------
	mapt%ntb = (/3,2,3,2,3 /); ! no of cases with N=[N-2,N-1,N,N+1,N+2]
	!mapt%grouptb ==> which itypes for these 3,2,3,2,3 
	mapt%grouptb= reshape( (/
     .   7,8,9,0,0,0,0,0,0,0,0,0,0,
     .   12,13,0,0,0,0,0,0,0,0,0,0,0,
     .   1,2,3,0,0,0,0,0,0,0,0,0,0,
     .   10,11,0,0,0,0,0,0,0,0,0,0,0,
     .   4,5,6,0,0,0,0,0,0,0,0,0,0
     .   /),(/5,13/), order=(/2,1/))

	return
	end subroutine

!------------------------------------------
	end module maps
	
