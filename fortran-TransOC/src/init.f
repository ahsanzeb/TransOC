
	module init
	use modmain
	use maps
	use modways
	implicit none
	contains
!-----------------------------------------
!	initialise some global variables
!-----------------------------------------
	subroutine initialise()
	implicit none
	! local
	integer :: i,ina, n0,n1,n2




	!-----------------------------------------
	allocate(basis(NBasisSets))
	allocate(eig(NHilbertSpaces))
	allocate(Hg(NHilbertSpaces))
	!------------------------------------------





	!------------------------------------------
	!	sys: nsites, occ; na, Asites
	!------------------------------------------
	sys%nsites = nsites
	allocate(sys%occ(nsites))

	n0 = 0;n1=0;n2=0;
	do i=1,nsites
		sys%occ(i) =  int(rand(0)*(3)) ! random 0,1,2; just to test
		if (sys%occ(i)==0) then
			n0 = n0 +1
		elseif(sys%occ(i)==1)then
			n1 = n1 +1
		elseif(sys%occ(i)==2)then
			n2 = n2 +1
		else
			write(*,*) "init: Error, sys%occ(i) != 0,1,2 "
			stop
		endif	
	enddo
	! set Asites
	allocate(Asites(n1))
	ina = 1;
	do i=1,nsites
		if (sys%occ(i)==1) then
			Asites(ina) = i;
			ina = ina + 1;
		endif
	enddo

	! set global var
	sys%n0=n0;
	sys%n1=n1;
	sys%n2=n2;

	na = n1;
	write(*,*)"init: na = ",na
	!-----------------------------------------
	! initialise ways, mapb, mapt
	!-----------------------------------------
	call UpdateWays()
	call UpdateReqType ! ReqType

	!write(*,*) "ReqType = ",mapt%req

	! does not exist, set to false, for use in mapt/mapb
	Hg(:)%xst = .false.
	basis(:)%xst = .false.
	
	call UpdateMapT ! mapt%map, mapt%cal

	!write(*,*) "mapt = ",mapt%map


	call UpdateGroupTB ! mapt%cal ===> ntb, grouptb
	call UpdateMapB

	!-----------------------------------------
	! calculate maphc, the map from hopping index to amplitude index,
	!	and channel to channel location, needed for channel 8 only,
	! but does not cost much to make map for all 26 x 4.
	call calmaphc()
	!-----------------------------------------


	!-----------------------------------------
	! set occupations
	






	!-----------------------------------------
	! sets hopping parameters for various hops, channels
	!-----------------------------------------
	!{th, tl, tlh, thl, JhR, JlR, JhL, JlL} = tpar;
	!th = tpar(1); tl = tpar(2);
	!tlh = tpar(3); thl = tpar(4);
	!JhR = tpar(5); JlR = tpar(6);
	!JhL = tpar(7); JlL = tpar(8);
	ts = reshape( (/
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . tl, th, tlh, thl,
     . tl, th, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . th, tl, tlh, thl,
     . JlR,0.0d0,0.0d0,0.0d0, JhR,0.0d0,0.0d0,0.0d0,
     . JlL,0.0d0,0.0d0,0.0d0, JhL,0.0d0,0.0d0,0.0d0,
     . JhR,0.0d0,0.0d0,0.0d0, JlR,0.0d0,0.0d0,0.0d0,
     . JhL,0.0d0,0.0d0,0.0d0, JlL,0.0d0,0.0d0,0.0d0,
     . JlR,0.0d0,0.0d0,0.0d0, JhR,0.0d0,0.0d0,0.0d0,
     . JhR,0.0d0,0.0d0,0.0d0, JlR,0.0d0,0.0d0,0.0d0,
     . JlL,0.0d0,0.0d0,0.0d0, JhL,0.0d0,0.0d0,0.0d0,
     . JhL,0.0d0,0.0d0,0.0d0, JlL,0.0d0,0.0d0,0.0d0,
     . kappa,0.0d0,0.0d0,0.0d0, gamma,0.0d0,0.0d0,0.0d0
     . /), (/ 26,4 /), order=(/2,1/) );
	!-----------------------------------------
	!Block injecton?? 
	!-----------------------------------------
	!ideally should also minimise computation of
	!corresponding transition amplitudes!
	!EBlock electron blocking layer; on left
	!HBlock hole blocking layer; on right
	! jump # 16,21 prob->0
	if(EBlock) then
		ts(16, 1) = 0.0d0;
		ts(21, 1) = 0.0d0;
	endif
	! jump # 10,18 prob->0
	if(HBlock) then
		ts(10, 1) = 0.0d0;
		ts(18, 1) = 0.0d0;
	endif
	!-----------------------------------------


	!-----------------------------------------
	! set dqc global vaiable 
	!-----------------------------------------
	call SetDQC()
	!-----------------------------------------





	! only bulk processes? nh=1-8
	!nh = 26;
	!if(nocontacts) then









	end subroutine initialise
!-------------------------------------------------------------
!	sets dqc array, energetic changes for various hops/channels
!-------------------------------------------------------------
	subroutine SetDQC() ! w0, Ebr, Ebl, Er 
  !local
	double precision, dimension(4):: dEbulk
	double precision, dimension(18):: dEcont
	integer, dimension(26):: signEr
	integer :: i,j
  	
	dEbulk = (/ 0.0d0, 0.0d0, -w0, w0 /)
	dEcont = (/
     .   -Ebr, -Ebr + w0, -Ebl, -Ebl + w0,
     .   Ebr - w0, Ebr, Ebl - w0, Ebl,
     .   Ebr, -Ebr + w0, Ebr - w0, -Ebr,
     .   Ebl, -Ebl + w0, Ebl - w0, -Ebl,
     .   -w0, -w0 /);
	signEr =(/
     .   1, -1, -1, 1, 1, -1, -1, 1,
     .   1, 1, -1, -1, -1, -1, 1, 1,
     .   -1, 1, -1, 1, 1, -1, 1, -1,
     .   0, 0 /);

		dqc(:,:) = 0.0d0
		do i=1,26
			do j=1,4
   			if (i .le. 8) then
					dqc(i,j) = dEbulk(j);
				else
					dqc(i,j) = dEcont(i - 8)
				endif
				dqc(i,j) = dqc(i,j) - signEr(i)*Er
			end do
		end do
	return
	end subroutine SetDQC
!-----------------------------------------

	end module init
