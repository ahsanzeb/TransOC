
	module init
	use modmain
	implicit none

	contains
!-----------------------------------------
!	initialise some global variables
!-----------------------------------------
	subroutine initialise()
	implicit none
	! local
	double precision::th, tl, tlh, thl, JhR, JlR, JhL, JlL
	logical:: EBL, HBL

	!-----------------------------------------
	! initialise mapb, mapt
	!-----------------------------------------
	allocate(mapb%map(5))
	allocate(mapb%cal(5))
	allocate(mapt%map(13))
	allocate(mapt%cal(13))
	mapb%nnu = 5; mapt%nnu = 13 ! calc all 
	mapb%map = (/ 1,2,3,4,5 /)
	mapb%cal = (/ 1,2,3,4,5 /)
	mapt%map = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	mapt%cal = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13 /)
	!-----------------------------------------


	!-----------------------------------------
	! calculate maphc, the map from hopping index to amplitude index,
	!	and channel to channel location, needed for channel 8 only,
	! but does not cost much to make map for all 26 x 4.
	call calmaphc()
	!-----------------------------------------



	!-----------------------------------------
	! sets hopping parameters for various hops, channels
	!-----------------------------------------
	!{th, tl, tlh, thl, JhR, JlR, JhL, JlL} = tpar;
	th = tpar(1); tl = tpar(2);
	tlh = tpar(3); thl = tpar(4);
	JhR = tpar(5); JlR = tpar(6);
	JhL = tpar(7); JlL = tpar(8);
	hpar = reshape( (/
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
	EBL= BlockInjection(1); ! electron blocking layer; on left
	HBL= BlockInjection(2); ! hole blocking layer; on right
	! jump # 16,21 prob->0
	if(EBL) then
		hpar(16, 1) = 0.0d0;
		hpar(21, 1) = 0.0d0;
	endif
	! jump # 10,18 prob->0
	if(HBL) then
		hpar(10, 1) = 0.0d0;
		hpar(18, 1) = 0.0d0;
	endif
	!-----------------------------------------


	!-----------------------------------------
	! set dqc global vaiable 
	!-----------------------------------------
	call SetDQC()
	!-----------------------------------------











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
