

	subroutine init()
	use modmain
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
     . JlR,0.0,0.0,0.0, JhR,0.0,0.0,0.0,
     . JlL,0.0,0.0,0.0, JhL,0.0,0.0,0.0,
     . JhR,0.0,0.0,0.0, JlR,0.0,0.0,0.0,
     . JhL,0.0,0.0,0.0, JlL,0.0,0.0,0.0,
     . JlR,0.0,0.0,0.0, JhR,0.0,0.0,0.0,
     . JhR,0.0,0.0,0.0, JlR,0.0,0.0,0.0,
     . JlL,0.0,0.0,0.0, JhL,0.0,0.0,0.0,
     . JhL,0.0,0.0,0.0, JlL,0.0,0.0,0.0,
     . kappa,0.0,0.0,0.0, gamma,0.0,0.0,0.0
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











	end subroutine init
	
