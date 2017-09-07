


	module readinput
	use modmain
	implicit none


	contains
!************************************************	
	subroutine input()
	implicit none
	! local
	character(40) :: block
	character(20):: Geometry
	integer :: iostat

	!--------------------------!
	!     default values       !
	!--------------------------!

	nsites = 5;
	nx = 1;
	niter = 20;
	NBasisSets = 5;
	NHilbertSpaces=13;
	smalln = 50;
	diagmaxitr = 150;
	sameg = .true.
	fixmap = .false.
	ndsec = 4;
	g = 0.3d0;
	th = 1.0d0; tl=1.0d0; tlh=0.005d0; thl=1.0d0;
	JhR=10*th; JlR=10*th; JhL=10*th; JlL=10*th;
	EBlock = .false.; HBlock = .false.;
	periodic = .true.; onlybulk = .false.;
	Ebr = 0.7d0; Ebl = 0.7d0; 
	Er = 1.0d0
	w0 = 2.0d0
	dw = 0.0d0
	detuning = .false.
	kappa = 0.005d0
	gamma = 0.005d0
	beta = 40.0d0
	crosshops = .true.
	AlwaysLP = .true.


	nokappa = .false.
	nogamma = .false.	

	!--------------------------!
	!     read from input.in   !
	!--------------------------!
	open(50,file='input.in',action='READ',
     . 					status='OLD',form='FORMATTED',iostat=iostat)
	if (iostat.ne.0) then
		write(*,*)
		write(*,'("Error(readinput): error opening input.in")')
		write(*,*)
		stop
	end if

10		continue

	read(50,*,end=30) block
	! check for a comment
	if ((scan(trim(block),'!').eq.1).or.
     .  (scan(trim(block),'#').eq.1)) goto 10

	select case(trim(block))

	case('Nsites')
		read(50,*,err=20) nsites

	case('Nexcitations')
		read(50,*,err=20) nx

	case('Niter')
		read(50,*,err=20) niter

	!case('SameLightMoleculeCoupling') !set always to true as diff g not implemented yet
	!	read(50,*,err=20) sameg

	case('DirectSolverSize')
		read(50,*,err=20) smalln

	case('IterSolverMaxIter')
		read(50,*,err=20) diagmaxitr

	case('NBasisSets')
		read(50,*,err=20) NBasisSets

	case('NHilbertSpaces')
		read(50,*,err=20) NHilbertSpaces

	case('NLowDegSectors') 
		read(50,*,err=20) ndsec
		! dynamics will be restricted to low lying degenerate
		! sectors only by energetic penalty, so no need to
		! calculate amplitudes for higher sectors

	case('DontReuseData')
		read(50,*,err=20) fixmap
		write(*,*) "DontReuseData = ",fixmap

	case('CavityMoleculeCoupling')
		read(50,*,err=20) g

	case('BulkHoppingParameters')
		read(50,*,err=20) th, tl, tlh, thl

	case('ContactsHoppingParameters')
		read(50,*,err=20) JhR, JlR, JhL, JlL

	case('BlockInjection')
		read(50,*,err=20) EBlock, HBlock

	case('Geometry')
		read(50,*,err=20) Geometry
		if (trim(Geometry) == 'periodic') then
			periodic = .true.
			onlybulk = .false.	
		elseif(trim(Geometry) == 'device') then
			periodic = .false.	
			onlybulk = .false.			
		elseif(trim(Geometry) == 'onlybulk') then
			periodic = .false.;
			onlybulk = .true.
		else
			write(*,'("Error(readinput): Geometry unknown... ")')
			stop
		endif

	case('ContactsBarriers')
		read(50,*,err=20) Ebr,Ebl

	case('EFieldNNSEnergy')
		read(50,*,err=20) Er

	case('ExcitonEnergy')
		read(50,*,err=20) w0

	case('Detuning')
		read(50,*,err=20) dw

	case('Kappa')
		read(50,*,err=20) kappa

	case('Gamma')
		read(50,*,err=20) gamma

	case('Beta')
		read(50,*,err=20) beta

	case('CrossHops')
		read(50,*,err=20) crosshops

	case('AlwaysLP')
		read(50,*,err=20) AlwaysLP

	case('')
		goto 10
	case default
	write(*,*)
	write(*,'("Error(readinput): invalid block: ",A)')trim(block)
	write(*,*)
	stop
	end select
	goto 10
20		continue
	write(*,*)
	write(*,'("Error(readinput): error reading from input.in")')
	write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
	write(*,'("Check input convention in manual")')
	write(*,*)
	stop
30			continue
	close(50)


	if(abs(dw) > 1.d-6) detuning=.true.

	if(kappa < 1.d-6) nokappa = .true.
	if(gamma < 1.d-6) nogamma = .true.
	if(AlwaysLP) PermSym = .true.; ! XXXXX add conditions on gi=g/ei=w0 

	return
	end subroutine
!************************************************
	end module readinput


