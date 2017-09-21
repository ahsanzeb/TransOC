


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
	integer :: iostat,i,j
	logical :: givenNel,givenNelrange,givenDwRange,givenNexcit

	givenNel = .false.; givenNelrange=.false.
	givenDwRange=.false.

	!--------------------------!
	!     default values       !
	!--------------------------!

	debug = .false.
	nsites = 5;
	nel = 1;
	nelmin=1; nelmax=1;dnelec=1
	nx = 1;
	niter = 500;
	ntraj = 10;
	NBasisSets = 5;
	NHilbertSpaces=13;
	smalln = 50;
	nog = .false.
	diagmaxitr = 150;
	sameg = .true.
	fixmap = .false.
	ndsec = 3;
	g = 0.3d0;
	th = 1.0d0; tl=1.0d0; tlh=0.005d0; thl=1.0d0;
	JhR=10*th; JlR=10*th; JhL=10*th; JlL=10*th;
	EBlock = .false.; HBlock = .false.;
	periodic = .true.; onlybulk = .false.;
	leads = .false.;
	Ebr = 0.7d0; Ebl = 0.7d0; 
	Er = 1.0d0
	w0 = 2.0d0
	dw = 0.0d0
	dwmin=0.0; dwmax=0.0; ddw=0.0
	detuning = .false.
	kappa = 0.1 !0.005d0 ! ref to th=1.0d0; scale times proportionally tp th.
	gamma = 0.1 !0.005d0
	beta = 40.0d0
	crosshops = .true.
	AlwaysLP = .true.
	PermSym = .false.
	nokappa = .false.
	nogamma = .false.	
	nolosses = .false.
	mincarriers = .false.
	givenNexcit = .false.;
	ztout= .true.;
	ntrapout= .true.;
	ratesout= .true.;
	!std= .true.;
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

	case('debug')
		read(50,*,err=20) debug

	case('ztout')
		read(50,*,err=20) ztout

	case('ntrapout')
		read(50,*,err=20)ntrapout

	case('ratesout')
		read(50,*,err=20) ratesout

	!case('stdout')
	!	read(50,*,err=20) std

	case('Nsites')
		read(50,*,err=20) nsites

	case('Doping','doping')
		read(50,*,err=20) nel
		givenNel = .true.

	case('DopingRange','Dopingrange','dopingrange')
		read(50,*,err=20) nelmin, nelmax,dnelec
		givenNelrange = .true.;

	case('NExcitationsRange','NexcitationsRange','nexcitationsrange')
		read(50,*,err=20) mexmin,mexmax,dmex
		givenNexcit = .true.;


	case('DetuningRange','Detuningrange','detuningrange')
		read(50,*,err=20) dwmin,dwmax,ddw
		givenDwRange = .true.;
		
	case('Nexcitations')
		read(50,*,err=20) nx

	case('Niter')
		read(50,*,err=20) niter

	case('Ntraj','ntraj')
		read(50,*,err=20) ntraj

	!case('SameLightMoleculeCoupling') !set always to true as diff g not implemented yet
	!	read(50,*,err=20) sameg

	case('DirectSolverSize','smalln')
		read(50,*,err=20) smalln

	case('IterSolverMaxIter','diagmaxitr')
		read(50,*,err=20) diagmaxitr

	case('NBasisSets')
		read(50,*,err=20) NBasisSets

	case('NHilbertSpaces')
		read(50,*,err=20) NHilbertSpaces

	case('NLowDegSectors','ndsec') 
		read(50,*,err=20) ndsec
		! dynamics will be restricted to low lying degenerate
		! sectors only by energetic penalty, so no need to
		! calculate amplitudes for higher sectors

	case('DontReuseData')
		read(50,*,err=20) fixmap
		write(*,*) "DontReuseData = ",fixmap

	case('CavityMoleculeCoupling','g')
		read(50,*,err=20) g

	case('NoCoupling','nog')
		read(50,*,err=20) nog

	case('MinCarriers','Mincarriers','mincarriers')
		read(50,*,err=20) mincarriers

	case('BulkHoppingParameters')
		read(50,*,err=20) th, tl, tlh, thl

	case('ContactsHoppingParameters')
		read(50,*,err=20) JhR, JlR, JhL, JlL

	case('BlockInjection')
		read(50,*,err=20) EBlock, HBlock

	case('Geometry','geometry')
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

	case('ContactsBarriers','Barriers','barriers')
		read(50,*,err=20) Ebr,Ebl

	case('EFieldNNSEnergy','Er','er')
		read(50,*,err=20) Er

	case('ExcitonEnergy','w0')
		read(50,*,err=20) w0

	case('Detuning','detuning','dw')
		read(50,*,err=20) dw

	case('Kappa','kappa')
		read(50,*,err=20) kappa

	case('Gamma','gamma')
		read(50,*,err=20) gamma

	case('NoLosses','Nolosses','nolosses')
		read(50,*,err=20) nolosses

	case('Beta','beta')
		read(50,*,err=20) beta

	case('CrossHops','Crosshops','crosshops')
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

	!write(*,*)"1. kappa,gamma; nokappa,nogamma",
  !   . kappa,gamma, nokappa,nogamma

	!write(*,*)"1. nolosses = ",nolosses




	! contact/leads present?
	leads = (.not. periodic) .and. (.not. onlybulk);
	!if (sum(abs(JhR+JlR+JhL+JlL)) < 1.0d-5) then
	!	leads = .false.
	!	onlybulk = .true.
	!endif




	if(.not. givenNexcit) then
		mexmin=nx;
		mexmax=nx;
		dmex=1;
	endif


	!--------------------------
	! keep the order 
	!--------------------------
	if(abs(g)<1.0d-2) nog = .true.
	if (nog) then
		AlwaysLP = .false. ! stay in the state (= ksub set) after a transition
		sameg = .true. ! g=0 for all molecules
	endif
	!--------------------------


	!--------------------------
	! keep the order 
	!--------------------------
	! if asked in input
	if(nolosses) then 
		nokappa = .true.;
		nogamma = .true.;
	endif
	! if kappa/gamma rates are too small
	if(kappa < 1.d-6) nokappa = .true.
	if(gamma < 1.d-6) nogamma = .true.
	if(nokappa .and. nogamma) then 
		nolosses = .true.
	else ! 
		nolosses = .false.
	endif
	!--------------------------


	!--------------------------
	! detuning?
	!--------------------------
	if(abs(dw) > 1.d-6 .and. (.not. givenDwRange)) then
		detuning=.true.
		dwmin=dw
		dwmax=dw
		ddw=dw
		ndw=1
	elseif(givenDwRange) then
		detuning=.true.
		if(abs(dwmax-dwmin) < 1.d-6)then
			write(*,*)"Warning(readinput): dwmax-dwmin too small!"
		endif
		ndw=int((dwmax-dwmin)/ddw) + 1;
	else ! abs(dw) < 1.d-6, range not given
		detuning=.false.
		dwmin=dw !dw=0
		dwmax=dw
		ddw=dw
		ndw=1
	endif
	!--------------------------


	
	if(AlwaysLP)then
		PermSym = .true.; ! XXXXX add conditions on gi=g/ei=w0 
	else
		PermSym = .false.; 
	endif
	
	if(thl < 1.d-6 .and. tlh < 1.d-6) then
		crosshops = .false.
		write(*,*) "main: tlh and thl are too small so no crosshops"
	endif

	if (.not. givenNelrange) then
		nelmin=nel; nelmax=nel; dnelec=1;
	endif

	if(nelmin <= -nsites .or. nelmax >= nsites) then
		nelmin = max(nelmin, -(nsites-1) );
		nelmax = min(nelmax, nsites-1);
		!dnelec=1;
		write(*,*)"Warning(readinput): DopingRange invalid!"
		write(*,*)" changed to sensible min, max values."
	endif




	!---------------------------------------------
	! write various parameters
	i = 1+int((mexmax-mexmin)/dmex);
	j	=	1+int((nelmax-nelmin)/dnelec)
	if (onlydoped) j=1+int((nelmax-nelmin-1)/dnelec);
	open(200,file='parameters.out',action='write')
		write(200,*) i,ndw,j,ntraj,niter
		write(200,*) mexmin,mexmax,dmex
		write(200,*) dwmin,dwmax, ddw
		write(200,*) nelmin,nelmax,dnelec
		write(200,*) ntraj, niter
	close(200)
	!---------------------------------------------


	return
	end subroutine
!************************************************
	end module readinput


