


	module readinput
	use modmain
	implicit none


	contains
!************************************************	
	subroutine input(node)
	implicit none
	integer, intent(in):: node
	! local
	character(40) :: block
	character(20):: Geometry
	integer :: iostat,i,j,idv
	logical :: givenNel,givenNelrange,givenDw,givenNexcit,givenEr
	logical :: OneDchains, giveExb
	integer:: imptype=0, lossgain=-1;
	double precision :: Eimp0=0.0d0
	logical :: givenT, giveng !, givenEblr

	
	giveExb = .false.
	givenNel = .true.; 
	givenNelrange=.false.
	givenDw=.false.
	givenEr = .false.;
	givenT =.false.; giveng=.false.; !givenEblr=.false.;
	OneDchains = .false.;
	!--------------------------!
	!     default values       !
	!--------------------------!
	epsr = 3.0d0; ! relative permitivity/dielectric constant
	a0 = 5.0d0; sigma0 = 0.25d0; nsigma=1.0; dinvl= 1.0d0; 
	! units: a0(nm), sigma0(a0), nsigma(sigma0), dinvl(1/a0)
	! dinvl= inverse localisation lenght:
	!	controls how fast hopping integeral decays; assumed tvrh=1.0*t_{h,l,hl,..}, etc, at dnns = a0;
	nproc = 34; ! number of hopping processes, useful for testing 1D chains by removing Up/down hops
	EqualDistr = .false.;! 1D chains: distribute carriers equally among chains
	debug = .false.
	nsites = 6;
	!J0=1.0d0;
	wcut=0.3d0;
	nel = -1;
	nelmin=1; nelmax=1;dnelec=1;
	Eq = 0.0;
	nx = 1;
	niter = 10;
	ntraj = 1;
	NBasisSets = 5;
	NHilbertSpaces=13;
	smalln = 50;
	nog = .false.
	diagmaxitr = 150;
	sameg = .true.
	fixmap = .false.
	ndsec = 3;
	g = 0.3d0;
	th = 1.0d0; tl=1.0d0; tlh=1.0d0; thl=1.0d0;
	JhR=10.0d0; JlR=10.0d0; JhL=10.0d0; JlL=10.0d0;
	EBlock = .false.; HBlock = .false.;
	periodic = .true.; onlybulk = .false.;
	leads = .false.;
	Ebr = 0.0d0; Ebl = 0.0d0; ! R/L contact barriers for electrons 
	!Eblmin=0.0d0; Eblmax=0.0d0; nEbl=1;
	!Ebrmin=0.0d0; Ebrmax=0.0d0; nEbr=1;
	Er = 0.0d0;
	Ermin=0.0; Ermax=0.0; dEr=0.0; ner = 1;
	w0 = 2.0d0; ! exciton energy: E_LUMO = w0+Exb;
	! we are ignoring electron-electron repulsion U : E_Double_occupied = E_LUMO = w0+Exb;
	!	We can absorb U in Exb term, since the only states that matter are excited and D. (?)
	Exb = 0.0; ! exciton binding energy
	dwmin=0.0; dwmax=0.0; ddw=0.0; ndw = 1;
	kappa = 0.1 !0.005d0 ! ref to th=1.0d0; scale times proportionally tp th.
	gamma = 0.1 !0.005d0
	beta = 40.0d0;
	crosshops = .true.
	AlwaysLP = .false.
	PermSym = .false.
	nokappa = .false.
	nogamma = .false.	
	nolosses = .false.
	mincarriers = .true.
	givenNexcit = .false.;
	ztout= .true.;
	ntrapout= .true.;
	ratesout= .true.;
	!onlydoped = .true.
	simplepf = .true.
	givenEr = .false.
	VRH = .false.
	coulomb = .false.
	impurity = .false.
	vqout = .false.
	!--------------
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

	case('OneDChains')
		read(50,*,err=20) OneDchains
		EqualDistr = .true.
		
	case('debug')
		read(50,*,err=20) debug

	case('vqout')
		read(50,*,err=20) vqout

	case('lossgain')
		read(50,*,err=20) lossgain

	case('ztout')
		read(50,*,err=20) ztout

	case('ntrapout')
		read(50,*,err=20)ntrapout

	case('ratesout')
		read(50,*,err=20) ratesout


	case('Nsites')
		read(50,*,err=20) nsites
		nsites = nsites + 1; ! maz: dhops; add a site so that
		! 'Nsites' will be used for total active sites in the model.
		if(1==0 .and. mod(nsites,3) .ne. 0) then
			write(*,*)"Warning(readinp): Nsites increased to multiple of 3"
			write(*,*)	"Warning(readinp): given Nsites=",Nsites
			nsites = 3*(nsites/3) + 3; ! integer arithmatic
			write(*,*)	"Warning(readinp): used Nsites=",Nsites
		endif
	!case('Doping','doping')
	!	read(50,*,err=20) nel
	!	givenNel = .true.

	case('DopingRange','Dopingrange','dopingrange')
		read(50,*,err=20) nelmin, nelmax,dnelec
		givenNelrange = .true.;


	case('SimplePenaltyFunction','simplepf')
		read(50,*,err=20) simplepf ! if false, use Ohmic bath spectral density
	case('BathSpectralDensityCutoff','wcut')
		read(50,*,err=20) wcut

	case('ChargingEnergy')
		read(50,*,err=20) Eq
		! Eq > 0 only
		if(Eq < 0.0d0) then
			write(*,*)"Warning(readinput): ChargingEnergy should be >= 0"
			write(*,*)"==> Finite e affinity cannot keep charging..."
			write(*,*)"==> Or change the form of dEQs to include it."	
		endif

	case('NExcitationsRange','NexcitationsRange','nexcitationsrange')
		read(50,*,err=20) mexmin,mexmax,dmex
		givenNexcit = .true.;

	case('Detuning')
		read(50,*,err=20) dwmin,dwmax,ndw
		givenDw = .true.;
		
	case('NExcitations','Nexcitations')
		read(50,*,err=20) nx

	case('Niter')
		read(50,*,err=20) niter

	!case('Ntraj','ntraj')
	!	read(50,*,err=20) ntraj

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
		read(50,*,err=20) gmin, gmax, ng
		giveng = .true.

	case('Temperature') ! in Kelvin
		read(50,*,err=20) Tmin, Tmax, nTemp
		givenT = .true.

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

	case('ExcitonEnergy','w0')
		read(50,*,err=20) w0

	case('ExcitonBindingEnergy','Exb')
		read(50,*,err=20) Exb
		giveExb = .true.

	case('DielectricConstant','Epsilonr')
		read(50,*,err=20) epsr

	case('Kappa','kappa')
		read(50,*,err=20) kappa

	case('Gamma','gamma')
		read(50,*,err=20) gamma

	case('NoLosses','Nolosses','nolosses')
		read(50,*,err=20) nolosses

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



		!'DevicesTypes'

		dev%ntyps=1
		!dev%ntyps: total number of types (eonly, honly, etc...)
		allocate(dev%typ(dev%ntyps))		! type of the device
		allocate(dev%ndev(dev%ntyps)) ! tot devices of a given type
		allocate(dev%X(dev%ntyps))
		do i=1,dev%ntyps
			dev%typ(i)=1; dev%ndev(i)=1;
			!------------------------------------
			!select case(trim(devtyp))
			!	case('e-only')
			!		dev%typ(i) = 1
			!	case('h-only')
			!		dev%typ(i) = 2
			!	case('e-h')
			!		dev%typ(i) = 3
			!end select
			!------------------------------------			
			allocate(dev%X(i)%Eb(dev%ndev(i),2))
			do j=1,dev%ndev(i)
				dev%X(i)%Eb(j,:) = 0.0
			enddo
		enddo
		dev%ntot = sum(dev%ndev(:)) ! total number of devices
		allocate(dev%idv(dev%ntot,2)) ! indexes table, to avoid two loops in main program.
		idv=0
		do i=1,dev%ntyps
			do j=1,dev%ndev(i)
				idv = idv + 1;
				dev%idv(idv,1) = i;
				dev%idv(idv,2) = j;
			enddo
		enddo



	lossgain=0
	if(lossgain==1)then
		nolosses = .false.
		crosshops=.true.
	elseif(lossgain==0)then
		nolosses = .true.
		crosshops=.false.
	endif

		if (.not. OneDchains) then
			if(.not. impurity) then
				nproc = 34
			elseif(impurity) then
				nproc=42
			endif
		else ! OneDChain
			if(.not. impurity) then
				nproc = 26
			elseif(impurity) then
				write(*,*)"Error(reainp): OneDchains, no impurity only..."
				write(*,*)"Error(reainp): Not implemented...."
				stop
			endif
		endif

	if(coulomb .and. nsites < 9) then
		write(*,*) "Error(readinput):"
		write(*,*) "Coulomb is not implemented for nsites < 9"
		stop
	elseif(coulomb .and. periodic) then
		write(*,*) "Error(readinput):"
		write(*,*) "Coulomb is not implemented for periodic "
		write(*,*) "TODO: Ewald method? "
		stop
	endif


	! Coulomb law constant K for the medium for rij in nm, q's in |e|.
	Kq = 1.4224d0/epsr; 	! default epsr = 3.0d0; ===> Kq = 0.4741; 
											! 0.47 eV.nm energy for two elementary charges
											! ~0.1 eV for two nnz at 5nm 
	if(coulomb .and. .not. giveExb) then ! set Exb assuming a radius of 1nm for the exciton
		! todo?: ask for exciton radius to compute Exb
		Exb = Kq; ! Kq(eV.nm/e)*1e/1nm ==> default = 0.4741 eV
	elseif(coulomb .and. Exb < Kq/a0) then !  .and. a0 > 1.0d0 ??
		write(*,*)"Serious Warning(readinp): Exb < Kq/a0 "
		write(*,*)"Electron and holes at nnz might have lower energy "
		write(*,*)"   compared to an exciton with this binding energy"
		write(*,*)"For exciton radius ~ 1nm, Exb should be ~ Kq "
	endif

! set energy of impurity level based on its type; after setting/reading Exb above
	if(impurity) then
		! impocc: neutral state configuraton/occupation of impurity site/level
		if(imptype==1) then
			Eimp = w0+Exb-Eimp0; ! Below LUMO
			impocc = 0; ! will capture an electron, hard to release it.
		elseif(imptype==2) then
			Eimp = Eimp0; ! Above HOMO
			impocc = 1; ! will capture a hole, will be hard to release its hole
		elseif(imptype==3) then
			Eimp = w0+Exb+Eimp0; ! Above LUMO
			impocc = 1; ! will release this and become +ve charged
		elseif(imptype==4) then
			Eimp = -Eimp; ! Below HOMO
			impocc = 0; ! will capture an electron and become -ve charged
		else
			write(*,*)"Error(readinp) : Impurity type 1-4 only"
			write(*,*)"1 for e & 2 for h trap, 3 for n & 4 for p dopant"
			stop
		endif	
		write(*,*)"impocc, Eimp0,Eimp=",impocc, Eimp0,Eimp
	endif
	
	
	call VarRange(dwmin,dwmax,ndw,ddw,givendw,1.0d-5,0.0d0)

	! last arg= default value if a single value is to be used as default
	call VarRange(Tmin,Tmax,nTemp,dtemp,givent,1.0d-5,1/(kb*beta))
	call VarRange(gmin,gmax,ng,dg,giveng,1.0d-5,g)
	!call VarRange(Eblmin,Eblmax,nEbl,dEbl,givenEblr,0.0d0,Ebl)
	!call VarRange(Ebrmin,Ebrmax,nEbr,dEbr,givenEblr,0.0d0,Ebr)

	! contact/leads present?
	leads = (.not. periodic) .and. (.not. onlybulk);
	if(leads) onlydoped = .false.
	!if (sum(abs(JhR+JlR+JhL+JlL)) < 1.0d-5) then
	!	leads = .false.
	!	onlybulk = .true.
	!endif

	!write(*,*)"readinput: periodic, onlybulk, leads = ",
  !   .       periodic, onlybulk, leads


	if(.not. givenNexcit) then
		mexmin=nx;
		mexmax=nx;
		dmex=1;
	endif


	!--------------------------
	! keep the order 
	!--------------------------
	!if(abs(g)<1.0d-2) nog = .true.
	if (nog) then
		AlwaysLP = .false. ! stay in the state (= ksub set) after a transition
		sameg = .true. ! g=0 for all molecules
	endif
	! nog ===> AlwaysLP=F ===> PermSym = F
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
	if(periodic) then
		nelmin=nel; nelmax=nel; dnelec=1;
	else
		nelmin=0; nelmax=0; dnelec=1;
	endif
	endif

	if(nelmin <= -nsites .or. nelmax >= nsites) then
		nelmin = max(nelmin, -(nsites-1) );
		nelmax = min(nelmax, nsites-1);
		!dnelec=1;
		write(*,*)"Warning(readinput): DopingRange invalid!"
		write(*,*)" changed to sensible min, max values."
	endif

	!write(*,*)"nelmin, nelmax, dnelec = ",nelmin, nelmax, dnelec

	!---------------------------------------------
	! write various parameters
	i = 1+int((mexmax-mexmin)/dmex);
	j	=	1+int((nelmax-nelmin)/dnelec)

	if (onlydoped) j=1+int((nelmax-nelmin-1)/dnelec);
	if(node==0) then
		open(200,file='parameters.out',action='write')
		write(200,*) i,ndw,j,ntraj,niter
		write(200,*) mexmin,mexmax,dmex
		write(200,*) dwmin,dwmax, ddw
		write(200,*) nelmin,nelmax,dnelec
		write(200,*) ntraj, niter
		close(200)
	endif
	!---------------------------------------------
	return
	end subroutine input
!************************************************
	subroutine giveinput(str)
	implicit none
	character(len=*) :: str
	write(*,'("Input variable ",a," not given, stopping!")') str
	stop ! stop the code
	return
	end 	subroutine giveinput
!************************************************
! VarRange sets the value of increment, 
!						and number of points/divisions
	subroutine VarRange(dwmin,dwmax,ndw,ddw,given,tol,def)
	implicit none
	double precision, intent(inout):: dwmin,dwmax
	double precision, intent(in):: tol,def
	logical, intent(in):: given
	integer, intent(inout):: ndw 
	double precision, intent(out)::ddw
	if(.not. given) then
		dwmin = def; dwmax=def;
		ddw=0.0; ndw=1;
	elseif(ndw==1 .or. abs(dwmax-dwmin) < tol)then
		ddw=0.0; ndw=1;
	else
		ddw = (dwmax-dwmin)/(ndw-1);
	endif

	return
	end 	subroutine VarRange
!************************************************

	end module readinput


