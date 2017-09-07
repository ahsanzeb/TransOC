
	module diag
	implicit none

	public :: iterdiag, ddiag
	private:: matvec

	contains
	!-----------------------------------
	subroutine iterdiag(it,n,nev,ncv)
	use modmain, only: diagmaxitr
	implicit none
	integer, intent(in) :: it ! which hamiltonian
	integer, intent(in) :: n ! n x n dimension of Hg
	integer, intent(in) :: nev ! number of eigenpairs to compute
	integer, intent(in) :: ncv ! leading dimensions for all arrays
	integer :: maxn, maxnev, maxncv,ldv
	! arrays
	double precision,allocatable, dimension(:,:) :: v
	double precision,allocatable, dimension(:) ::  workl
	double precision,allocatable, dimension(:) :: workd
	double precision,allocatable, dimension(:,:) :: d
	double precision,allocatable, dimension(:) :: resid

	logical, allocatable, dimension(:) :: selectt
	integer :: iparam(11), ipntr(11)

	! scalars
	character :: bmat*1, which*2
	integer :: ido, lworkl, info, ierr, j
	integer :: nconv, maxitr, mode, ishfts
	!logical :: rvec
	double precision :: tol, sigma, zero = 0.0d+0

c	dsaupd documentation:
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Lanczos vectors are generated, the algorithm generates 
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the 
c          matrix-vector product OP*x. (See remark 4 there).

	!n = nx*nx
	maxn = n;
	!ncv = 2*nev; ! ncv <= n,  

	maxnev = nev; 
	maxncv= ncv;
	ldv = maxn;

	allocate(v(ldv,maxncv))
	allocate(workl(maxncv*(maxncv+8)))
	allocate(workd(3*maxn))
	allocate(d(maxncv,2))
	allocate(resid(maxn))
	allocate(selectt(maxncv))

	! irrelevant messages 	
	if ( n .gt. maxn ) then
		print *, ' ERROR with _SDRV1: N is greater than MAXN '
		go to 9000
	else if ( nev .gt. maxnev ) then
		print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
		go to 9000
	else if ( ncv .gt. maxncv ) then
		print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
		go to 9000
	end if

	! requirements:
	!	NEV + 1 <= NCV <= MAXNCV	 

	which = 'SM'; ! 'SA'
					! SA ???? How to find lowest eigenvalues
					! is Hg positive definite??? 
					! unless the diagonal detuning makes it negative	
					! make Hg positive definite by subtracting its norm etc??


	lworkl = ncv*(ncv+8) ! at least ncv*(ncv+8)
	tol = zero 
	info = 0
	ido = 0

	iparam(1) = 1 ! ishfts 
	iparam(3) = diagmaxitr ! max iterations 
	iparam(7) = 1 ! mode, standard eigenvalue problem

	!write(*,*)" ---- started ------- "

		!-------------------------------------
		!Beginning of reverse communication
		!call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
		!&	 		 	 	  ncv, v, ldv, iparam, ipntr, workd, workl,
		!&	 		 	 	  lworkl, info )
10		call dsaupd ( ido, 'I', n, which, nev, tol, resid, 
     .   ncv, v, ldv, iparam, ipntr, workd, workl,
     .   lworkl, info )

	!write(*,*)" ---- ido = ",ido


		if (ido .eq. -1 .or. ido .eq. 1) then
		
			!write(*,*)"shape(workd(ipntr(1))),shape(workd(ipntr(2)))"			
			!write(*,*)shape(workd(ipntr(1):ipntr(1)+n-1))
			!write(*,*)shape(workd(ipntr(2):ipntr(2)+n-1))		
			call matvec(it, n, workd(ipntr(1):ipntr(1)+n-1),
     .                    workd(ipntr(2):ipntr(2)+n-1))
			!write(*,*)" ---- matvec iterations .... "

			!write(*,*)"diag: testing -  stop"
			!stop
			
			go to 10
		end if 
		!End of Reverse communication
		!-------------------------------------

	write(*,*)" =======>>>>> info = ",info


	if ( info .lt. 0 ) then ! see documentation of DSAUPD
		write(*,*)'Error(diag): Error with _saupd, info = ', info
		write(*,*)'Error(diag): Check documentation in _saupd '
		stop
	else
		!call dseupd (.true., 'All', selectt, d, v, ldv, sigma, 
		! bmat, n, which, nev, tol, resid, ncv, v, ldv, 
		! iparam, ipntr, workd, workl, lworkl, ierr )
		call dseupd (.true., 'All', selectt, d, v, ldv, sigma, 
     .   'I', n, which, nev, tol, resid, ncv, v, ldv, 
     .   iparam, ipntr, workd, workl, lworkl, ierr )

	endif

	if ( ierr .lt. 0 ) then ! see documentation of DSEUPD
		write(*,*)'Error(diag): Error with _seupd, info = ', ierr
		write(*,*)'Error(diag): Check documentation in _seupd '
		stop
	endif


	write(*,*)" diag : -------------- done--------------------"
	
	
	
	! set eig
	! eig(it)%ntot = 
	! eig(it)%n1 = 
	! eig(it)%n2 = 
	! eig(it)%evec = 
	! eig(it)%eval = 

9000		continue

	deallocate(v,workl,workd,d,resid,selectt)
	
	return
	end subroutine iterdiag

	!-----------------------------------
	subroutine matvec(it,n,x,y)
	! multiplies Hamiltonian stored in slot it with input vector x
	! gives y as output: y = H.x
	! Hamiltonians are stored in CSR format
	use modmain, only: Hg ! Hg dim = n x n
	implicit none
	integer, intent(in) :: it,n
	double precision, dimension(n), intent(in):: x
	double precision, dimension(n), intent(out):: y
	integer :: i,j,k

	write(*,*)"matvec: Hg(it)%xst = ",Hg(it)%xst
	write(*,*)"matvec: Hg(it)%ntot = ",Hg(it)%ntot
	write(*,*)"matvec: Hg(it)%nnz = ",Hg(it)%nnz
	write(*,*)" maxval col: "	, maxval(Hg(it)%col)
	write(*,*)" maxval rowptr: "	, maxval(Hg(it)%rowpntr)
	write(*,*)Hg(it)%rowpntr
	! Hg only upper triangular
	! below the diagonal elements are double counted
	! to fix this we store half of the diagonal elements in Hg(it)%dat
	! instead of using an if statement here in the loop.
	y = 0.0d0
 	do i=1,Hg(it)%srptr-1
		do j=Hg(it)%rowpntr(i), Hg(it)%rowpntr(i+1)-1
			k = Hg(it)%col(j); ! icol
			y(i) = y(i) + x(k)*Hg(it)%dat(j)
			y(k) = y(k) + x(i)*Hg(it)%dat(j) 
		enddo
 	enddo
 	
	return
	end subroutine matvec
!-----------------------------------

	subroutine ddiag(H,W,ntot)
	use modmain, only: clock !
	implicit none
	integer, intent(in) :: ntot
	double precision, dimension(ntot,ntot), intent(inout):: H	
	double precision, dimension(ntot), intent(inout):: W
	! local
	integer :: lwork != 3*ntot ! LWORK >= max(1,3*N-1)
	double precision, dimension(3*ntot):: WORK
	integer :: info
	double precision ::starttime, endtime,starttime1, endtime1
	EXTERNAL	 DSYEV ! LAPACK/BLAS

	write(*,'(/,a)') " diagonal: diagonalising H ... "
	starttime1 = clock()

	lwork = 3*ntot ! LWORK >= max(1,3*N-1)
	
	CALL DSYEV('Vectors','Upper',ntot,H,ntot,W,WORK,LWORK,INFO)
	! Check for convergence.
	IF( INFO.GT.0 ) THEN
		WRITE(*,'(/,a)')
     .   'diagonal: The algorithm failed to compute eigenvalues.'
		STOP
	END IF

	! ******** set global variable eig *********


	endtime1 = clock()
	write(*,'(/,a,E8.2)')" diagonal: total time taken	= ",
     . 			endtime1-starttime1

	return
	END subroutine ddiag
!-----------------------------------------------

	
	end 	module diag
