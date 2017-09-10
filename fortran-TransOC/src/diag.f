
	module diag
	implicit none

	public :: diagonalise
	private	:: iterdiag, ddiag
	private	:: matvec,matvec1,DegenSectors

	contains
!------------------------------------------
!	call to diagonalise() will diagonalise
! all newly calculated Hamiltonians
	subroutine diagonalise()
	use modmain, only: mapt,Hg,eig
	implicit none
	! local
	integer:: i,ntot,ib, nt,it,itype,nev,ncv

		! use better storage format (and matvec routines for iter diag)
		do ib=1,5
			nt = mapt%ntb(ib) ! itypes (for same n, not imp here)
												! 		for which H is calculated & eig is needed
			do i=1,nt
				itype = mapt%grouptb(ib,i); ! which of 13 types?
				it = mapt%map(itype);! map for the location of itype
				ntot = Hg(it)%ntot;
				!write(*,*)"diag: it = ",it
				! calculate eigenpairs
				if(Hg(it)%dense) then
					! find all eigenpairs
					eig(it)%ntot = ntot ! dim of final hilbert space	
					eig(it)%n1 = ntot ! dim of hilbert space
					eig(it)%n2 = ntot; ! number of vectors
					nev = ntot;
					!write(*,*)"diag: direct: it, Hg(it)%dense",it,Hg(it)%dense

					! upper triangular H stored in evec, & eval allocated
					call ddiag(eig(it)%evec,eig(it)%eval,ntot)
				else
					!write(*,*)"diag: iter: it, Hg(it)%dense",it, Hg(it)%dense

					! iterdiag will find the storage format of H
					!	and use appropriate matvec routines
					! it will set the global variables 
					!	eig(it)%evec and eig(it)%eval
					nev = Hg(it)%nev;
					ncv = Hg(it)%ncv;					
					call iterdiag(it, ntot, nev, ncv)
					write(*,*) "diag: iter done.... "
					write(*,*) "=====> it, ntot, nev,ncv =",it,ntot, nev,ncv
				endif

				!write(*,*)"============ it = ",it
				! calculate degenerate sectors
				if(allocated(eig(it)%esec))deallocate(eig(it)%esec)
				if(allocated(eig(it)%ind))deallocate(eig(it)%ind)
				allocate(eig(it)%esec(nev)) ! max nev sec if no degeneracy
				allocate(eig(it)%ind(nev+1)) ! indexed for sectors

				eig(it)%ind = -10
				!write(*,*)"diag: nev",nev
				call DegenSectors(eig(it)%eval,nev,
     .   eig(it)%nsec,eig(it)%esec,eig(it)%ind) ! make degenerate sectors
				
			enddo
		end do

	return
	end	subroutine diagonalise
!------------------------------------------
	subroutine iterdiag(it,n,nev,ncv)
	use modmain, only: diagmaxitr, sameg, detuning,eig
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
	logical :: better
	double precision :: tol, sigma, zero = 0.0d+0


	! use more efficient storage and matvec multiplication?
	better = (.not. detuning) .and. sameg;



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

	! requirements:
	!	NEV + 1 <= NCV	 <= n 

	which = 'SA'; ! lowest eigenvalues?

	lworkl = ncv*(ncv+8) ! at least ncv*(ncv+8)
	tol = zero 
	info = 0
	ido = 0

	iparam(1) = 1 ! ishfts 
	iparam(3) = diagmaxitr ! max iterations 
	iparam(7) = 1 ! mode, standard eigenvalue problem

	!write(*,*)" ---- started ---- n,nev,ncv ",n,nev,ncv

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
			!write(*,*)" -------- 1 ------ "
			!write(*,*)"shape(workd(ipntr(1))),shape(workd(ipntr(2)))"			
			!write(*,*)shape(workd(ipntr(1):ipntr(1)+n-1))
			!write(*,*)shape(workd(ipntr(2):ipntr(2)+n-1))	
			if(better) then
				call matvec1(it, n, workd(ipntr(1):ipntr(1)+n-1),
     .                    workd(ipntr(2):ipntr(2)+n-1))
			else
				!write(*,*)" -------- 2 ------ "

				call matvec(it, n, workd(ipntr(1):ipntr(1)+n-1),
     .                    workd(ipntr(2):ipntr(2)+n-1))
     	endif
     
			!write(*,*)" ---- matvec iterations .... "

			!write(*,*)"diag: testing -  stop"
			!stop
			
			go to 10
		end if 
		!End of Reverse communication
		!-------------------------------------

	!write(*,*)" =======>>>>> info = ",info


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
c        %-------------------  dseupd  -----------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NEV columns of the two dimensional |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
	endif

	if ( ierr .lt. 0 ) then ! see documentation of DSEUPD
		write(*,*)'Error(diag): Error with _seupd, info = ', ierr
		write(*,*)'Error(diag): Check documentation in _seupd '
		stop
	endif


	!write(*,*)" diag : -------------- done--------------------"
	
	! set eig
	eig(it)%ntot = n
	eig(it)%n1 = n
	eig(it)%n2 = nev
	if(allocated(eig(it)%evec))deallocate(eig(it)%evec)
	if(allocated(eig(it)%eval))deallocate(eig(it)%eval)
	allocate(eig(it)%evec(n,nev))
	allocate(eig(it)%eval(nev))
	eig(it)%evec = v(:,1:nev)
	eig(it)%eval = d(1:nev,1) ! dim of d = (maxncv,2), so first nev values?

	!write(*,*) "diag: n1,n2 = ",n,nev

	deallocate(v,workl,workd,d,resid,selectt)

9000		continue

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

	!write(*,*)"matvec: Hg(it)%xst = ",Hg(it)%xst
	!write(*,*)"matvec: Hg(it)%ntot = ",Hg(it)%ntot
	!write(*,*)"matvec: Hg(it)%nnz = ",Hg(it)%nnz
	!write(*,*)" maxval col: "	, maxval(Hg(it)%col)
	!write(*,*)" maxval rowptr: "	, maxval(Hg(it)%rowpntr)
	!write(*,*)"-----dat------"
	!write(*,*)Hg(it)%dat
	!write(*,*)"-----col------"
	!write(*,*)Hg(it)%col
	!write(*,*)"-----rowptr------"
	!write(*,*)Hg(it)%rowpntr
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

	!-----------------------------------
	subroutine matvec1(it,n,x,y)
	! multiplies Hamiltonian stored in slot it with input vector x
	! gives y as output: y = H.x
	! Hamiltonians are stored in CSR format
	! used the fact that matrix elements in
	! a ksub sector has the same value for all non-zero elements
	! when {g_i}=g; all molecules light couplings same
	use modmain, only: Hg ! Hg dim = n x n
	implicit none
	integer, intent(in) :: it,n
	double precision, dimension(n), intent(in):: x
	double precision, dimension(n), intent(out):: y
	integer :: i,j,k,s
	double precision:: yi, val
	
	!write(*,*)"matvec1: Hg(it)%xst = ",Hg(it)%xst
	!write(*,*)"matvec1: Hg(it)%ntot = ",Hg(it)%ntot
	!write(*,*)"matvec1: Hg(it)%nnz = ",Hg(it)%nnz
	!write(*,*)"matvec1: spntr = ",Hg(it)%spntr
	!write(*,*)" maxval col: "	, maxval(Hg(it)%col)
	!write(*,*)" maxval rowptr: "	, maxval(Hg(it)%rowpntr)
	!write(*,*)Hg(it)%rowpntr

	! Hg only upper triangular
	! below the diagonal elements are double counted
	! to fix this we store half of the diagonal elements in Hg(it)%dat
	! instead of using an if statement here in the loop.
	y = 0.0d0

	do s=1,Hg(it)%m1-1
		val = Hg(it)%dat(s);
 		do i= Hg(it)%spntr(s), Hg(it)%spntr(s+1)-1 !1,Hg(it)%srptr-1
			yi = 0.0d0;
			do j=Hg(it)%rowpntr(i), Hg(it)%rowpntr(i+1)-1
				k = Hg(it)%col(j); ! icol
				yi = yi + x(k); ! multiply val later
				y(k) = y(k) + x(i)*val
			enddo
			y(i) = y(i) + yi*val;
 		enddo
 	enddo
	return
	end subroutine matvec1
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

	!write(*,'(/,a)') " diagonal: diagonalising H ... "
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
	!write(*,'(/,a,E8.2)')" diagonal: total time taken	= ",
  !   . 			endtime1-starttime1

	return
	END subroutine ddiag
!-----------------------------------------------
	subroutine DegenSectors(es,ne,nsec,esec,ind)
	! importnat:only if es are in ascending order
	integer, intent(in) :: ne
	double precision,dimension(ne),intent(in):: es
	double precision,dimension(ne),intent(out):: esec ! 1:nsec contains esec
	integer, intent(out) :: nsec
	integer,dimension(ne+1),intent(out):: ind !start index of sectors; 1:nsec
	! local
	double precision:: tol = 1.0d-3 ! tolerance
	integer:: i,j

	!write(*,*)"diag: e = ",es

	esec(1) = es(1); j = 1; ! start with the lowest energy
	ind(1) = 1;
	do i=1,ne
		if (es(i) > esec(j) + tol ) then
			j = j + 1; 
			esec(j) = es(i);
			ind(j) = i;
		endif		
	end do

	nsec = j; ! total number of sectors
	!if(j<ne)
	ind(j+1) = ne+1; ! set last element of ind
	return
	end subroutine DegenSectors
!---------------------------------------


	
	end 	module diag
