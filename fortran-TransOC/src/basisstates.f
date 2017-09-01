
	module basisstates
	use modmain
	
	implicit none

	contains

	subroutine mkbasis(na,nx)
	! this routine calculates the index pointers
	! and subsets for all 5 basis defined in modmain

	! if a range of N above and below the N-2,N+2
	! are used to save processing on cost of memory, then 
	! the routine might build up extra sectors that are rarely used
	! otherwise at every iteration of dyanamics, some of N will go out of the range
	!	and basis will be calculated for them, so overall this will avoid building up
	!	extra ksub sectors.
	
	!use modmain, only mapb
	implicit none
	integer(kind=1), intent (in) :: na,nx
	integer(kind=4) :: ntot
	integer(kind=1) :: nxmax,r,m1,i,j,n,l
	integer(kind=1), dimension(5):: nalist5
	integer  :: i1,i2, maxk,ibl,ib
	logical :: calc
	type(BSectors), allocatable :: sec(:)

	! list of N for 5 cases
	nalist5 = na + dns
	! max nx possible in 13 cases => nx+1 or nx+2
	nxmax = nx+1;
	if(crosshops) nxmax = nx+2;
	write(*,*) "nxmax, nx = ",nxmax, nx
	! indexes pointers for basis set sectors with diff no of up spins

	do i=1,5
		n = nalist5(i);
		m1 = min(nalist5(i),nxmax); ! max up spin possible
		ibl = mapb%map(i); ! localtion of this ib=i
		! calculate full basis or update for some extra k-subsets?
		calc = .false.	
		do ib=1,5
			if (mapb%cal(ib) == i) then
				calc = .true.
				exit
			endif
		end do
		!-------------------------------------------
		if (calc) then
			! update maxk
			basis(ibl)%maxk = m1;
			! indexed pointers
			if(allocated(basis(ibl)%pntr)) deallocate(basis(ibl)%pntr)
			allocate(basis(ibl)%pntr(m1+2))
			call pointerslist(nalist5(i),m1,basis(ibl)%pntr)
			! subsets
			if(allocated(basis(ibl)%sec)) deallocate(basis(ibl)%sec)
			allocate(basis(ibl)%sec(m1)) ! m1 sectors
			!	ignore the 0-up (all down) sector (a single basis state)
			!	start j=1:m1 for sectors with j up spins (2nd onwards)
			do j=1,m1,1	
				ntot = basis(ibl)%pntr(j+2) - basis(ibl)%pntr(j+1); ! for config of this type
				if(allocated(basis(ibl)%sec(j)%sets))
     .								deallocate( basis(ibl)%sec(j)%sets )
				allocate(basis(ibl)%sec(j)%sets(ntot,j))
				call mksets(n,j,ntot,basis(ibl)%sec(j)%sets)
			end do
		!-------------------------------------------
		else ! update
			maxk = basis(ibl)%maxk ! subsets up to maxk available
			!	update the subsets etc only if m1 > maxk
			if (m1 > maxk) then
				! update maxk
				basis(ibl)%maxk = m1; 
				! pntr not costly, so compute full array
				deallocate(basis(ibl)%pntr)
				allocate(basis(ibl)%pntr(m1+2))
				call pointerslist(nalist5(i),m1,basis(ibl)%pntr)
				! subsets
				!if (allocated(sec)) deallocate(sec)
				allocate(sec(maxk))
				sec = basis(ibl)%sec ! save current value
				deallocate(basis(ibl)%sec)
				allocate(basis(ibl)%sec(m1)) ! m1 sectors
				basis(ibl)%sec(1:maxk) = sec ! set to prev value
				deallocate(sec)
				! compute extra sectors
				do j=maxk+1,m1,1
					ntot = basis(ibl)%pntr(j+2) - basis(ibl)%pntr(j+1); ! for config of this type
					allocate(basis(ibl)%sec(j)%sets(ntot,j))
					call mksets(n,j,ntot,basis(ibl)%sec(j)%sets)
				end do
			end if

		endif ! calc or update
		!-------------------------------------------
	end do ! i [ib]

	end subroutine mkbasis


!------------------------------------------

	subroutine mksets(n,k,ntot,comb)
	! genrates subsets
	implicit none
	integer(kind=1), intent (in) :: n,k
	integer(kind=4), intent (in) :: ntot
	integer(kind=1), dimension(ntot,k),intent(out):: comb
	integer :: ind,i
	integer(kind=1) :: m, m2
	logical mtc
	integer(kind=1), dimension(k):: a
	
	comb(:,:)=0;
	!----------- k=n ------------
	if (k==n) then
		comb(1,:) = (/ (ind,ind=1,n,1) /);
		return
	endif
	!----------- k<n ------------
	ind = 1;	
	mtc=.true.
	do while (mtc)
		if (ind==1) 	mtc=.false.
		call ksub_next(n,k,a,mtc,m, m2) ! ( n, k, a, more, m, m2 )
		!write(*,*) (a(i),i=1,k)
		comb(ind,:) = a(:)
		ind = ind + 1;
	end do
	return
	end subroutine
	!------------------------------

	subroutine ksub_next ( n, k, a, more, m, m2 )

!*****************************************************************************80
!
!! KSUB_NEXT generates the subsets of size K from a set of size N.
!
!	Licensing:
!
!		This code is distributed under the GNU LGPL license.
!
!	Modified:
!
!		26 May 2015
!
!	Author:
!
!		Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!		FORTRAN90 version by John Burkardt.
!
!	Reference:
!
!		Albert Nijenhuis, Herbert Wilf,
!		Combinatorial Algorithms for Computers and Calculators,
!		Second Edition,
!		Academic Press, 1978,
!		ISBN: 0-12-519260-6,
!		LC: QA164.N54.
!
!	Parameters:
!
!		Input, integer ( kind = 4 ) N, the size of the set from which subsets
!		are drawn.
!
!		Input, integer ( kind = 4 ) K, the desired size of the subsets.	K must
!		be between 0 and N.
!
!		Input/output, integer ( kind = 4 ) A(K).	A(I) is the I-th element of the
!		subset.	Thus A(I) will be an integer between 1 and N.
!		Note that the routine will return the values in A
!		in sorted order: 1 <= A(1) < A(2) < ... < A(K) <= N
!
!		Input/output, logical MORE.	Set MORE = FALSE before first call
!		for a new sequence of subsets.	It then is set and remains
!		TRUE as long as the subset computed on this call is not the
!		final one.	When the final subset is computed, MORE is set to
!		FALSE as a signal that the computation is done.
!
!		Input/output, integer ( kind = 4 ) M, M2, two variables used by this
!		procedure for bookkeeping.	The user must declare these variables,
!		and the output values from one call must be used as the input values
!		on the next.	The user should not change these values.
!
	implicit none

	integer ( kind = 1 ) k

	integer ( kind = 1 ) a(k)
	integer ( kind = 1 ) j
	integer ( kind = 1 ) m
	integer ( kind = 1 ) m2
	logical more
	integer ( kind = 1 ) n

	if ( k < 0 .or. n < k ) then
		write ( *, '(a)' ) ''
		write ( *, '(a)' ) 'KSUB_NEXT - Fatal error!'
		write ( *, '(a,i8)' ) 'N = ', n
		write ( *, '(a,i8)' ) 'K = ', k
		write ( *, '(a)' ) 'but 0 <= K <= N is required!'
		stop 1
	end if

	if ( .not. more ) then
		m2 = 0
		m = k
	else
		if ( m2 < n - m ) then
			m = 0
		end if
		m = m + 1
		m2 = a(k+1-m)
	end if

	do j = 1, m
		a(k+j-m) = m2 + j
	end do

	more = ( a(1) /= (n-k+1) )

	return
	end

!------------------------------------------
	subroutine pointerslist(n,k,pntr)
	! pointer indexed to the subsets with different number of up spins
	! pntr = Insert[Table[Sum[Bimonial[n, j], {j, 0, i}], {i, 0, k}], 0, 1];
	integer(kind=1), intent(in) :: n,k
	integer(kind=4), intent(out) :: pntr(k+2)
	integer(kind=1):: i,r
	
	pntr(1) = 0; pntr(2) = 1;
	do i=3,k+2,1
		r = i-2;
		pntr(i)=pntr(i-1)+ nCr(n,r)
		!write(*,*) "k,i, pntr(i) = ",k,i, pntr(i)
	end do

	return
	end subroutine pointerslist
!------------------------------------------
	integer function nCr(n,r)
	! ^nC_r
	integer(kind=1) n,r,i,large,small
	integer(kind=8) Numer,Denom 

	 !define =0 for n<r ??? only needed in LexicoIndex
	 !or just dont call it when n<r ????
	if (n<r) then
		nCr=0
		return
	endif

	if (n==r .or. r==0) then
		nCr=1
		return
	elseif (r.lt.(n-r)) then
		small=r
		large=n-r
	else
		large=r
		small=n-r
	end if
          
	Numer=n
	do i=n-1,(large+1),-1
		Numer = Numer*i
	end do
          
	Denom=small
	do i=small-1,2,-1
		Denom = Denom* i
	end do
          
	nCr=Numer/Denom

	return   
	end function

!------------------------------------------
	integer(kind=4) function LexicoIndex(list,n,k)
	implicit none
	integer(kind=1), intent(in) :: n,k
	integer(kind=1), dimension(k), intent(in)::list
	integer(kind=1) :: p,np,kp
	!x1 = c[N, k] - Sum[c[N - listf[[p + 1]], k - p], {p, 0, k - 1}];
	LexicoIndex = nCr(n,k);
	do p=0,k-1,1
		np = n - list(p+1);
		kp = k - p;
		LexicoIndex = LexicoIndex - nCr(np,kp) ! only if (np .ge. kp )	
	end do
	return
	end function
!------------------------------------------


	!subroutine Shift(set,k,l)
	!implicit none
	!integer, intent(in):: k,l
	!integer, dimension(k), intent(inout) :: set
	! local
	!integer :: i
	!set(l:k) = 	set(l:k) + 1;
	!end subroutine Shift


!------------------------------------------
	end module
