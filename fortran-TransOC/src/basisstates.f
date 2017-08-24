	!program combinations
	module basisstates
	use modmain
	
	implicit none

	contains

	subroutine mkbasis(na,nx)
	! this routine calculates the index pointers
	! and subsets for all 5 basis defined in modmain
	implicit none
	integer, intent (in) :: na,nx
	integer :: r,m1,ntot

	!nalist = na + dna
	!nxlist = nx + dnx

	!nainds
	! list of N for 5 cases
	nalist5 = na + dnalist5
	! max nx possible in 13 cases => nx+1 or nx+2
	nxmax = nx+1;
	if(crosshops) nxmax = nx+2;

	! indexes pointers for basis set sectors with diff no of up spins
	do i=1,5
		m1 = min(nalist5(i),nxmax); ! max up spin possible
		! indexed pointers
		if(allocated(basis(i)%pntr)) deallocate(basis(i)%pntr)
		allocate(basis(i)%pntr(m1+2))
		call pointerslist(nalist5(i),m1,basis(i)%pntr)
		!write(*,*) "i =",i,", basis(i)%pntr = ", basis(i)%pntr

		! subsets
		ntot = basis(i)%pntr(m1+2); ! configurations of this type
		if(allocated(basis(i)%sec)) deallocate(basis(i)%sec)
		allocate(basis(i)%sec(ntot))
		!	ignore the 0-up (all down) sector (a single basis state)
		!	start j=1-m1 for sectors with j up spins (2nd onwards)
		do j=1,m1,1	
			ntot = basis(i)%pntr(j+2) - basis(i)%pntr(j+1); ! for config of this type
			if(allocated(basis(i)%sec(j)%sets)) deallocate(basis(i)%sec(j)%sets)
			allocate(basis(i)%sec(j)%sets(ntot,j))
			call mksets(n,j,ntot,basis(i)%sec(j)%sets(ntot,j))
		end do
	end do

	end subroutine mkbasis


!------------------------------------------

	subroutine mksets(n,k,ntot,combset)
	! genrates subsets
	implicit none
	integer, intent (in) :: n,k,ntot
	integer(kind=1), dimension(ntot,k),intent(out):: combset
	integer(kind=1), dimension(k)::comb
	integer :: ind
	
	ind = 1; combset(:,:)=0;
	call subsets(n,k,1)

	contains
	recursive subroutine subsets(n,k,r)
	! genrates subsets
	implicit none
	integer, intent (in) :: r,n,k
	integer :: l

	if (r > k) then
	 !write (*,*) comb
	 combset(ind,:) = comb;
	 ind = ind + 1;
	else
	 do l = 1, n
		if ((r == 1) .or. (l > comb(r - 1))) then
			comb(r) = l
			call subsets (n,k,r + 1)
		end if
	 end do
	end if

	end subroutine subsets
	end subroutine mksets

!------------------------------------------
	subroutine pointerslist(n,k,pntr)
	! pointer indexed to the subsets with different number of up spins
	! pntr = Insert[Table[Sum[Bimonial[n, j], {j, 0, i}], {i, 0, k}], 0, 1];
	integer*4, intent(in) :: n,k
	integer*4, intent(out) :: pntr(k+2)
	integer i
	
	pntr(1) = 0; pntr(2) = 1;
	do i=3,k+2,1
		pntr(i)=pntr(i-1)+ nCr(n,i-2)
		!write(*,*) "k,i, pntr(i) = ",k,i, pntr(i)
	end do

	return
	end subroutine pointerslist
!------------------------------------------
	integer function nCr (n,r)
	! ^nC_r
	integer n,r,i,large,small
	integer*8 Numer,Denom 

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
	integer, intent(in) :: n,k
	integer(kind=1), dimension(k), intent(in)::list
	integer :: p,np,kp

	
	!x1 = c[N, k] - Sum[c[N - listf[[p + 1]], k - p], {p, 0, k - 1}]; 
	LexicoIndex = nCr(n,k);
	do p=0,k-1
		np = n - list(p+1);
		kp = k -p;
		LexicoIndex = LexicoIndex - nCr(np,kp)
	end do
	return
	end function
!------------------------------------------
	integer(kind=4) function LexicoIndexSort(list,n,k)
	use lists, only: Sort
	! using LexicoIndex() will avoid dependence on module lists.
	!	BETTER TO REMOVE THIS FUNCTION AND USE LexicoIndex() INSTEAD.
	implicit none
	integer, intent(in) :: n,k
	integer(kind=1), dimension(k), intent(in)::list
	integer :: p,np,kp
	! sort the list first
	call Sort(list,k)
	!x1 = c[N, k] - Sum[c[N - listf[[p + 1]], k - p], {p, 0, k - 1}]; 
	LexicoIndex = nCr(n,k);
	do p=0,k-1
		np = n - list(p+1);
		kp = k -p;
		LexicoIndex = LexicoIndex - nCr(np,kp)
	end do
	return
	end function
!------------------------------------------





!------------------------------------------

	!end program
	end module
