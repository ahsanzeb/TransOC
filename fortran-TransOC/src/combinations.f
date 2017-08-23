	!program combinations
	module combinations
	use modmain
	
	implicit none

	contains
!------------------------------------------
	recursive subroutine subsets (n,k,r)
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

	!end program
	end module
