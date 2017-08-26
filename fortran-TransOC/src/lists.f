	module lists
	! some functions to manipulate lists (1-d arrays)
	implicit none

	contains
!-------------------------------------
	subroutine Drop(a,la,x,b)
	implicit none
	integer(kind=1), intent(in):: x,la
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(la-1), intent(out):: b
	integer(kind=1):: i,j
	j = 1;
	do i=1,la
		if (a(i) .ne. x) then
			b(j)=a(i);
			j = j + 1;
		endif
	end do
	return
	end subroutine
!-------------------------------------
	subroutine Append(a,la,x,b)
	implicit none
	integer(kind=1), intent(in):: x,la
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(la+1), intent(out):: b
	b(1:la) = a(:); b(la+1) = x;
	return
	end subroutine
!-------------------------------------
	subroutine SortedInsert(a,la,x,b)
	! only use when a is already sorted
	implicit none
	integer(kind=1), intent(in):: x,la
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(la+1), intent(out):: b
	integer(kind=1) :: ix,j,i

	j = 1;
	do i=1,la
		if (a(i) .lt. x) then
			b(j)=a(i); j = j + 1;
		else
			b(j) = x; ix = i;
			exit
		endif		
	end do

	if( j > la ) then ! a(:) > x and x not inserted in above loop
		b(j) = x
	else
		b(j:la) = a(ix:la) ! x was inserted, the rest of a here
	endif
	
	return
	end subroutine
!-------------------------------------
	subroutine Substitute(a,la,x,y)
	implicit none
	integer(kind=1), intent(in):: x,la,y
	integer(kind=1), dimension(la), intent(inout):: a
	integer(kind=1):: i
	do i=1,la
		if (a(i) == x) then
			a(i)=y;
			exit
		endif
	end do
	return
	end subroutine
!-------------------------------------
	subroutine Complement0(a,la,b,lb,c,lc)
	!	output array size unknown to calling routine
	implicit none
	integer(kind=1), intent(in):: la,lb
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(lb), intent(in):: b
	integer(kind=1), intent(out) :: lc
	integer(kind=1), dimension(la), intent(out):: c ! dim=la, calling routine will fix the size later
	integer i,j
	logical inb

	lc=0;
	do i=1,la
		inb = .false.
		do j=1,lb	
		if (a(i) == b(j)) then
			inb = .true.
			exit
		endif
		end do
		if(.not. inb) lc = lc+1; c(lc) = a(i); 
	end do

	return
	end subroutine
!-------------------------------------
	subroutine Complement(a,la,b,lb,c)
	!	b is a subset of a: output array size known lc = la-lb
	implicit none
	integer(kind=1), intent(in):: la,lb
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(lb), intent(in):: b
	integer(kind=1), dimension(la-lb), intent(out):: c
	integer:: i,j,lc
	logical:: inb

	!write(*,*) "a=",a
	!write(*,*) "la-lb, b=", la-lb, b

	lc=0;
	do i=1,la
	
		inb = .false.
		
		do j=1,lb	
			if (a(i) == b(j)) then
				inb = .true.
				exit
			endif
		end do
		
		if(.not. inb) then
			lc = lc+1;
			c(lc) = a(i);
		endif

	end do

	!write(*,*) "====> c=",c
	return
	end subroutine
!-------------------------------------
	logical function MemberQ(a,la,x)
	implicit none
	integer(kind=1), intent(in):: la,x
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1):: i
	
	MemberQ = .false.
	do i=1,la
		if (a(i) == x) then
			MemberQ = .true.
			exit
		endif
	end do
	return
	end function
!-------------------------------------
	logical function FreeQ(a,la,x)
	implicit none
	integer(kind=1), intent(in):: la,x
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1):: i
	
	FreeQ = .true.
	do i=1,la
		if (a(i) == x) then
			FreeQ = .false.
			exit
		endif
	end do
	return
	end function
!-------------------------------------
	subroutine Intersection(a,la,b,lb,c,lc)
	!	output array size unknown to calling routine
	implicit none
	integer(kind=1), intent(in):: la,lb
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(lb), intent(in):: b
	integer(kind=1), intent(out) :: lc
	integer(kind=1), dimension(la), intent(out):: c ! dim=la, calling routine will fix the size later
	integer i,j
	logical inb

	lc=0;
	do i=1,la
		inb = .false.
		do j=1,lb	
		if (a(i) == b(j)) then
			inb = .true.
			exit
		endif
		end do
		if(inb) then
			lc = lc+1;
			c(lc) = a(i);
		endif 
	end do
	
	return
	end subroutine
!-------------------------------------
	subroutine Join(a,la,b,lb,c)
	!	output array size known lc = la+lb
	implicit none
	integer(kind=1), intent(in):: la,lb
	integer(kind=1), dimension(la), intent(in):: a
	integer(kind=1), dimension(lb), intent(in):: b
	integer(kind=1), dimension(la+lb), intent(out):: c

	c(1:la) = a(:);
	c(la+1:la+lb) = b(:);
	
	return
	end subroutine
!-------------------------------------
 ! a modified version of https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
	recursive subroutine Sort(a,la)
 	implicit none
	integer(kind=1), intent(in):: la
	integer(kind=1), dimension(la), intent(inout):: a
 
	! LOCAL VARIABLES
	integer(kind=1) :: left, right
	real :: random
	real :: pivot
	integer(kind=1) :: temp
	integer(kind=1) :: marker,marker2
 
	if (la > 1) then
 
		call random_number(random)
		pivot = A(int(random*real(la-1))+1)   ! random pivor (not best performance, but avoids worst-case)
		left = 0
		right = la + 1
 
		do while (left < right)
			right = right - 1
			do while (A(right) > pivot)
				right = right - 1
			end do
			left = left + 1
			do while (A(left) < pivot)
				left = left + 1
			end do
			if (left < right) then
				temp = A(left)
				A(left) = A(right)
				A(right) = temp
			end if
	end do
 
	if (left == right) then
		marker = left + 1
	else
		marker = left
	end if

	marker2 = marker-1;
	call Sort(A(:marker2),marker2);
	marker2 = la-marker+1;
	call Sort(A(marker:),marker2)
 
	end if
	end subroutine Sort
!-------------------------------------

	end module
	
