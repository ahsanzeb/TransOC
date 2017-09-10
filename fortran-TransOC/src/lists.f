	module lists
	use modmain, only: isk
	! some functions to manipulate lists (1-d arrays)
	implicit none

	!public::
	!private::

	contains
!-------------------------------------
	subroutine Drop(a,la,x,b)
	implicit none
	integer, intent(in):: x,la
	integer(kind=isk), dimension(la), intent(in):: a
	integer(kind=isk), dimension(la-1), intent(out):: b
	integer:: i,j
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
	subroutine Drop4(a,la,x,b)
	implicit none
	integer, intent(in):: x,la
	integer, dimension(la), intent(in):: a
	integer, dimension(la-1), intent(out):: b
	integer:: i,j
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
	subroutine SortedInsert(a,la,x,b)
	! only use when a is already sorted
	implicit none
	integer, intent(in):: x,la
	integer(kind=isk), dimension(la), intent(in):: a
	integer(kind=isk), dimension(la+1), intent(out):: b
	integer:: ix,j,i

	j = 1;
	do i=1,la
		if (x > a(i)) then
			b(j)=a(i); j = j + 1;
		else
			b(j) = x; ix = i;
			exit
		endif		
	end do

	if( j > la ) then ! a(:) > x and x not inserted in above loop
		b(j) = x
	else
		b(j+1:la+1) = a(ix:la) ! x was inserted, the rest of a here
	endif
		!write(*,*)"a, x=",a,x
		!write(*,*)"b=",b
			
	return
	end subroutine
!-------------------------------------
	subroutine Substitute(a,la,x,y)
	implicit none
	integer, intent(in):: x,la,y
	integer, dimension(la), intent(inout):: a
	integer:: i
	do i=1,la
		if (a(i) == x) then
			a(i)=y;
			exit
		endif
	end do
	return
	end subroutine
!-------------------------------------
	subroutine Complement(a,la,b,lb,c)
	!	b is a subset of a: output array size known lc = la-lb
	implicit none
	integer, intent(in):: la,lb
	integer(kind=isk), dimension(la), intent(in):: a
	integer(kind=isk), dimension(lb), intent(in):: b
	integer(kind=isk), dimension(la-lb), intent(out):: c
	integer:: i,j,lc
	logical:: inb

	!write(*,*) "a=",a
	!write(*,*) "b=",b
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
	integer, intent(in):: la,x
	integer(kind=isk), dimension(la), intent(in):: a
	integer:: i
	
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
	logical function MemberQ4(a,la,x)
	implicit none
	integer, intent(in):: la,x
	integer, dimension(la), intent(in):: a
	integer:: i
	
	MemberQ4 = .false.
	do i=1,la
		if (a(i) == x) then
			MemberQ4 = .true.
			exit
		endif
	end do
	return
	end function
!-------------------------------------
	logical function FreeQ(a,la,x)
	implicit none
	integer, intent(in):: la,x
	integer, dimension(la), intent(in):: a
	integer:: i
	
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
	subroutine Join(a,la,b,lb,c)
	!	output array size known lc = la+lb
	implicit none
	integer, intent(in):: la,lb
	integer, dimension(la), intent(in):: a
	integer, dimension(lb), intent(in):: b
	integer, dimension(la+lb), intent(out):: c

	c(1:la) = a(:);
	c(la+1:la+lb) = b(:);
	
	return
	end subroutine
!-------------------------------------
	integer function GetPosition(Asites,na,l)
	integer, intent(in):: na,l
	integer, dimension(na), intent(in):: Asites
	! local
	integer :: i

	GetPosition = 0;
	do i=1,na
		if(l == ASites(i)) then
			GetPosition = i;
			exit
		endif
	end do
	if(GetPosition==0) then
		write(*,*) "GetPosition: element not found!"
		write(*,*) " input array =  ",Asites
		write(*,*) "length, finding position of element = ",na,l		
		stop
	endif
	return
	end function GetPosition
!----------------------------------









	end module
	
