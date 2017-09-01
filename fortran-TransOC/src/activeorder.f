!---------------------------------------
	subroutine ActiveOrder(ASites,na,l,l1,l2,order)
	implicit none
	integer(kind=1), intent(in) :: na,l
	integer, dimension(na), intent(in):: Asites
	integer(kind=1), intent(out) :: l1,l2
	logical, intent(out) :: order
	! local
	integer(kind=1) :: i,x1,x2,lp1
	

	lp1 = l + 1;
	! order in the list of active sites
	x1 = GetPosition(l);
	x2 = GetPosition(lp1);
	if (x1<x2) then
		l1=x1; l2=x2; order=.true.;
	else
		l2=x1; l1=x2; order=.false.;
	endif

	return
	
	contains
	!----------------------------------
	integer function GetPosition(x)
	integer(kind=1), intent(in):: x
	GetPosition = 0;
	do i=1,na
		if(l == ASites(i)) then
			GetPosition = i;
			exit
		endif
	end do
	if(GetPosition==0) then
		write(*,*) "GetPosition: element not found!"
		stop
	endif
	return
	end function GetPosition
	!----------------------------------
	end subroutine ActiveOrder
!---------------------------------------

