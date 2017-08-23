
!Dropp[Ds_, elem_] := Complement[Ds, {elem}];
!pos = DPhiL[[ws]];
!Ds = Dropp[Ds, pos];

	program dropelement
	implicit none

	integer, dimension(10):: a	
	integer, dimension(9):: b	
	a=[1,2,3,4,5,10,9,8,7,6];

	write(*,*) "---------------"
 	call drop(a,10,8,b)
 	write(*,*) "drop 8 from ",a
 	write(*,*) " we get ",b
	write(*,*) "---------------"
 	call append(b,9,100,a)
 	write(*,*) "append 100 to b "
 	write(*,*) " we get ",a
	write(*,*) "---------------"
 	call substitute(a,10,100,8)
 	write(*,*) "substitute 100 by 8 in b"
 	write(*,*) " we get ",a
	write(*,*) "---------------"



	contains
!-------------------------------------
	subroutine drop(a,la,x,b)
	implicit none
	integer, intent(in):: x,la
	integer, dimension(la), intent(in):: a
	integer, dimension(la-1), intent(out):: b
	integer i,j,k
	j = 1;
	do i=1,la
		if (a(i) .ne. x) b(j)=a(i); j = j + 1;
	end do
	return
	end subroutine
!-------------------------------------
	subroutine append(a,la,x,b)
	implicit none
	integer, intent(in):: x,la
	integer, dimension(la), intent(in):: a
	integer, dimension(la+1), intent(out):: b
	b(1:la) = a(:); b(la+1) = x;
	return
	end subroutine
!-------------------------------------
	subroutine substitute(a,la,x,y)
	implicit none
	integer, intent(in):: x,la,y
	integer, dimension(la), intent(inout):: a
	integer i
	do i=1,la
		if (a(i) == x) then
			a(i)=y;
			exit
		endif
	end do
	return
	end subroutine
!-------------------------------------



	end program
	
