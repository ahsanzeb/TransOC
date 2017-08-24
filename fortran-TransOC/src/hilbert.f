

	program xxx
	implicit none
	integer i,j,jj,k
	
	type :: HilbertSpace
		integer :: ntot
		integer, allocatable :: pntr(:)
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
	end type HilbertSpace




	type(HilbertSpace), dimension(13):: hs

	do i=1,13
		allocate(hs(i)%pntr(i))
		hs(i)%pntr = (/(j**2, j=1,i, 1)/)


		allocate(hs(i)%evec(i,i))
		do k=1,i,1
			hs(i)%evec(:,k) = (/ (j+k, j=1,i, 1)/)
		end do
		!if(allocated(x))deallocate(x)
		!allocate(x(i))
		!x = pntr(i)%dat
		!write(*,*) "----------------------"
		!write(*,*)	hs(i)%pntr
		!write(*,*)	hs(i)%evec
	end do

	!write(*,*)	hs(3)%evec

	
	end program
