
	program dropelement

	use modmain
	use combinations
	
	implicit none


	type(BasisSet), dimension(5) basis
		!integer, allocatable :: pntr(:)
		!integer, allocatable :: sets(:)


	basis

	do i=0,n
		if(allocated(pntr)) deallocate(pntr)
		allocate(pntr(i+2))
		call pointerslist(n,i,pntr)
		write(*,*) "i =",i,", pntr = ", pntr
	end do

	allocate(combset(100,5))
	k = min(3,n);
	!allocate(comb(k))
	ind = 1; combset(:,:)=0;
	call subsets(n,k,1)
	write (*, *)"====================="
	do i=1,20
	 write (*, *) combset(i,:)
	end do




	end program
	
