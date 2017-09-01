
	module rates
	use modmain, only: qt,eig,itypes,mapb,mapt,maph,mapc
	implicit none

	public:: 
	private::

	contains



	! Right and left contact energy barriers
	{Ebl, Ebr} = Eb;
	! classical + quantum changes
	DCQ = DEQCsum[w0, Ebr, Ebl, Er];

	





	!---------------------------------------
	double precision function penalty(de)
	use modmain, only: beta
	implicit none
		penalty = 1.0d0;
		if(de > 0.0d0 ) penalty = dexp(-de*beta);
	return
	end function penalty
	!---------------------------------------
	end module rates













	subroutine CalRates(ih,ic,is,r,nr)
	use modmain, only: qt,eig,itypes,psi,mapt,maph
	implicit none

	integer, intent(in) :: ih,ic,is
	double precision, intent(out) :: r 
	! local
	integer(kind=4) :: nc,ns 		
	double precision :: r 
	double precision, allocatable :: rcs(:,:)

	rates(ih)%nc
	rates(ih)%ns
	rates(ih)%r
	rates(ih)%rcs()

	qt(ia)%cs(ic,is)%namp != n2
	qt(ia)%cs(ic,is)%amp
	
	qt(ia)%cs(ic,is)%nsec
	qt(ia)%cs(ic,is)%amp2


	rs = 0.0d0
	
	do i=1,nrs
		rs(i) = penalty(de) * amp2
	end do
	end subroutine CalRates









