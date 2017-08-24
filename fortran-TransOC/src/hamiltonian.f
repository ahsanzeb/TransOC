
	module hamiltonian
	use modmain
	implicit none




!------------------------------------------
	subroutine HgMap(ib,n,k,ntot,map)
	use modmain, only: basis, sites,na
	use lists, only: Complement,Append,Sort,LexicoIndex
	implicit none
	integer, intent(in) :: ib,n,k,ntot
	integer(kind=4), dimension(ntot,n-k), intent(out):: map
	! local
	integer:: nak,flip,i,j
	integer(kind=1), dimension(k):: set1
	integer(kind=1), dimension(k+1):: set2
	integer(kind=1), dimension(n-k):: flipsites
		
	if (k==0) then
		map(1,:) = (/ (i,i=1,n,1)/)
	elseif(k==n) then
		! all map to the single state with all up
		map(:,1) = 1;
	else
		nak = na-k;
		do i=1,ntot
		! sec(k) are defined for k>0,
		! so index match with k
			set1=basis(ib)%sec(k)%sets(i,:)
			call Complement(sites,na,set1,k,flipsites)
			do j=1,nak
				flip = flipsites(j);
				call Append(set1,k,flip,set2)
				call Sort(set2,k+1);
				map(i,j) = LexicoIndex(set2,n,k+1)
			end do
		end do
	endif
	return
	end subroutine HgMap
!------------------------------------------




	end module


	
