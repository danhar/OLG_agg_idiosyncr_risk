module markov_station_distr
    implicit none
contains

pure function f_markov_statdist(pi,nt) result(statdist)
! calculates the stationary distribution of a markov chain by largest eigenvalue,
! or by simulation if optional number of steps nt is given.

    use kinds
!    use evesf_int	! IMSL Math.pdf, p. 563f: Extreme eigenvalues and eigenvectors of real symmetric matrix

    real(dp) ,dimension(:,:) ,intent(in) :: pi	      ! Markov chain
    integer  ,optional       ,intent(in) :: nt 	      ! optional number of simulation steps
    real(dp) ,dimension(size(pi,1))      :: statdist
	real(dp) ,dimension(1)				 :: eigenval  ! return from evesf
	real(dp) ,dimension(size(pi,1),1)	 :: eigenvec
	real(dp) ,dimension(size(pi,1),size(pi,2)) :: identity
	integer							     :: i, tc

	if(.not.present(nt)) then
	! calculate by largest eigenvalue, NEED TO CHECK!!!
		identity=0.0
		forall(i=1:size(identity,1)) identity(i,i)=1.0
!		call evesf(1,(identity-transpose(pi)),.false.,eigenval,eigenvec)
		if ((abs(eigenval(1))-1.0) > 0.0001_dp) return			! Need to change!
		statdist=eigenvec(:,1)/sum(eigenvec(:,1))

	else
	! calculate by iterating on transition matrix
	    statdist=1.0/real(size(pi,1),dp)
	    do tc=1,nt
	    	statdist=matmul(transpose(pi), statdist)
	    end do
	endif

end function f_markov_statdist
end module markov_station_distr
