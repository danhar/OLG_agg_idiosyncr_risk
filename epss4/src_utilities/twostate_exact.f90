module twostate_exact
    use kinds
    implicit none
contains

pure subroutine calibridiorisk2(rho,sigma,transition,states,stat_dist, converged)
! Calculate an 'exact 2-state-approximation' to an AR(1)-process of the form
! y(t) = rho y(t-1)+ e(t),   e(t)~N(sigma,1)
! See Harenberg_JMpaper.pdf, App. B3 for the formulae below

    use fun_zbrent
    use sub_zbrac
    use markov_station_distr

    real(dp)                  ,intent(in)  :: rho, sigma(:) ! correlation coefficient and standard deviation of innovation term
    real(dp) ,dimension(2,2)  ,intent(out) :: transition    ! Markov chain transition probabilities
    real(dp) ,dimension(2,size(sigma)) ,intent(out) :: states ! symmetric and evenly spaced Markov state spabe
    real(dp) ,dimension(2)    ,intent(out) :: stat_dist     ! stationary distribution
    logical                   ,intent(out) :: converged     ! true if solution found, else false
    real(dp), parameter                    :: tol= 1.0e-06
    real(dp)                               :: pi1, brackl, bracku
    real(dp), dimension(size(sigma))       :: sig_uncond, pheye, tildsig
    integer                                :: nz, zc

    converged = .true.
    nz = size(sigma)

    if (any(sigma == 0.0)) then
        states     = 1.0
        pi1        = (1+rho)/2.0 ! for lack of a better guess, I take the approximate value
        transition = reshape([    pi1, 1.0-pi1, &
                              1.0-pi1,     pi1] ,shape(transition), order=[2,1])
        stat_dist  = 0.5_dp
        return
    endif

	! unconditional std
	sig_uncond = sigma/sqrt((1.0-rho**2))

	do zc =1, nz
	    brackl=0.5*sig_uncond(zc)
	    bracku=1.5*sig_uncond(zc)
	    call s_zbrac(func_tildsig,brackl,bracku, converged)
	    if (.not. converged) return
	    tildsig(zc) = f_zbrent(func_tildsig,brackl,bracku)
	    pheye(zc)=calc_pheye(tildsig(zc))
	    if (sign(func_tildsig(tildsig(zc)),1.0)>tol) converged = .false.
	    states(:,zc)=ettav(tildsig(zc))
	enddo
	! entries of transition matrix:
	pi1=(1.0+rho)/2.0
	pi1=f_zbrent(func_pi1,0.0_dp,1.0_dp)
	if (sign(func_pi1(pi1),1.0)>tol) converged = .false.

	transition = reshape([    pi1, 1.0-pi1, &
	                      1.0-pi1,     pi1] ,shape(transition), order=[2,1])

    stat_dist = f_markov_statdist(transition,5000)
contains

	pure real(dp) function func_tildsig(tildsig) result(df)
	    real(dp), intent(in) :: tildsig
	    df=calc_pheye(tildsig)**2+tildsig**2-sig_uncond(zc)**2
	end function func_tildsig
	!----------------------------------------------------------------------------------
	pure real(dp) function calc_pheye(tildsig)
	   real(dp), intent(in) :: tildsig
	    real(dp), dimension(2) :: temp
        temp(1)=1.0-tildsig
        temp(2)=1.0+tildsig
        temp=exp(temp)
        calc_pheye=log(2.0/sum(temp))+1.0
	end function calc_pheye

	!----------------------------------------------------------------------------------
	pure real(dp) function func_pi1(pi)
	real(dp), intent(in) :: pi
	real(dp) :: etemp, temp(2)
	integer :: zc
	temp=0.0
	do zc=1,nz
	    temp(1)=temp(1)+1.0/real(nz,dp)*(pi*(pheye(zc)-tildsig(zc))**2.0+&
	        (1.0-pi)*(pheye(zc)-tildsig(zc))*(pheye(zc)+tildsig(zc)))
	    temp(2)=temp(2)+1.0/real(nz,dp)*(pi*(pheye(zc)+tildsig(zc))**2.0+&
	        (1.0-pi)*(pheye(zc)+tildsig(zc))*(pheye(zc)-tildsig(zc)))
    enddo
	etemp=sum(temp)/2.0

	func_pi1=etemp/(sum(sig_uncond**2.0)/nz)-rho
	end function func_pi1

end subroutine calibridiorisk2
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

pure function ettav(msig)
	real(dp), dimension(2) :: ettav
	real(dp), intent(in) :: msig
	real(dp) :: temp
	temp=exp(1.0-msig)+exp(1.0+msig)
	ettav(1)=1.0-msig
	ettav(2)=1.0+msig
	ettav=2.0*exp(ettav)/temp
end function ettav

end module twostate_exact
