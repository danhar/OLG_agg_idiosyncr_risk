module meanshock_wrapper
    implicit none
contains

subroutine SolveMeanShock(coeffs, grids, policies, simvars, lifecycles, Phi, value, err, output_path)

    use kinds
    use types
    use aggregate_grids ,only: tAggGrids
    use laws_of_motion  ,only: tCoeffs, Initialize
    use error_class
    use params_mod      ,only: n_coeffs,nj,nx,nz,n_eta,alpha,etagrid,stat_dist_z, partial_equilibrium
    use income
	use meanshock_equilib ,only: ms_equilib_set, ms_equilib_get, ms_equilib
	use sub_broyden
    use fun_locate
    use fun_aggregate_diff

    type(tCoeffs)   ,intent(out)    :: coeffs       ! coefficients for laws of motion
    type(tAggGrids) ,intent(inout) :: grids        ! grids for aggregate states k and mu
    type(tPolicies) ,intent(out)   :: policies
    type(tSimvars)  ,intent(out)   :: simvars
    type(tLifecycle),intent(out)   :: lifecycles
    real(dp) ,allocatable ,intent(out)   :: Phi(:,:,:), value(:,:,:,:,:,:)  ! distribution, valuefunction
    type(tErrors)   ,intent(out)   :: err
    character(len=*), intent(in)   :: output_path
	real(dp) ,dimension(2)		   :: xvars, fvals	! input/output of ms_equilib, passed to sub_broyden
	real(dp) ,dimension(n_eta)     :: m_etagrid
	real(dp)                       :: mean_zeta, mean_delta ! mean shocks
	real(dp)					   :: wz, wd, w(nz)	! weights (distance to mean zeta, mean delta)
	real(dp) ,dimension(nx,n_eta,nj) :: apgrid_ms, stock_ms, kappa_ms, xgrid_ms, value_ms ! mean shock projections
	integer         	           :: i		        ! index

    coeffs          = Initialize('msge',n_coeffs,nz)
	mean_zeta	    = dot_product(stat_dist_z, zeta)
	mean_delta		= dot_product(stat_dist_z, delta)
	i				= f_locate(zeta,mean_zeta)	  ! In 'default', f_locate returns ju-1 if x>xgrid(ju-1) !
	if (zeta(i+1)==zeta(i)) then   ! this happens e.g. if scale_AR = -1.0
        wz          = 0.5_dp
    else
	    wz			= (mean_zeta-zeta(i))/(zeta(i+1)-zeta(i))
    endif
	i				= f_locate(delta,mean_delta)
	!wd				= (mean_deltam-delta(i))/(delta(i+1)-delta(i))
	wd				= 0.5_dp ! WARNING: overwriting wd manually, coz delta is not in increasing order
	w(1)			= (1.0-wd)*(1.0-wz)
	w(2)			= wd*(1.0-wz)
	w(3)			= (1.0-wd)*wz
	w(4)			= wd*wz
	m_etagrid	    = wz*etagrid(:,1)+(1-wz)*etagrid(:,nz)

    call ms_equilib_set(coeffs, grids, mean_zeta, mean_delta, w, m_etagrid)

    ! Initial guesses
    xvars(1)   = grids%K (1)
    xvars(2)   = grids%mu(1)

	if (partial_equilibrium) then
	    fvals = ms_equilib(xvars)
	else ! Find capital and mu for the mean shock general equilibrium
		call s_broyden(ms_equilib,xvars,fvals,err%not_converged, get_fd_jac_o=.true.,tolf_o=1e-2_dp)  ! ,maxstp_o=.8_dp,tolf_o=1e-10_dp
		if (err%not_converged) call err%write2file(fvals, output_path)
    endif

	grids%K (1) = xvars(1)
	grids%mu(1) = xvars(2)
    call ms_equilib_get(Phi, policies, value, err)
    call err%print2stderr

    ! Calc mean shock policies for lc profiles and for saving
    xgrid_ms =0.0; apgrid_ms =0.0; stock_ms =0.0; value_ms = 0.0
    do i=1,nz
        xgrid_ms    = xgrid_ms  + w(i)* policies%xgrid(:,:,i,:,1,1)
        apgrid_ms   = apgrid_ms + w(i)* policies%apgrid(:,:,i,:,1,1)
        stock_ms    = stock_ms  + w(i)* policies%apgrid(:,:,i,:,1,1)*policies%kappa(:,:,i,:,1,1)
        value_ms    = value_ms  + w(i)* value(:,:,i,:,1,1)
    enddo
    where (apgrid_ms .ne. 0.0)
        kappa_ms = stock_ms/apgrid_ms
    elsewhere
        kappa_ms = 0.0
    end where

    call simulate_ms(simvars)

    call ms_lc_profiles(lifecycles)

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - pure subroutine simulate_ms(simvars)
! - pure subroutine ms_lc_profiles(lifecycles)
!-------------------------------------------------------------------------------
    pure subroutine simulate_ms(simvars)
    ! Variable simvars is used to get a very rough approximation of standard deviations.
    ! At index 1, store mean shock values of aggregates.
    ! At index i+1, store aggregates for z=i, calculated by assuming:
    ! - agents live in a world where z(i) always realizes, use LOMs of mean shock (i.e. kp =k, mup = mu)
    ! - Phi remains constant at Phi_ms

        use params_mod   ,only: L_N_ratio, n, g, pi_z, de_ratio
        use fun_aggregate_diff
        use distribution ,only: CheckPhi
        use partial_sorting     ! function valnth

        type(tSimvars) ,intent(out) :: simvars
        real(dp) ,dimension(nx,n_eta,nj):: r_pf
        integer                     :: i

        call AllocateType(simvars,nz+1)   ! allocate all simulation variables
        simvars%z(1)    = 0     ! to indicate that these are mean-shock-results
        simvars%K(1)    = xvars(1)
        simvars%mu(1)   = xvars(2)

        simvars%output(1)= mean_zeta*simvars%K(1)**alpha
        simvars%stock(1) = sum(stock_ms*Phi)/ L_N_ratio ! different from K(t+1) since it is in today's per capita terms. Thus no (1+g)(1+n) in denominator.
        simvars%bonds(1) = sum((apgrid_ms-stock_ms)*Phi)/ L_N_ratio
        simvars%invest(1)= simvars%stock(1) + simvars%bonds(1) -simvars%K(1)*(1.0-mean_delta)
        simvars%C(1)    = sum((xgrid_ms-apgrid_ms)*Phi)/L_N_ratio       ! Consumption per worker
        simvars%rf(1)   = f_riskfree_rate(simvars%K(1),simvars%mu(1),stat_dist_z)
        simvars%r(1)    = f_stock_return(simvars%K(1), mean_zeta, mean_delta, simvars%rf(1))
        r_pf = sign(1.0,apgrid_ms)*(simvars%rf(1) + kappa_ms*simvars%mu(1))/(1.0+g)
        simvars%r_pf_median(1) = valnth(pack(r_pf, Phi/=0.0), ceiling(size(pack(r_pf, Phi/=0.0))/2.0))
        ! The next calculation is neglecting sign(1.0,apgridt), but that would become unnecessarily tedious
        simvars%r_pf_kappa_med(1)=(simvars%rf(1) + valnth(pack(kappa_ms,Phi/=0.0), ceiling(size(pack(kappa_ms, Phi/=0.0))/2.0)) *simvars%mu(1))/(1.0+g)
        simvars%wage(1) = f_netwage (simvars%K(1), mean_zeta)
        simvars%pens(1) = f_pensions(simvars%K(1), mean_zeta)
        simvars%tau(1)  = f_tau     (simvars%K(1), mean_zeta)
        simvars%welf(1) = sum(value_ms(:,:,1)*Phi(:,:,1))

        call CheckPhi(Phi, simvars%Phi_1(1), simvars%Phi_nx(1))
        simvars%bequests(1)   = f_bequests(simvars%rf(1), simvars%r(1), kappa_ms, apgrid_ms, Phi)
        simvars%err_aggr(1)   = f_aggregate_diff(simvars%output(1), simvars%invest(1), simvars%C(1), simvars%bequests(1))
        simvars%B(1)          = fvals(2)
        simvars%err_income(1) = f_income_diff(simvars%K(1), mean_zeta, simvars%r(1), simvars%rf(1), mean_delta)

        do i= 1,nz  ! can't I just save z and then call simulate_economy?
            simvars%z(i+1)     = i
            simvars%K(i+1)     = sum(policies%apgrid(:,:,i,:,1,1)*Phi)/(L_N_ratio*(1.0+n)*(1.0+g))
            simvars%mu(i+1)    = simvars%mu(1)            ! Keep mu constant, could keep rf constant instead
	        simvars%output(i+1)= zeta(i)*simvars%K(i+1)**alpha
	        simvars%stock(i+1) = sum((policies%kappa(:,:,i,:,1,1)*policies%apgrid(:,:,i,:,1,1))*Phi)/ L_N_ratio ! different from K(t+1) since it is in today's per capita terms. Thus no (1+g)(1+n) in denominator.
	        simvars%bonds(i+1) = sum(((1.0-policies%kappa(:,:,i,:,1,1))*policies%apgrid(:,:,i,:,1,1))*Phi)/ L_N_ratio
	        simvars%invest(i+1)   = simvars%stock(i+1) + simvars%bonds(i+1) -simvars%K(i+1)*(1.0-delta(i))
            simvars%C(i+1)     = sum((policies%xgrid(:,:,i,:,1,1)-policies%apgrid(:,:,i,:,1,1)) *Phi)/ L_N_ratio
            simvars%rf(i+1)    = f_riskfree_rate(simvars%K(i+1),simvars%mu(i+1),pi_z(i,:))
            simvars%r(i+1)     = f_stock_return (simvars%K(i+1), zeta(i), delta(i), simvars%rf(i+1))
            r_pf = sign(1.0,policies%apgrid(:,:,i,:,1,1))*(simvars%rf(i+1) + policies%kappa(:,:,i,:,1,1)*simvars%mu(i+1))/(1.0+g)
            simvars%r_pf_median(i+1) = valnth(pack(r_pf, Phi/=0.0), ceiling(size(pack(r_pf, Phi/=0.0))/2.0))
            ! The next calculation is neglecting sign(1.0,apgridt), but that would become unnecessarily tedious
            simvars%r_pf_kappa_med(i+1)=(simvars%rf(i+1) + valnth(pack(policies%kappa(:,:,i,:,1,1),Phi/=0.0), ceiling(size(pack(policies%kappa(:,:,i,:,1,1), Phi/=0.0))/2.0)) *simvars%mu(i+1))/(1.0+g)
	        simvars%wage(i+1)  = f_netwage (simvars%K(i+1), zeta(i))
	        simvars%pens(i+1)  = f_pensions(simvars%K(i+1), zeta(i))
	        simvars%tau (i+1)  = f_tau     (simvars%K(i+1), zeta(i))
	        simvars%welf(i+1)  = sum(value(:,:,i,1,1,1)*Phi(:,:,1))

            simvars%Phi_1(i+1)      = simvars%Phi_1(1)
            simvars%Phi_nx(i+1)     = simvars%Phi_nx(1)
	        simvars%bequests(i+1)   = f_bequests(simvars%rf(i+1), simvars%r(i+1), policies%kappa(:,:,i,:,1,1), policies%apgrid(:,:,i,:,1,1), Phi)
            simvars%err_aggr(i+1)   = f_aggregate_diff(simvars%output(i+1), simvars%invest(i+1), simvars%C(i+1), simvars%bequests(i+1))
	        simvars%B(i+1)          = simvars%K(i+1) * de_ratio/(1.0 + de_ratio) - simvars%bonds(i+1)
	        simvars%err_income(i+1) = f_income_diff(simvars%K(i+1), zeta(i), simvars%r(i+1), simvars%rf(i+1), delta(i))
        enddo

        simvars%K (nz+2) = simvars%K (1) ! K and rf have one more index in the real simulations.
        simvars%rf(nz+2) = simvars%rf(1)
    end subroutine simulate_ms

    pure subroutine ms_lc_profiles(lifecycles)
    ! Calculates the policies in the mean shock and then the corresponding life-cycle profiles
        use params_mod   ,only        :  g
        type(tLifecycle) ,intent(out) :: lifecycles
        integer                       :: jc
        call AllocateType(lifecycles,nj)

        lifecycles%ap      = sum(sum(apgrid_ms * Phi,1),1)
        lifecycles%cons    = sum(sum((xgrid_ms-apgrid_ms) * Phi,1),1)
        lifecycles%stock   = sum(sum(stock_ms * Phi,1),1)
        lifecycles%return  = sum(sum(Phi*apgrid_ms*(1.0 + simvars%rf(1) + kappa_ms*simvars%mu(1))/(1.0+g),1),1)
        do jc=1,nj
            lifecycles%cons_var(jc)  = sum((((xgrid_ms(:,:,jc)-apgrid_ms(:,:,jc)) - lifecycles%cons(jc)))**2 * Phi(:,:,jc))
            lifecycles%return_var(jc)= sum((apgrid_ms(:,:,jc)*(1.0 + simvars%rf(1) + kappa_ms(:,:,jc)*simvars%mu(1))/(1.0+g) - lifecycles%return(jc))**2 * Phi(:,:,jc))
        enddo

        where (lifecycles%ap .ne. 0.0)
            lifecycles%kappa = lifecycles%stock/lifecycles%ap
        elsewhere
            lifecycles%kappa = 0.0
        end where
    end subroutine ms_lc_profiles

end subroutine SolveMeanShock
end module meanshock_wrapper
