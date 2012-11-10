module mean_shock_mod
! This module conains all procedures to solve for the 'mean shock equilibrium' (MSE).
! The MSE is a degenerate stochastic steady state which we compute to get initial estimates
! for aggregate capital k, the equity premium mu, and the distribution Phi.
    implicit none
contains

subroutine solve_meanshock(coeffs, grids, policies, simvars, lifecycles, Phi, xgrid_ms, value, err, output_path)
    ! Set up environment to use rootfinder for means shock k and mean shock mu,
    ! then pass the function ms_equilibrium as a function argument to a root finder.
    ! In this version, the function argument is an internal procedure, which is a thread-safe Fortran 2008 feature implemented
    ! in the Intel Fortran Compiler >= 11.0 and in gfortran >= 4.5

    use kinds
    use types
    use policies_class  ,only: tPolicies
    use aggregate_grids ,only: tAggGrids
    use laws_of_motion  ,only: tCoeffs, Initialize
    use error_class
    use params_mod      ,only: n_coeffs,nj,nx,n_eta,alpha,etagrid,stat_dist_z, partial_equilibrium
    use income
	use sub_broyden
    use fun_locate
    use fun_aggregate_diff

    type(tCoeffs)   ,intent(out)   :: coeffs       ! coefficients for laws of motion
    type(tAggGrids) ,intent(inout) :: grids        ! grids for aggregate states k and mu
    type(tPolicies) ,intent(out)   :: policies
    type(tSimvars)  ,intent(out)   :: simvars
    type(tLifecycle),intent(out)   :: lifecycles
    real(dp) ,allocatable ,intent(out)   :: Phi(:,:,:), value(:,:,:,:,:,:), xgrid_ms(:,:,:)  ! distribution, valuefunction, mean shock xgrid
    type(tErrors)   ,intent(out)   :: err
    character(len=*), intent(in)   :: output_path
	real(dp) ,dimension(2)		   :: xvars, fvals	! input/output of ms_equilib, passed to sub_broyden
	real(dp) ,dimension(:), allocatable  :: m_etagrid, w
	real(dp)                       :: mean_zeta, mean_delta ! mean shocks
	real(dp)					   :: wz, wd ! weights (distance to mean zeta, mean delta)
	real(dp) ,dimension(:,:,:), allocatable :: apgrid_ms, stocks_ms, kappa_ms, value_ms ! mean shock projections
	integer         	           :: i, nz		        ! index

    nz = size(stat_dist_z)
    allocate(w(nz))
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

    ! Initial guesses
    xvars(1)   = grids%k (1)
    xvars(2)   = grids%mu(1)

	if (partial_equilibrium) then
	    fvals = ms_equilibrium(xvars)
	else ! Find capital and mu for the mean shock general equilibrium
		call s_broyden(ms_equilibrium,xvars,fvals,err%not_converged, get_fd_jac_o=.true.,tolf_o=1e-2_dp)  ! ,maxstp_o=.8_dp,tolf_o=1e-10_dp
		if (err%not_converged) call err%write2file(fvals, output_path)
    endif

	grids%k (1) = xvars(1)
	grids%mu(1) = xvars(2)

    call get_equilibrium_values(policies,value,apgrid_ms, stocks_ms, xgrid_ms, kappa_ms, value_ms, Phi, err)
    call err%print2stderr

    call simulate_ms(simvars)

    call ms_lc_profiles(lifecycles)

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - pure function ms_equilibrium(msvars) result(distance)
! - pure subroutine get_equilibrium_values(fine,v_fine,apgrid_ms, stocks_ms, xgrid_ms, kappa_ms, value_ms, Phi, errs)
! - pure subroutine simulate_ms(simvars)
! - pure subroutine ms_lc_profiles(lifecycles)
!-------------------------------------------------------------------------------


    function ms_equilibrium(msvars) result(distance)
    ! Solve for the MSE, i.e. where k'=k and mu'=mu given that all realizations are at the mean_z.
    ! The procedure is would be pure pure but for the OMP directives in olg_backwards_solution (but it does read access host variables).
        use params_mod             ,only: L_N_ratio, n, g, stat_dist_z, de_ratio, nx_factor
        use error_class            ,only: tErrors
        use aggregate_grids        ,only: AllocateType
        use household_solution_mod ,only: olg_backwards_recursion
        use distribution           ,only: TransitionPhi
        use interpolate_xgrid      ,only: InterpolateXgrid
        use income

        implicit none
        real(dp) ,dimension(:),intent(in) :: msvars             ! k_ms, mu_ms
        real(dp) ,dimension(size(msvars)) :: distance           ! excess demands
        type(tPolicies) :: policies, fine
        type(tAggGrids) :: grid
        type(tErrors)   :: errs
        real(dp) ,dimension(:,:,:,:,:,:) ,allocatable :: value, v_fine
        real(dp) ,dimension(:,:,:) ,allocatable :: apgrid_ms, stocks_ms, xgrid_ms, Phi ! mean shock projections and distribution
        real(dp)                          :: netwage_ms, pens_ms, r_ms, rf_ms, kp_ms, agg_bond_demand
        integer                           :: i

        call AllocateType(grid, size(grids%k), size(grids%mu))
        grid%k  = msvars(1)
        grid%mu = msvars(2)

        call olg_backwards_recursion(policies,coeffs, grid, value, errs)

        call InterpolateXgrid(nx_factor, policies, value, fine, v_fine)

        allocate(apgrid_ms(size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)), &
                 stocks_ms(size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)), &
                 xgrid_ms (size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4))  )

        ! Projection of policies / grid on mean shock
        xgrid_ms =0.0; apgrid_ms =0.0; stocks_ms =0.0
        do i=1,size(fine%xgrid,3)
            xgrid_ms  = xgrid_ms  + w(i)* fine%xgrid (:,:,i,:,1,1)
            apgrid_ms = apgrid_ms + w(i)* fine%apgrid(:,:,i,:,1,1)
            stocks_ms = stocks_ms + w(i)* fine%stocks(:,:,i,:,1,1)
        enddo

        ! Prices in mean shock path
        netwage_ms    = f_netwage (grid%k(1), mean_zeta)
        pens_ms       = f_pensions(grid%k(1), mean_zeta)
        rf_ms         = f_riskfree_rate(grid%k(1),grid%mu(1),stat_dist_z)
        r_ms          = f_stock_return(grid%k(1), mean_zeta, mean_delta, rf_ms)

        ! Get distribution in mean shock path
        Phi = TransitionPhi(rf_ms,r_ms,netwage_ms,pens_ms,xgrid_ms,apgrid_ms,stocks_ms,m_etagrid)

        ! Aggregate

        kp_ms           = sum(apgrid_ms *Phi)/(L_N_ratio*(1.0+n)*(1.0+g))
        agg_bond_demand = sum((apgrid_ms-stocks_ms)*Phi)
        ! Exess demands
        distance(1)     = kp_ms - msvars(1)
        distance(2)     = kp_ms * de_ratio/(1.0 + de_ratio) - agg_bond_demand/(L_N_ratio*(1.0+n)*(1.0+g))

    end function ms_equilibrium


    subroutine get_equilibrium_values(fine,v_fine,apgrid_ms, stocks_ms, xgrid_ms, kappa_ms, value_ms, Phi, errs)
    ! Solve for the MSE one time given the MSE value for k and mu in order to get the (other) equilibrium objects.
    ! The procedure is would be pure pure but for the OMP directives in olg_backwards_solution (but it does read access host variables).
        use params_mod             ,only: L_N_ratio, n, g, stat_dist_z, de_ratio, nx_factor
        use error_class            ,only: tErrors
        use household_solution_mod ,only: olg_backwards_recursion
        use distribution           ,only: TransitionPhi
        use interpolate_xgrid      ,only: InterpolateXgrid
        use income

        implicit none
        type(tPolicies) ,intent(out) :: fine
        type(tErrors)   ,intent(out) :: errs
        real(dp) ,dimension(:,:,:,:,:,:) ,allocatable ,intent(out) :: v_fine
        real(dp) ,dimension(:,:,:)       ,allocatable ,intent(out) :: apgrid_ms, stocks_ms, xgrid_ms, kappa_ms, value_ms, Phi ! mean shock projections and distribution
        type(tPolicies)       :: policies
        real(dp) ,allocatable :: value(:,:,:,:,:,:)
        real(dp)              :: netwage_ms, pens_ms, r_ms, rf_ms
        integer               :: i

        call olg_backwards_recursion(policies,coeffs, grids, value, errs)

        call InterpolateXgrid(nx_factor, policies, value, fine, v_fine)

        allocate(apgrid_ms(size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)), &
                 stocks_ms(size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)), &
                 xgrid_ms (size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)), &
                 kappa_ms (size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)), &
                 value_ms (size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4))  )

        ! Projection of policies / grids on mean shock
        xgrid_ms =0.0; apgrid_ms =0.0; stocks_ms =0.0
        do i=1,size(fine%xgrid,3)
            xgrid_ms  = xgrid_ms  + w(i)* fine%xgrid (:,:,i,:,1,1)
            apgrid_ms = apgrid_ms + w(i)* fine%apgrid(:,:,i,:,1,1)
            stocks_ms = stocks_ms + w(i)* fine%stocks(:,:,i,:,1,1)
            value_ms  = value_ms  + w(i)* v_fine(:,:,i,:,1,1) ! I don't need it right here, but in meanshock_wrapper
        enddo
        where (apgrid_ms .ne. 0.0)
            kappa_ms = stocks_ms/apgrid_ms
        elsewhere
            kappa_ms = 0.0
        end where

        ! Prices in mean shock path
        netwage_ms    = f_netwage (grids%k(1), mean_zeta)
        pens_ms       = f_pensions(grids%k(1), mean_zeta)
        rf_ms         = f_riskfree_rate(grids%k(1),grids%mu(1),stat_dist_z)
        r_ms          = f_stock_return(grids%k(1), mean_zeta, mean_delta, rf_ms)

        ! Get distribution in mean shock path
        Phi = TransitionPhi(rf_ms,r_ms,netwage_ms,pens_ms,xgrid_ms,apgrid_ms,stocks_ms,m_etagrid)

    end subroutine get_equilibrium_values

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
        simvars%stock(1) = sum(stocks_ms*Phi)/ L_N_ratio ! different from K(t+1) since it is in today's per capita terms. Thus no (1+g)(1+n) in denominator.
        simvars%bonds(1) = sum((apgrid_ms-stocks_ms)*Phi)/ L_N_ratio
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
        simvars%bequests(1)   = f_bequests(simvars%rf(1), simvars%r(1), stocks_ms, apgrid_ms, Phi)
        simvars%err_aggr(1)   = f_aggregate_diff(simvars%output(1), simvars%invest(1), simvars%C(1), simvars%bequests(1))
        simvars%B(1)          = fvals(2)
        simvars%err_income(1) = f_income_diff(simvars%K(1), mean_zeta, simvars%r(1), simvars%rf(1), mean_delta)

        do i= 1,nz  ! can't I just save z and then call simulate_economy?
            simvars%z(i+1)     = i
            simvars%K(i+1)     = sum(policies%apgrid(:,:,i,:,1,1)*Phi)/(L_N_ratio*(1.0+n)*(1.0+g))
            simvars%mu(i+1)    = simvars%mu(1)            ! Keep mu constant, could keep rf constant instead
	        simvars%output(i+1)= zeta(i)*simvars%K(i+1)**alpha
	        simvars%stock(i+1) = sum(policies%stocks(:,:,i,:,1,1)*Phi)/ L_N_ratio ! different from K(t+1) since it is in today's per capita terms. Thus no (1+g)(1+n) in denominator.
	        simvars%bonds(i+1) = sum((policies%apgrid(:,:,i,:,1,1)-policies%stocks(:,:,i,:,1,1))*Phi)/ L_N_ratio
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
	        simvars%bequests(i+1)   = f_bequests(simvars%rf(i+1), simvars%r(i+1), policies%stocks(:,:,i,:,1,1), policies%apgrid(:,:,i,:,1,1), Phi)
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
        lifecycles%stock   = sum(sum(stocks_ms * Phi,1),1)
        lifecycles%return  = sum(sum(Phi*sign(1.0,apgrid_ms)*(1.0 + simvars%rf(1) + kappa_ms*simvars%mu(1))/(1.0+g),1),1)
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

end subroutine solve_meanshock
end module mean_shock_mod
