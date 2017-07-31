!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

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
    use classes_mod ,only: tPolicies, tAggGrids, tErrors, tSimvars, tLifecycle, tCoeffs
    use laws_of_motion  ,only: Initialize
    use params_mod      ,only: alpha,etagrid,stat_dist_z, partial_equilibrium, surv_rates, scale_AR, tol_coeffs
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
    type(tPolicies)                :: policies_ms
	real(dp) ,dimension(:), allocatable  :: xvars, fvals, m_etagrid	! input/output of ms_equilib, passed to sub_broyden. Mean etagrid.
	real(dp)                       :: mean_zeta, mean_delta, bequests_ms ! mean shocks
	real(dp) ,dimension(:,:,:), allocatable :: value_ms ! mean shock
	logical ,parameter :: normalize_xvars_to_unity = .true.

    coeffs          = Initialize()
	mean_zeta	    = dot_product(stat_dist_z, zeta)
	mean_delta		= dot_product(stat_dist_z, delta)
	m_etagrid	    = matmul(etagrid,stat_dist_z)

    ! Initial guesses
    if (scale_AR == -1.0) then
        allocate(xvars(1), fvals(1))
        xvars(1)   = grids%k (1)
        if (normalize_xvars_to_unity) xvars(1) = 1.0
    else
        allocate(xvars(2), fvals(2))
        xvars(1)   = grids%k (1)
        xvars(2)   = grids%mu(1)
        if (normalize_xvars_to_unity) xvars = 1.0
    endif

    if(surv_rates) then
        bequests_ms = 0.03_dp   ! This guess was the mean shock value in previous runs
    else
        bequests_ms = 0.0
    endif

	if (partial_equilibrium .or. (grids%fixed .and. .not. scale_AR == -1.0)) then
	! if grids%fixed then we do not need a guess for the aggregate grid in GE.
	! However, grids%fixed is source-set. For the no AR and case, we want to calc the GE.
	    fvals = ms_equilibrium(xvars)
	else ! Find capital and mu for the mean shock general equilibrium
		call s_broyden(ms_equilibrium,xvars,fvals,err%not_converged, get_fd_jac_o=.true.,tolf_o=tol_coeffs, maxlnsrch_o=5)  ! ,maxstp_o=.8_dp,tolf_o=1e-6_dp
		if (err%not_converged) call err%write2file(fvals, output_path)
    endif

    if (normalize_xvars_to_unity) then
        grids%k(1) = xvars(1)*grids%k(1)
    else
	    grids%k(1) = xvars(1)
    endif

	if (size(xvars)==1) then
	    grids%mu(1) = 0.0
	else
	    if (normalize_xvars_to_unity) then
	        grids%mu(1) = xvars(2)*grids%mu(1)
	    else
	        grids%mu(1) = xvars(2)
        endif
    endif

    call get_equilibrium_values(policies,value,policies_ms, value_ms, Phi, err)
    xgrid_ms = policies_ms%xgrid(:,:,1,:,1,1)
    call err%print2stderr

    call simulate_ms(simvars)

    call ms_lc_profiles(lifecycles)

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - (pure) function ms_equilibrium(msvars) result(distance)
! - (pure) subroutine get_equilibrium_values(fine,v_fine,apgrid_ms, stocks_ms, xgrid_ms, kappa_ms, value_ms, Phi, errs)
! - pure subroutine simulate_ms(simvars)
! - pure subroutine ms_lc_profiles(lifecycles)
!-------------------------------------------------------------------------------


    function ms_equilibrium(msvars) result(distance)
    ! Solve for the MSE, i.e. where k'=k and mu'=mu given that all realizations are at the mean_z.
    ! The procedure is would be pure pure but for the OMP directives in olg_backwards_solution (but it does read access host variables).
    ! Update: if I update the bequests_ms this is also non-pure. A more correct (but superfluous) way would be to find the fixed-point in (bequests,Phi).
        use params_mod             ,only: L_N_ratio, n, g, stat_dist_z, de_ratio, nx_factor
        use error_class            ,only: tErrors
        use household_solution_mod ,only: olg_backwards_recursion
        use distribution           ,only: TransitionPhi
        use interpolate_xgrid      ,only: InterpolateXgrid
        use income

        implicit none
        real(dp) ,dimension(:),intent(in) :: msvars             ! k_ms, mu_ms
        real(dp) ,dimension(size(msvars)) :: distance           ! excess demands
        type(tPolicies) :: policies, fine, policies_ms
        type(tAggGrids) :: grid
        type(tErrors)   :: errs
        real(dp) ,dimension(:,:,:,:,:,:) ,allocatable :: value, v_fine
        real(dp) ,dimension(:,:,:) ,allocatable :: Phi ! mean shock projections and distribution
        real(dp)                          :: netwage_ms, pens_ms, r_ms, rf_ms, kp_ms, agg_bond_demand
        integer                           :: i, nx_factor_ms

        call grid%allocate(size(grids%k), size(grids%mu))

        grid%k  = msvars(1)
        if (normalize_xvars_to_unity) grid%k  = msvars(1) *grids%k(1)

        if (size(msvars) == 1) then
            grid%mu = 0.0
        else
            grid%mu = msvars(2)
            if (normalize_xvars_to_unity) grid%mu = msvars(2) *grids%mu(1)
        endif

        call olg_backwards_recursion(policies,coeffs, grid, value, errs)

        nx_factor_ms = nx_factor * 20 ! use higher nx_factor than in KS, though it didn't help much
        call InterpolateXgrid(nx_factor_ms, policies, value, fine, v_fine)

        policies_ms = fine%mean(3,stat_dist_z)

        ! Prices in mean shock path
        netwage_ms    = f_netwage (grid%k(1), mean_zeta)
        pens_ms       = f_pensions(grid%k(1), mean_zeta)
        rf_ms         = f_riskfree_rate(grid%k(1),grid%mu(1),stat_dist_z)
        r_ms          = f_stock_return(grid%k(1), mean_zeta, mean_delta, rf_ms)

        ! Get distribution in mean shock path
        Phi = TransitionPhi(rf_ms,r_ms,netwage_ms,pens_ms,bequests_ms,policies_ms%xgrid(:,:,1,:,1,1),policies_ms%apgrid(:,:,1,:,1,1),policies_ms%stocks(:,:,1,:,1,1),m_etagrid)

        ! Update the guess for bequests (instead of finding a fixed-point in (Phi,bequests)
        ! bequests_ms   = f_bequests(rf_ms, r_ms, stocks_ms, apgrid_ms, Phi) ! This update seemed to create problems for the rootfinder.

        ! Aggregate

        kp_ms           = sum(policies_ms%apgrid(:,:,1,:,1,1)*Phi)/(L_N_ratio*(1.0+n)*(1.0+g))
        agg_bond_demand = sum((policies_ms%apgrid(:,:,1,:,1,1)-policies_ms%stocks(:,:,1,:,1,1))*Phi)
        ! Exess demands
        distance(1)     = kp_ms - grid%k(1)
        if (size(distance) > 1) distance(2) = kp_ms * de_ratio/(1.0 + de_ratio) - agg_bond_demand/(L_N_ratio*(1.0+n)*(1.0+g))

    end function ms_equilibrium


    subroutine get_equilibrium_values(fine,v_fine,policies_ms, value_ms, Phi, errs)
    ! Solve for the MSE one time given the MSE value for k and mu in order to get the (other) equilibrium objects.
    ! The procedure is would be pure pure but for the OMP directives in olg_backwards_solution (but it does read access host variables).
    ! Update: the update of bequests_ms is also non-pure. A more correct (but superfluous) way would be to find the fixed-point in (bequests,Phi).
        use params_mod             ,only: L_N_ratio, n, g, stat_dist_z, de_ratio, nx_factor
        use error_class            ,only: tErrors
        use household_solution_mod ,only: olg_backwards_recursion
        use distribution           ,only: TransitionPhi
        use interpolate_xgrid      ,only: InterpolateXgrid
        use income

        implicit none
        type(tPolicies) ,intent(out) :: fine, policies_ms
        type(tErrors)   ,intent(out) :: errs
        real(dp) ,dimension(:,:,:,:,:,:) ,allocatable ,intent(out) :: v_fine
        real(dp) ,dimension(:,:,:)       ,allocatable ,intent(out) :: value_ms, Phi ! mean shock projections and distribution
        type(tPolicies)       :: policies, pol_ms
        real(dp) ,allocatable :: value(:,:,:,:,:,:)
        real(dp)              :: netwage_ms, pens_ms, r_ms, rf_ms
        integer               :: i, nx_factor_ms

        call olg_backwards_recursion(policies,coeffs, grids, value, errs)

        nx_factor_ms = nx_factor * 20 ! use higher nx_factor than in KS, though it didn't help much
        call InterpolateXgrid(nx_factor_ms, policies, value, fine, v_fine)

        policies_ms = fine%mean(3,stat_dist_z)
        call policies_ms%calc_kappa()
        allocate(value_ms(size(fine%xgrid,1),size(fine%xgrid,2),size(fine%xgrid,4)))
        value_ms = 0.0
        do i=1,size(fine%xgrid,3)
            value_ms  = value_ms  + stat_dist_z(i)* v_fine(:,:,i,:,1,1)
        enddo

        ! Prices in mean shock path
        netwage_ms    = f_netwage (grids%k(1), mean_zeta)
        pens_ms       = f_pensions(grids%k(1), mean_zeta)
        rf_ms         = f_riskfree_rate(grids%k(1),grids%mu(1),stat_dist_z)
        r_ms          = f_stock_return(grids%k(1), mean_zeta, mean_delta, rf_ms)

        ! Get distribution in mean shock path
        Phi = TransitionPhi(rf_ms,r_ms,netwage_ms,pens_ms,bequests_ms,policies_ms%xgrid(:,:,1,:,1,1),policies_ms%apgrid(:,:,1,:,1,1),policies_ms%stocks(:,:,1,:,1,1),m_etagrid)

    end subroutine get_equilibrium_values

    pure subroutine simulate_ms(simvars)
    ! Variable simvars is used to get a very rough approximation of standard deviations.
    ! At index 1, store mean shock values of aggregates.
    ! At index i+1, store aggregates for z=i, calculated by assuming:
    ! - agents live in a world where z(i) always realizes, use LOMs of mean shock (i.e. kp =k, mup = mu)
    ! - Phi remains constant at Phi_ms

        use params_mod   ,only: L_N_ratio, n, g, pi_z, de_ratio, calc_euler_errors
        use fun_aggregate_diff
        use distribution ,only: CheckPhi
        use simulation_mod, only: calc_inequality_measures, f_euler_errors
        use sorting_partial_mod     ! function valnth

        type(tSimvars) ,intent(out) :: simvars
        real(dp) ,allocatable       :: r_pf(:,:,:)
        real(dp) ,dimension(2)      :: eul_err_temp
        integer                     :: i, nz

        associate(xgrid_ms  => policies_ms%xgrid (:,:,1,:,1,1), &
                  apgrid_ms => policies_ms%apgrid(:,:,1,:,1,1), &
                  stocks_ms => policies_ms%stocks(:,:,1,:,1,1), &
                  kappa_ms  => policies_ms%kappa (:,:,1,:,1,1)  )

        nz = size(stat_dist_z)
        call simvars%allocate(nz+1)   ! allocate all simulation variables
        simvars%z(1)    = 0     ! to indicate that these are mean-shock-results
        simvars%K(1)    = grids%k(1)
        simvars%mu(1)   = grids%mu(1)

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
        if (size(fvals)>1) then
            simvars%B(1)      = fvals(2)
        else
            simvars%B(1)      = simvars%K(1) * de_ratio/(1.0 + de_ratio) - simvars%bonds(1) ! differs from fvals(2)
        endif
        simvars%err_income(1) = f_income_diff(simvars%K(1), mean_zeta, simvars%r(1), simvars%rf(1), mean_delta)
        ! the euler errors work only for the no AR economy, because then the zt does not matter. for MS we would need a zt for mean shock
        if (.not. calc_euler_errors) then
            simvars%eul_err_max(1)=0.0
            simvars%eul_err_avg(1)=0.0
        else
            eul_err_temp = f_euler_errors(1, simvars%rf(1), simvars%mu(1),simvars%K(1),coeffs, grids, policies_ms, value, xgrid_ms, apgrid_ms, kappa_ms, Phi)
            simvars%eul_err_max(1)=eul_err_temp(1)
            simvars%eul_err_avg(1)=eul_err_temp(2)
        endif

        call calc_inequality_measures(simvars, xgrid_ms, apgrid_ms, stocks_ms, Phi, m_etagrid, simvars%pens(1), simvars%wage(1), 1)

        do i= 1,nz  ! can't I just save z and then call simulate_economy?
            simvars%z(i+1)     = i
            simvars%K(i+1)     = sum(policies%apgrid(:,:,i,:,1,1)*Phi)/(L_N_ratio*(1.0+n)*(1.0+g))
            if (simvars%K(i+1) < 0.1_dp) simvars%K(i+1) = 0.1_dp  ! This can happen particularly in the partial equilibrium experiments. It corresponds to the aggregate grid bounds check in the true simulations.
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
            simvars%eul_err_max(i+1)= 0.0
            simvars%eul_err_avg(i+1)= 0.0
            call calc_inequality_measures(simvars,policies%xgrid(:,:,i,:,1,1), policies%apgrid(:,:,i,:,1,1), policies%stocks(:,:,i,:,1,1), Phi, etagrid(:,i), simvars%pens(i+1), simvars%wage(i+1), i+1)
        enddo

        simvars%K (nz+2) = simvars%K (1) ! K and rf have one more index in the real simulations.
        simvars%rf(nz+2) = simvars%rf(1)

        end associate
    end subroutine simulate_ms

    pure subroutine ms_lc_profiles(lifecycles)
    ! Calculates the policies in the mean shock and then the corresponding life-cycle profiles
        use params_mod   ,only        :  g, pop_frac
        type(tLifecycle) ,intent(out) :: lifecycles
        integer                       :: jc, nx, n_eta, nj
        real(dp) ,dimension(:,:,:) ,allocatable :: cons_ms, weight

        associate(xgrid_ms  => policies_ms%xgrid (:,:,1,:,1,1), &
                  apgrid_ms => policies_ms%apgrid(:,:,1,:,1,1), &
                  stocks_ms => policies_ms%stocks(:,:,1,:,1,1), &
                  kappa_ms  => policies_ms%kappa (:,:,1,:,1,1)  )

        call lifecycles%allocate(size(apgrid_ms,3))
        ! Intel Fortran Compiler XE 13.0 Update 1 (and previous) has a bug on realloc on assignment. If that is corrected, I think I can remove the following allocation block.
        nx = size(xgrid_ms,1); n_eta = size(xgrid_ms,2); nj = size(xgrid_ms,3);
        allocate(cons_ms(nx,n_eta,nj), weight(nx,n_eta,nj))


        do jc=1,nj
            weight(:,:,jc) = Phi(:,:,jc)/pop_frac(jc)
        enddo

        cons_ms = xgrid_ms-apgrid_ms
        lifecycles%ap      = sum(sum(apgrid_ms * weight,1),1)
        lifecycles%cons    = sum(sum(cons_ms * weight,1),1)
        lifecycles%log_cons= sum(sum(log(cons_ms) * weight,1),1)
        lifecycles%stock   = sum(sum(stocks_ms * weight,1),1)
        lifecycles%return  = sum(sum(weight*sign(1.0,apgrid_ms)*(1.0 + simvars%rf(1) + kappa_ms*simvars%mu(1))/(1.0+g),1),1)
        do jc=1,size(apgrid_ms,3)
            lifecycles%cons_var(jc)     = sum(((cons_ms(:,:,jc) - lifecycles%cons(jc)))**2 * weight(:,:,jc))
            lifecycles%var_log_cons(jc) = sum(((log(cons_ms(:,:,jc)) - lifecycles%log_cons(jc)))**2 * weight(:,:,jc))
            lifecycles%return_var(jc)   = sum((apgrid_ms(:,:,jc)*(1.0 + simvars%rf(1) + kappa_ms(:,:,jc)*simvars%mu(1))/(1.0+g) - lifecycles%return(jc))**2 * weight(:,:,jc))
        enddo
        lifecycles%exp_value = value_ms * Phi ! recall that this is not a lifecycle profile.
        lifecycles%xgrid     = xgrid_ms ! recall that this is not a lifecycle profile.

        where (lifecycles%ap .ne. 0.0)
            lifecycles%kappa = lifecycles%stock/lifecycles%ap
        elsewhere
            lifecycles%kappa = 0.0
        end where

        end associate
    end subroutine ms_lc_profiles

end subroutine solve_meanshock
end module mean_shock_mod
