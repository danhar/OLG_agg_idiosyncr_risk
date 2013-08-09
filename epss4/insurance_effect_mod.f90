module insurance_effect_mod
! In this module, the pure insurance effect is computed, excluding the mean effect.
! This is achieved by averaging the existing policy function over the risk that is thrown out.
! In particular, agents do not re-optimize in the world without the risk.
! However, the resulting numbers are hard to understand. I do not think there is a bug, I rather think that this experiment is conceptually flawed.
! So typically, this module will not be entered, because in params_mod, the default value ist calc_insurance_effects = .false.
! This module was added very late and implements experiments that are hard to insert into the existing code base.
! As a consequence, a lot of code is duplicated here.
    implicit none
contains

subroutine calc_insurance_effect(policies, value, agg_grid, simvars_in, Phi_in, coeffs, calib_name, projectname, welfare)
    use kinds            ,only: dp
    use classes_mod      ,only: tSimvars, tLifecycle, tPolicies, tAggGrids, tStats, tCoeffs
    use laws_of_motion   ,only: Initialize
    use lifecycles_class ,only: average
    use simvars_class    ,only: print_error
    use fun_locate
    use params_mod       ,only: exogenous_xgrid, nx_factor, stat_dist_z, stat_dist_eta, ccv, scale_IR, scale_AR, scale_IR_orig, etagrid, &
                                params_set, params_set_thisrun, CheckParams
    use simulation_mod   ,only: simulate
    use interpolate_xgrid

    type(tPolicies)  ,intent(in)  :: policies
    real(dp)         ,intent(in)  :: value(:,:,:,:,:,:)
    type(tAggGrids)  ,intent(in)  :: agg_grid
    type(tSimvars)   ,intent(in)  :: simvars_in(:)    ! (zt, kt, mut, bt,...), first element contains starting values
    real(dp)         ,intent(in)  :: Phi_in(:,:,:) ! distribution. Returns: average Phi if (exogenous_xgrid), else Phi in nt
    type(tCoeffs)    ,intent(in)  :: coeffs
    character(len=*) ,intent(in)  :: calib_name, projectname
    real(dp)         ,intent(out) :: welfare(5)
    type(tLifecycle)              :: lifecycles, lifecycles_array(size(simvars_in))         ! lifecycle profiles
    type(tPolicies)               :: pol_fine, pol_minus_risk    ! if (exogenous_xgrid) then this will hold interpolated policies
    type(tStats)                  :: K, mu, netwage, pens, bequests, r, rf
    type(tSimvars) ,allocatable   :: simvars(:)
    type(tAggGrids)               :: agg_grid_noAR
    type(tCoeffs)                 :: coeffs_noAR
    real(dp)       ,allocatable   :: Phi(:,:,:), xgrid_mean_old(:,:,:), xgrid_mean_new(:,:,:), Phi_spread(:,:,:,:)
    real(dp), allocatable ,dimension(:,:,:,:,:,:)   :: val_newx, val_minus_risk
    character(len=7)          :: runchar
    logical                   :: ccv_thisrun
    real(dp)                  :: scale_IR_thisrun, scale_AR_thisrun, w
    integer                   :: i

    ccv_thisrun = ccv
    scale_IR_thisrun = scale_IR
    scale_AR_thisrun = scale_AR

    K%name ='K'; call K%calc_stats(simvars_in)
    mu%name='mu'; call mu%calc_stats(simvars_in)
    netwage%name='netwage'; call netwage%calc_stats(simvars_in)
    pens%name='pension'; call pens%calc_stats(simvars_in)
    bequests%name='bequests'; call bequests%calc_stats(simvars_in)
    r%name='r'; call r%calc_stats(simvars_in)
    rf%name='rf'; call rf%calc_stats(simvars_in)

    call InterpolateXgrid(nx_factor, policies, value, pol_fine, val_newx)
    xgrid_mean_old = sum(pol_fine%xgrid(:,:,1,:,:,1),4)/size(pol_fine%xgrid,5) ! only an approximation of the grid over which Phi is defined

    print *,'- insurance_effect: simulating all risks'
    write(runchar,'(a)') 'all'
    pol_minus_risk = policies
    call sim_pe(welfare(1))

    print *,'- insurance_effect: simulating no CCV'
    call params_set('ccv', .false.)
    write(runchar,'(a)') 'noCCV'
    pol_minus_risk = policies
    ! not possible to interpolate policy functions, since ccv doesn't have its own dimension.
    ! So here only remove ccv from income in simulation, use the re-optimized policies and adjust them with aggregate consumption for behavioral response.
    call sim_pe(welfare(2))

    print *,'- insurance_effect: simulating no IR'
    call params_set('scale_IR', -1.0_dp)
    pol_minus_risk = policies%mean(2,stat_dist_eta)
    write(runchar,'(a)') 'noIR'
    call sim_pe(welfare(3))

    call agg_grid_noAR%allocate(1,1)
    agg_grid_noAR%k  =  K%avg_exerr_()
    agg_grid_noAR%mu = mu%avg_exerr_()
    coeffs_noAR = Initialize() ! mean shock initialization

    print *,'- insurance_effect: simulating no AR'
    call params_set('scale_AR', -1.0_dp)
    pol_minus_risk = policies%mean(3,stat_dist_z)
    pol_minus_risk = pol_minus_risk%interpolate(5,agg_grid%k , K%avg_exerr_())
    pol_minus_risk = pol_minus_risk%interpolate(6,agg_grid%mu,mu%avg_exerr_())
    if (scale_IR_orig .ne. -1.0) call params_set('scale_IR', 0.0_dp)
    write(runchar,'(a)') 'noAR'
    call sim_pe(welfare(4))

    print *,'- insurance_effect: simulating no risk'
    call params_set('scale_AR', -1.0_dp)
    call params_set('scale_IR', -1.0_dp)
    pol_minus_risk = pol_minus_risk%mean(2,stat_dist_eta)
    write(runchar,'(a)') 'norisk'
    call sim_pe(welfare(5))

    call params_set('ccv', ccv_thisrun)
    call params_set('scale_IR', scale_IR_thisrun)
    call params_set('scale_AR', scale_AR_thisrun)
    call params_set_thisrun

contains
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine sim_pe(welf)
        use valuefunction_mod ,only: value_f => value
        use distribution           ,only: TransitionPhi
        real(dp) ,intent(out) ::welf
        Phi     = Phi_in    ! initial values, makes clear that Phi is not output
        simvars = simvars_in

        call params_set_thisrun
        call CheckParams

        if (runchar == 'noCCV' .or. runchar=='noIR' .or. runchar=='all') then
            val_minus_risk = value_f(pol_minus_risk,coeffs,agg_grid)

            if (exogenous_xgrid) then
                ! This is the standard case which should always be used, because we make the xgrid much finer
                call InterpolateXgrid(nx_factor, pol_minus_risk, val_minus_risk, pol_fine, val_newx)
                ! We also want the initial Phi that we take from previous simulations to be defined over the new xgrid
                xgrid_mean_new = sum(pol_fine%xgrid(:,:,1,:,:,1),4)/size(pol_fine%xgrid,5) ! only an approximation of the grid over which Phi is defined
                call InterpolateXgrid(Phi, xgrid_mean_old, xgrid_mean_new)
                Phi_spread = spread(Phi,4,size(simvars))
                !$OMP  PARALLEL DO IF(size(simvars)>1)
                do i=1,size(simvars)
                    call simulate(pol_fine, val_newx, agg_grid, simvars(i), Phi_spread(:,:,:,i), lifecycles_array(i))
                enddo
                !$OMP END PARALLEL DO
            else
                ! In this case, we simultaneously solve for the rf(t+1) (actually the mu(t+1)) and the corresponding xgrid(since it depends on mu).
                ! That is very costly, so we do not refine xgrid. While more correct theoretically, the coarse xgrid makes the solution less precise.
                Phi_spread = spread(Phi,4,size(simvars))
                !$OMP  PARALLEL DO IF(size(simvars)>1)
                do i=1,size(simvars)
                    call simulate(pol_minus_risk, val_minus_risk, agg_grid, simvars(i), Phi_spread(:,:,:,i), lifecycles_array(i))
                enddo
                !$OMP END PARALLEL DO
            endif
            lifecycles=average(lifecycles_array)
            Phi = sum(Phi_spread,4)/size(Phi_spread,4)

            welf = calc_average_welfare(simvars)
        else
            val_minus_risk = value_f(pol_minus_risk,coeffs_noAR,agg_grid_noAR)
            call InterpolateXgrid(nx_factor, pol_minus_risk, val_minus_risk, pol_fine, val_newx)

            ! Alternatively, one could take Phi from the noCCV simulations, which is an average, but I think the following is more correct
            Phi = TransitionPhi(rf%avg_exerr_(),r%avg_exerr_(),netwage%avg_exerr_(),pens%avg_exerr_(),bequests%avg_exerr_(),pol_fine%xgrid(:,:,1,:,1,1),pol_fine%apgrid(:,:,1,:,1,1),pol_fine%stocks(:,:,1,:,1,1),etagrid(:,1))

            if (allocated(simvars)) deallocate(simvars)
            allocate(simvars(1))
            call simulate_ms(simvars(1))
            call ms_lc_profiles(lifecycles)


            welf = sum(val_newx(:,:,1,1,1,1)*Phi(:,:,1))

        endif

        call save_and_plot_results('insurance/'//trim(adjustl(runchar)),agg_grid)
    end subroutine sim_pe

    !-------------------------------------------------------------------------------
    pure real(dp) function calc_average_welfare(simvars)
        use statistics        ,only: tStats
        type(tSimvars), intent(in) :: simvars(:)
        type(tStats)      :: welfare_stats
        welfare_stats%name = 'welfare'
        call welfare_stats%calc_stats(simvars)
        calc_average_welfare = welfare_stats%avg_exerr_() ! without the periods where simvars hit gridbounds
    end function calc_average_welfare
    !-------------------------------------------------------------------------------

    pure subroutine simulate_ms(simvars)
    ! Variable simvars is used to get a very rough approximation of standard deviations.
    ! At index 1, store mean shock values of aggregates.
    ! At index i+1, store aggregates for z=i, calculated by assuming:
    ! - agents live in a world where z(i) always realizes, use LOMs of mean shock (i.e. kp =k, mup = mu)
    ! - Phi remains constant at Phi_ms

        use params_mod   ,only: L_N_ratio, n, g, pi_z, de_ratio, alpha, zeta, delta
        use income
        use fun_aggregate_diff
        use distribution ,only: CheckPhi
        use partial_sorting     ! function valnth

        type(tSimvars) ,intent(out) :: simvars
        real(dp) ,dimension(:,:,:) ,allocatable       :: r_pf, stocks_ms, apgrid_ms, xgrid_ms, kappa_ms, value_ms
        real(dp)                    :: mean_zeta, mean_delta
        integer                     :: i, nz

        nz = size(pol_fine%apgrid,3)
        mean_zeta       = dot_product(stat_dist_z, zeta)
        mean_delta      = dot_product(stat_dist_z, delta)


        stocks_ms = pol_fine%stocks(:,:,1,:,1,1)
        apgrid_ms = pol_fine%apgrid(:,:,1,:,1,1)
        xgrid_ms  = pol_fine%xgrid(:,:,1,:,1,1)
        where (apgrid_ms .ne. 0.0)
            kappa_ms = stocks_ms/apgrid_ms
        elsewhere
            kappa_ms = 0.0
        end where
        value_ms = val_newx(:,:,1,:,1,1)

        call simvars%allocate(nz+1)   ! allocate all simulation variables
        simvars%z(1)    = 0     ! to indicate that these are mean-shock-results
        simvars%K(1)    = K%avg_exerr_()
        simvars%mu(1)   = mu%avg_exerr_()

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
        simvars%B(1)      = simvars%K(1) * de_ratio/(1.0 + de_ratio) - simvars%bonds(1) ! differs from fvals(2)
        simvars%err_income(1) = f_income_diff(simvars%K(1), mean_zeta, simvars%r(1), simvars%rf(1), mean_delta)

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
        enddo

        simvars%K (nz+2) = simvars%K (1) ! K and rf have one more index in the real simulations.
        simvars%rf(nz+2) = simvars%rf(1)
    end subroutine simulate_ms
!-------------------------------------------------------------------------------

    pure subroutine ms_lc_profiles(lifecycles)
    ! Calculates the policies in the mean shock and then the corresponding life-cycle profiles
        use params_mod   ,only        :  g
        type(tLifecycle) ,intent(out) :: lifecycles
        real(dp) ,dimension(:,:,:) ,allocatable       :: stocks_ms, apgrid_ms, xgrid_ms, kappa_ms
        integer                       :: jc

        stocks_ms = pol_fine%stocks(:,:,1,:,1,1)
        apgrid_ms = pol_fine%apgrid(:,:,1,:,1,1)
        xgrid_ms  = pol_fine%xgrid(:,:,1,:,1,1)
        where (apgrid_ms .ne. 0.0)
            kappa_ms = stocks_ms/apgrid_ms
        elsewhere
            kappa_ms = 0.0
        end where

        call lifecycles%allocate(size(apgrid_ms,3))

        lifecycles%ap      = sum(sum(apgrid_ms * Phi,1),1)
        lifecycles%cons    = sum(sum((xgrid_ms-apgrid_ms) * Phi,1),1)
        lifecycles%stock   = sum(sum(stocks_ms * Phi,1),1)
        lifecycles%return  = sum(sum(Phi*sign(1.0,apgrid_ms)*(1.0 + simvars(1)%rf(1) + kappa_ms*simvars(1)%mu(1))/(1.0+g),1),1)
        do jc=1,size(apgrid_ms,3)
            lifecycles%cons_var(jc)  = sum((((xgrid_ms(:,:,jc)-apgrid_ms(:,:,jc)) - lifecycles%cons(jc)))**2 * Phi(:,:,jc))
            lifecycles%return_var(jc)= sum((apgrid_ms(:,:,jc)*(1.0 + simvars(1)%rf(1) + kappa_ms(:,:,jc)*simvars(1)%mu(1))/(1.0+g) - lifecycles%return(jc))**2 * Phi(:,:,jc))
        enddo

        where (lifecycles%ap .ne. 0.0)
            lifecycles%kappa = lifecycles%stock/lifecycles%ap
        elsewhere
            lifecycles%kappa = 0.0
        end where
    end subroutine ms_lc_profiles
!-------------------------------------------------------------------------------

    subroutine save_and_plot_results(dir, grids)
    ! Calc stats and save all results for each run, then plot and save to pdf
        use ifport             ,only: system  ! Intel Fortran portability library
        use params_mod         ,only: construct_path
        use coefficients_class ,only: tCoeffs
        use error_class        ,only: tErrors
        use save_results_mod
        character(len=*) ,intent(in)   :: dir
        type(tAggGrids)  ,intent(in)   :: grids
        character(:),allocatable  :: output_path
        type(tErrors)    :: err_temp
        real(dp)      :: secs_temp
        integer       :: it_temp, syserr

        print*, '- insurance_effects_mod: Saving results and plots to folder ./model_output '
        output_path=construct_path(calib_name,dir)
        syserr = system('mkdir -p '//output_path//' > /dev/null 2>&1') ! Creates directory for output files, suppresses error if dir exists

        secs_temp = 0.0
        it_temp   = 0
        call    err_temp%allocate(size(pol_minus_risk%xgrid,5), size(pol_minus_risk%xgrid,6))

        call save_results(Phi, simvars, coeffs, grids,lifecycles,&
                             pol_minus_risk, secs_temp, it_temp, projectname, calib_name, dir, err_temp)

        call plot_results(output_path, 'plot_all')

    end subroutine save_and_plot_results


end subroutine calc_insurance_effect

end module insurance_effect_mod
