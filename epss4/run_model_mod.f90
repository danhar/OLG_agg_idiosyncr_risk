module run_model_mod
    implicit none
contains

subroutine run_model(projectname, calib_name, welfare, welfare_ins, simvars_o, cal_iter_o, agg_cons_o)
    use classes_mod
    use ifport            ,only: system  ! Intel Fortran portability library
    use omp_lib           ,only: OMP_get_max_threads
    use laws_of_motion    ,only: Initialize
	use markovchain_mod   ,only: MarkovChain
	use distribution      ,only: CheckPhi
	use mean_shock_mod    ,only: solve_meanshock
    use krusell_smith_mod ,only: solve_krusellsmith
    use simvars_class     ,only: read_unformatted, write_unformatted
    use insurance_effect_mod ,only: calc_insurance_effect
	use params_mod        ,only: construct_path, set_apmax, SaveParams, & ! procedures
	                             partial_equilibrium, estimate_from_simvars, mean_return_type, calc_insurance_effects, welfare_decomposition,& ! logicals and characters
	                             dp, nk,nmu, nz, nt, ms_guess, factor_k, factor_mu,cover_k, cover_mu, k_min,k_max,mu_min,mu_max,pi_z, seed, scale_AR

	character(len=*) ,intent(in)  :: projectname, calib_name
	character(len=*) ,intent(in) ,optional :: cal_iter_o
	real(dp)         ,intent(out) :: welfare, welfare_ins(:)! expected ex-ante utility of a newborn
	type(tSimvars)   ,intent(out) ,optional ,allocatable :: simvars_o(:)
	real(dp)         ,intent(out) ,optional :: agg_cons_o
    type(tCoeffs)     :: coeffs, coeffs_old  ! coefficients for laws of motion
    type(tPolicies)   :: policies
    type(tAggGrids)   :: grids, ms_grids
    type(tLifecycle)  :: lifecycles
    type(tSimvars), allocatable:: simvars(:), simvars_old(:)    ! simulation variables, e.g. zt, kt,...
    type(tStats)      :: K, mu, rf, cons
    type(tErrors)     :: err
    real(dp),allocatable :: value(:,:,:,:,:,:), Phi(:,:,:), xgrid_ms(:,:,:) ! value function, distribution, and mean shock xgrid
    real(dp)          :: ms_rf_temp
	integer           :: start_time, it, i, syserr ! 'it' cointains total iterations of Krusell-Smith loop
	logical           :: calibrating
	character(:),allocatable :: dir, output_path

    call system_clock(start_time)
    if (present(simvars_o)) then
        calibrating = .true.
    else
        calibrating =.false.
    endif

    call SaveParams(projectname, calib_name, cal_iter_o)

!-------------------------------------------------------------------------------
! Mean Shock (partial or general) Equilibrium (to generate good initial guesses)
!-------------------------------------------------------------------------------
    print*, ' '
    if (partial_equilibrium) then
        print*,'- run_model: mean shock PARTIAL equilibrium'
        dir    = 'mspe'
        call ms_grids%read_unformatted('ms')
        if (scale_AR == -1.0) then
            print*,'- run_model: setting mu = 0.0, mean return to type '//mean_return_type
            ms_grids%mu =0.0
            ms_grids%k = inverted_mean_return(mean_return_type)
        endif
    else
        print*,'- run_model: mean shock GENERAL equilibrium'
        dir    = 'msge'
        ms_grids = ms_guess
        ! call ms_grids%read_unformatted('ms')
    endif

    if (calibrating) then
        output_path = construct_path(calib_name)
    else
        output_path = construct_path(calib_name,dir)
        syserr = system('mkdir '//output_path//' > /dev/null 2>&1') ! Creates directory for output files, suppresses error if dir exists
    endif

    it = 0  ! no Krusell-Smith iterations in Mean shock (but variable still needed for saving results)
    allocate(simvars(1)) ! only one core used for mean shock
    call solve_meanshock(coeffs, ms_grids, policies, simvars(1), lifecycles, Phi, xgrid_ms, value, err, output_path)
    if (err%not_converged) call err%print2stderr(dir)
    welfare = calc_average_welfare(simvars)

    ! Check distribution and save results
    print*, ' '
    if (.not. calibrating) call CheckPhi(Phi,output_path) ! writes errors to file
    if (.not. partial_equilibrium  .and. .not. err%not_converged) call ms_grids%write_unformatted('ms')
    if (.not. calibrating) call save_and_plot_results(dir, ms_grids, err)

    if (scale_AR == -1.0) then
        if (present(simvars_o)) simvars_o = simvars
        if (.not. calibrating) then
            print*,' **** Completed solution of calibration ', calib_name, ' **** '
            print*, ' '
            print*, ' '
        else
            print*,' *** Model solution successful during calibration of ', calib_name
            call save_and_plot_results(dir, ms_grids, err)
        endif
        if (present(agg_cons_o)) agg_cons_o = simvars(1)%C(1)
        return
    endif
    ms_rf_temp = simvars(1)%rf(1)
    deallocate(simvars)

!stop 'Stop after ms'
!-------------------------------------------------------------------------------
! Non-stationary (partial or general) Equilibrium (i.e. Krusell/Smith algorithm)
!-------------------------------------------------------------------------------
    ! set new apmax using the results of mean shock
    call set_apmax(ms_grids%k(1)*factor_k)

    if(partial_equilibrium) then
        print*,'- run_model: Krusell-Smith PARTIAL equilibrium'
        dir    = 'pe'
        call read_unformatted_ks(grids, coeffs, simvars)
        if (scale_AR == -1.0) then
            ! This is never executed at the moment because of the conditional return in line 85
            print*,'scale_AR = -1.0, i.e. no aggregate risk'
            print*,'setting size(grids%mu)=1 with grids%mu= 0, i.e. zero expected equity premium'
            deallocate(grids%mu); allocate(grids%mu(1))
            grids%mu = 0.0
            print*,'setting simvars%mu= 0'
            do i=1,size(simvars)
                simvars(i)%mu = 0.0
            enddo
        endif
        simvars_old = simvars
        coeffs_old  = coeffs

    else
        print*,'- run_model: Krusell-Smith GENERAL equilibrium'
        dir    = 'ge'
        call grids%allocate(nk,nmu)
        call grids%set_params(k_min,k_max,mu_min,mu_max)
        call grids%construct(ms_grids,factor_k,factor_mu,cover_k, cover_mu, nk,nmu)

        allocate(simvars(OMP_get_max_threads()))
        call simvars%allocate(nt)
        call random_seed(put=seed) ! so that same sequence for different experiments

        if (estimate_from_simvars) then
            print*, '- run_model: using previous simvars to initialize'
            call read_unformatted(simvars_old)
            K%name ='K' ; call K%calc_stats(simvars_old)
            mu%name='mu'; call mu%calc_stats(simvars_old)
            rf%name='rf'; call rf%calc_stats(simvars_old)

            do i=1,size(simvars)
                simvars(i)%z     = MarkovChain(pi_z,nt)
                simvars(i)%K(1)  = K%avg_()
                simvars(i)%mu(1) = mu%avg_()
                simvars(i)%rf(1) = rf%avg_()
            enddo
            call grids%update(K%min_(), K%max_(), mu%min_(), mu%max_())
            ! call grids%update(K%avg_(), mu%avg_(), K%std_(), mu%std_())

            coeffs = Initialize(simvars_old)
        else
            print*, '- run_model: using mean shock values (or hard-coded guesses) to initialize'
            do i=1,size(simvars)
                simvars(i)%z     = MarkovChain(pi_z,nt)
                simvars(i)%K(1)  = ms_grids%k(1)    ! starting value for simulation
                simvars(i)%mu(1) = ms_grids%mu(1)   ! starting value for simulation
                simvars(i)%rf(1) = ms_rf_temp       ! starting value for simulation
            enddo

            coeffs = Initialize(ms_grids)
        endif

    endif

    if (calibrating) then
        output_path = construct_path(calib_name)
    else
        output_path = construct_path(calib_name,dir)
        syserr = system('mkdir '//output_path//' > /dev/null 2>&1') ! Creates directory for output files, suppresses error if dir exists
    endif

    call solve_krusellsmith(grids, projectname, calib_name, output_path, it, coeffs, simvars, Phi, xgrid_ms, policies, value, lifecycles, err, calibrating)
    if (err%not_converged) call err%print2stderr(dir)
    welfare = calc_average_welfare(simvars)

    ! Check distribution and save results
    print*, ' '
    if (.not. calibrating) call CheckPhi(Phi,output_path)
    if (.not. partial_equilibrium .and. .not. err%not_converged) call save_unformatted(grids, coeffs, simvars)
    call save_and_plot_results(dir, grids, err)
    if (present(simvars_o)) simvars_o = simvars

    if (present(agg_cons_o)) then
        cons%name='cons'
        call cons%calc_stats(simvars)
        agg_cons_o = cons%avg_()
    endif

    welfare_ins =1.0 ! not zero to avoide divide by zero
    if (.not. calibrating .and. calc_insurance_effects .and. (index(calib_name,'GE1')>0 .or. .not. welfare_decomposition)) then
        if (partial_equilibrium) then
            simvars = simvars_old
            coeffs  = coeffs_old
        endif
        call calc_insurance_effect(policies, value, grids, simvars, Phi, coeffs, calib_name, projectname, welfare_ins)
    endif

    if (.not. calibrating) then
        print*,' **** Completed solution of calibration ', calib_name, ' **** '
        print*, ' '
        print*, ' '
    else
        print*,' *** Model solution successful during calibration of ', calib_name
    endif

    ! Solution of mean shock OLG model during transition ?
    ! Solution of model with stochastic transition ?

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - subroutine save_and_plot_results(dir, grids, err)
! - pure real(dp) function calc_average(simvars%welfare)
! - real(dp) function average_of_simulations()
! - subroutine save_unformatted(grids, coeffs, simvars)
! - subroutine read_unformatted(grids, coeffs, simvars)
!-------------------------------------------------------------------------------

    subroutine save_and_plot_results(dir, grids, err)
    ! Calc stats and save all results for each run, then plot and save to pdf
        use save_results_mod
	    character(len=*) ,intent(in)   :: dir
	    type(tAggGrids)  ,intent(in)   :: grids
        type(tErrors)    ,intent(in)   :: err
        real(dp) :: secs
        integer  :: end_time, count_rate

        if (.not. calibrating) print*, '- run_model: Saving results and plots to folder ./model_output '
        call system_clock(end_time,count_rate)
        secs= real(end_time-start_time,dp)/real(count_rate,dp) ! excludes time for saving and plotting results

        call save_results(Phi, simvars, coeffs, grids,lifecycles,&
                             policies, secs, it, projectname, calib_name, dir, err, cal_iter_o)

        if (.not. calibrating) call plot_results(output_path, 'plot_all')

    end subroutine save_and_plot_results

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
    real(dp) function inverted_mean_return(mean_return_type)
        use statistics        ,only: tStats
        use simvars_class ,only: read_unformatted
        use params_mod, only: del_mean, de_ratio
        use income, only: f_net_mpk, alpha
        character(len=*) ,intent(in)  :: mean_return_type
        type(tSimvars) ,allocatable :: simvars_ge(:)
        type(tAggGrids)   :: ms_grids_temp
        type(tStats)   :: mpk, r, rf, r_pf_median, r_pf_kappa_med
        real(dp) :: mean_return
        integer :: tc

        call read_unformatted(simvars_ge)

        select case(mean_return_type)
        case('mean_mpk')
            mpk%name='mpk'; call mpk%calc_stats(simvars_ge)
            mean_return = mpk%avg_exerr_()
!           mean_return = real(3.07358605E-02,dp) ! old value
        case('weighted_aggregate_return')
            r%name='r';   call r%calc_stats(simvars_ge)
            rf%name='rf'; call rf%calc_stats(simvars_ge)
!            rf = real(1.84343749E-02,dp)   ! old value
!            r  = real(7.48861757E-02,dp)   ! old value
            mean_return = (r%avg_exerr_()+rf%avg_exerr_()*de_ratio)/(1.0+de_ratio)
        case('median_portf_return')
            r_pf_median%name='r_pf_median';   call r_pf_median%calc_stats(simvars_ge)
            mean_return = r_pf_median%avg_exerr_()
            ! mean_return = real(1.23087113E-01,dp)  ! This is the value for pfrw: portfolio return unweighted
        case('median_portf_share') ! not sure this one is smart or correct
            r_pf_kappa_med%name='r_pf_kappa_med';   call r_pf_kappa_med%calc_stats(simvars_ge)
            r%name='r';   call r%calc_stats(simvars_ge)
            rf%name='rf'; call rf%calc_stats(simvars_ge)
            mean_return = r_pf_kappa_med%avg_exerr_() * r%avg_exerr_()+ (1.0 - r_pf_kappa_med%avg_exerr_())* rf%avg_exerr_()

            ! mean_return = real(1.06801078E-01,dp) ! this is the value for pfr: portfolio return (with pop weights)
        case('Siegel2002') ! empirical estimate from Siegel (2002)
            mean_return = 0.042_dp

        case('meanshock_eq_r') ! rf from previous mean shock equilibrium
            call ms_grids_temp%read_unformatted('ms')
            inverted_mean_return = ms_grids_temp%k(1)
            return
        end select

        inverted_mean_return = (alpha/(mean_return+del_mean))**(1.0/(1.0-alpha))

    end function inverted_mean_return
!-------------------------------------------------------------------------------

    subroutine save_unformatted(grids, coeffs, simvars)
        type(tAggGrids) ,intent(in) :: grids
        type(tCoeffs)   ,intent(in) :: coeffs
        type(tSimvars)  ,intent(in) :: simvars(:)

        call grids%write_unformatted('ge')
        call coeffs%write_unformatted

        call write_unformatted(simvars)

    end subroutine save_unformatted
!-------------------------------------------------------------------------------

    subroutine read_unformatted_ks(grids, coeffs, simvars)
        type(tAggGrids) ,intent(out) :: grids
        type(tCoeffs)   ,intent(out) :: coeffs
        type(tSimvars)  ,allocatable ,intent(out) :: simvars(:)

        call grids%read_unformatted('ge')
        call coeffs%read_unformatted

        call read_unformatted(simvars)
    end subroutine read_unformatted_ks
!-------------------------------------------------------------------------------

end subroutine run_model

end module run_model_mod
