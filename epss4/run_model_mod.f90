module run_model_mod
    implicit none
contains

subroutine run_model(projectname, calib_name, welfare, simvars_o)
    use classes_mod
    use ifport            ,only: system  ! Intel Fortran portability library
    use omp_lib           ,only: OMP_get_max_threads
    use laws_of_motion    ,only: Initialize
	use markovchain_mod   ,only: MarkovChain
	use distribution      ,only: CheckPhi
	use mean_shock_mod    ,only: solve_meanshock
    use krusell_smith_mod ,only: solve_krusellsmith
    use simvars_class     ,only: read_unformatted, write_unformatted
	use params_mod        ,only: construct_path, set_apmax, & ! procedures
	                             partial_equilibrium, estimate_from_simvars, & ! logicals
	                             dp, nk,nmu, nz, n_coeffs, nt, ms_guess, factor_k, factor_mu,cover_k, cover_mu, pi_z, seed, scale_AR

	character(len=*) ,intent(in)  :: projectname, calib_name
	real(dp)         ,intent(out) :: welfare ! expected ex-ante utility of a newborn
	type(tSimvars)   ,intent(out) ,optional, allocatable :: simvars_o(:)
    type(tCoeffs)     :: coeffs  ! coefficients for laws of motion
    type(tPolicies)   :: policies
    type(tAggGrids)   :: grids, ms_grids
    type(tLifecycle)  :: lifecycles
    type(tSimvars), allocatable:: simvars(:)    ! simulation variables, e.g. zt, kt,...
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
!-------------------------------------------------------------------------------
! Mean Shock (partial or general) Equilibrium (to generate good initial guesses)
!-------------------------------------------------------------------------------
    print*, ' '
    if (partial_equilibrium) then
        print*,'- run_model: mean shock PARTIAL equilibrium'
        dir    = 'mspe'
        call ms_grids%read_unformatted('ms')
        if (scale_AR == -1.0) then
            print*,'- run_model: setting ms_grids%mu = 0.0, ms_grids%k = average over simulations'
            ms_grids%mu =0.0
            ms_grids%k = inverted_average_mpk()
        endif
    else
        print*,'- run_model: mean shock GENERAL equilibrium'
        dir    = 'msge'
        ms_grids = ms_guess
    endif
    output_path = construct_path(dir,calib_name)
    if (.not. calibrating) syserr = system('mkdir '//output_path) ! Creates directory for output files
    it = 0  ! no Krusell-Smith iterations in Mean shock (but variable still needed for saving results)

    allocate(simvars(1)) ! only one core used for mean shock
    call solve_meanshock(coeffs, ms_grids, policies, simvars(1), lifecycles, Phi, xgrid_ms, value, err, output_path)
    if (err%not_converged) call err%print2stderr(dir)
    welfare = calc_average_welfare(simvars)

    ! Check distribution and save results
    if (.not. calibrating) call CheckPhi(Phi,output_path) ! writes errors to file
    if (.not. partial_equilibrium) then
        syserr = system('cp model_input/last_results/*.unformatted model_input/last_results/previous/')
        call ms_grids%write_unformatted('ms')
    endif
    if (.not. calibrating) call save_and_plot_results(dir, ms_grids, err)

    ms_rf_temp = simvars(1)%rf(1)
    deallocate(simvars)
    print*, ' '

    if (scale_AR == -1.0) then
        print*,' **** Completed solution of calibration ', calib_name, ' **** '
        print*, ' '
        print*, ' '
        return
    endif
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
            ! This is never executed at the moment because of the conditional return in line 70
            print*,'scale_AR = -1.0, i.e. no aggregate risk'
            print*,'setting size(grids%mu)=1 with grids%mu= 0, i.e. zero expected equity premium'
            deallocate(grids%mu); allocate(grids%mu(1))
            grids%mu = 0.0
            print*,'setting simvars%mu= 0'
            do i=1,size(simvars)
                simvars(i)%mu = 0.0
            enddo
        endif

    else
        print*,'- run_model: Krusell-Smith GENERAL equilibrium'
        dir    = 'ge'
        allocate(simvars(OMP_get_max_threads()))
        call simvars%allocate(nt)
        call random_seed(put=seed) ! so that same sequence for different experiments
        do i=1,size(simvars)
            simvars(i)%z     = MarkovChain(pi_z,nt)
            simvars(i)%K(1)  = ms_grids%k(1)    ! starting value for simulation
            simvars(i)%mu(1) = ms_grids%mu(1)   ! starting value for simulation
            simvars(i)%rf(1) = ms_rf_temp       ! starting value for simulation
        enddo
        ! Not simvars%mu(1) = ms_grids%mu(1), because overwritten in simulations. Instead, calc mu0 from agg_grid.
        coeffs = Initialize(dir, n_coeffs,nz, estimate_from_simvars, ms_grids)
        call grids%construct(ms_grids,factor_k,factor_mu,cover_k, cover_mu, nk,nmu)
    endif
    output_path = construct_path(dir,calib_name)
    if (.not. calibrating) syserr = system('mkdir '//output_path)

    call solve_krusellsmith(grids, projectname, calib_name, output_path, it, coeffs, simvars, Phi, xgrid_ms, policies, value, lifecycles, err)
    if (err%not_converged) call err%print2stderr(dir)
    welfare = calc_average_welfare(simvars)

    ! Check distribution and save results
    print*, ' '
    if (.not. calibrating) call CheckPhi(Phi,output_path)
    if (.not. partial_equilibrium) call save_unformatted(grids, coeffs, simvars)
    if (.not. calibrating) call save_and_plot_results(dir, grids, err)
    if (present(simvars_o)) simvars_o = simvars
    if (.not. calibrating) then
        print*,' **** Completed solution of calibration ', calib_name, ' **** '
        print*, ' '
        print*, ' '
    else
        print*,' * Model solution succsessful during calibration of ', calib_name
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

        print*, '- run_model: Saving results and plots to folder ./model_output '
        call system_clock(end_time,count_rate)
        secs= real(end_time-start_time,dp)/real(count_rate,dp) ! excludes time for saving and plotting results

        call save_results(Phi, simvars, coeffs, grids,lifecycles,&
                             policies, secs, it, projectname, calib_name, dir, err)

        call plot_results(output_path, 'plot_all')

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
    real(dp) function inverted_average_mpk() !average_of_simulations()
        use statistics        ,only: tStats
        use params_mod, only: del_mean, de_ratio
        use income, only: f_net_mpk, alpha, delta, zeta
        type(tSimvars) ,allocatable :: s_temp(:)
        type(tStats)   :: k_stats, mpk_stats
        real(dp), dimension(nt) :: mpk
        real(dp) :: rf, r
        integer :: tc

        call read_unformatted(s_temp)

        !call k_stats%calc_stats(s_temp%K, s_temp%err_mu, s_temp%err_K)

!        do tc = 1,nt
!            mpk(tc) = f_net_mpk(s_temp%K(tc), zeta(s_temp%z(tc)), delta(s_temp%z(tc)))
!        enddo
!        call mpk_stats%calc_stats(mpk, s_temp%err_mu, s_temp%err_K)
!        inverted_average_mpk = (alpha/(mpk_stats%avg_()+del_mean))**(1.0/(1.0-alpha))

        rf = real(1.84343749E-02,dp)
        r  = real(7.48861757E-02,dp)
        mpk(1) = (r+rf*de_ratio)/(1.0+de_ratio)
        mpk(1) = real(3.07358605E-02,dp)
        !real(1.23087113E-01,dp) ! This is the value for pfrw: portfolio return unweighted
        !real(1.06801078E-01,dp) ! this is the value for pfr: portfolio return (with pop weights)
        inverted_average_mpk = (alpha/(mpk(1)+del_mean))**(1.0/(1.0-alpha))

    end function inverted_average_mpk !average_of_simulations
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
        use params_mod  ,only: n_coeffs, nz
        type(tAggGrids) ,intent(out) :: grids
        type(tCoeffs)   ,intent(out) :: coeffs
        type(tSimvars)  ,allocatable ,intent(out) :: simvars(:)

        call grids%read_unformatted('ks')
        call coeffs%read_unformatted

        call read_unformatted(simvars)
    end subroutine read_unformatted_ks
!-------------------------------------------------------------------------------

end subroutine run_model

end module run_model_mod
