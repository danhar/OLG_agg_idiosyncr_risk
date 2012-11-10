module run_model_mod
    implicit none
contains

subroutine run_model(projectname, calib_name, welfare)
    use types
    use policies_class
    use error_class
    use ifport            ,only: system  ! Intel Fortran portability library
    use unformatted_io    ,only: SaveUnformatted, ReadUnformatted
    use aggregate_grids
    use laws_of_motion    ,only: tCoeffs, Initialize, MakeVector, MakeType
	use markovchain_mod   ,only: MarkovChain
	use distribution      ,only: CheckPhi
	use mean_shock_mod    ,only: solve_meanshock
    use krusell_smith_mod ,only: solve_krusellsmith
	use params_mod        ,only: construct_path, set_apmax, & ! procedures
	                             partial_equilibrium, estimate_from_simvars, save_all_iterations, & ! logicals
	                             nk,nmu, nz, n_coeffs, nt, ms_guess, factor_k, factor_mu,cover_k, cover_mu, pi_z, seed, scale_AR

	character(len=*) ,intent(in)  :: projectname, calib_name
	real(dp)         ,intent(out) :: welfare ! expected ex-ante utility of a newborn
    type(tCoeffs)     :: coeffs  ! coefficients for laws of motion
    type(tPolicies)   :: policies
    type(tAggGrids)   :: grids, ms_grids
    type(tLifecycle)  :: lifecycles
    type(tSimvars)    :: simvars    ! simulation variables, e.g. zt, kt,...
    type(tErrors)     :: err
    real(dp),allocatable :: value(:,:,:,:,:,:), Phi(:,:,:), xgrid_ms(:,:,:) ! value function, distribution, and mean shock xgrid
	integer           :: start_time, it, syserr ! 'it' cointains total iterations of Krusell-Smith loop
	character(:),allocatable :: dir, output_path

    call system_clock(start_time)
!-------------------------------------------------------------------------------
! Mean Shock (partial or general) Equilibrium (to generate good initial guesses)
!-------------------------------------------------------------------------------
    print*, ' '
    if (partial_equilibrium) then
        print*,'- main: mean shock PARTIAL equilibrium'
        dir    = 'mspe'
        call ReadUnformatted(ms_grids)
        if (scale_AR == -1.0) then
            print*,'- main: setting ms_grids%mu = 0.0, ms_grids%k = average over simulations'
            ms_grids%mu =0.0
            ms_grids%k = inverted_average_mpk()
        endif
    else
        print*,'- main: mean shock GENERAL equilibrium'
        dir    = 'msge'
        ms_grids = ms_guess
    endif
    output_path = construct_path(dir,calib_name)
    syserr = system('mkdir '//output_path) ! Creates directory for output files
    it = 0  ! no Krusell-Smith iterations in Mean shock (but variable still needed for saving results)

    call solve_meanshock(coeffs, ms_grids, policies, simvars, lifecycles, Phi, xgrid_ms, value, err, output_path)
    if (err%not_converged) call err%print2stderr(dir)
    welfare = calc_average_welfare(simvars)

    ! Check distribution and save results
    call CheckPhi(Phi,output_path) ! writes errors to file
    if (.not. partial_equilibrium) then
        syserr = system('cp model_input/last_results/*.unformatted model_input/last_results/previous/')
        call SaveUnformatted(Phi, ms_grids, simvars%pens(1))
    endif
    call save_and_plot_results(dir, ms_grids, err)
    print*, ' '

    if (scale_AR == -1.0) then
        print*,' **** Completed solution of calibration ', calib_name, ' **** '
        print*, ' '
        print*, ' '
        return
    endif
!stop 'Stop after ms'
!-------------------------------------------------------------------------------
! Full (partial or general) Equilibrium (i.e. Krusell/Smith algorithm)
!-------------------------------------------------------------------------------
    ! set new apmax using the results of mean shock
    call set_apmax(ms_grids%k(1)*factor_k)

    if(partial_equilibrium) then
        print*,'- main: Krusell-Smith PARTIAL equilibrium'
        dir    = 'pe'
        call ReadUnformatted(grids, coeffs, simvars)
        if (scale_AR == -1.0) then
            ! This is never executed at the moment because of the conditional return in line 70
            print*,'scale_AR = -1.0, i.e. no aggregate risk'
            print*,'setting size(grids%mu)=1 with grids%mu= 0, i.e. zero expected equity premium'
            deallocate(grids%mu); allocate(grids%mu(1))
            grids%mu = 0.0
            print*,'setting simvars%mu= 0'
            simvars%mu = 0.0
        endif

    else
        print*,'- main: Krusell-Smith GENERAL equilibrium'
        dir    = 'ge'
        ! realizations of aggregate shock z for simulation (keep constant over all simulations)
        call AllocateType(simvars,nt)
        simvars%z     = MarkovChain(pi_z,nt, seed=seed)
        simvars%K(1)  = ms_grids%k(1) ! starting value for simulation
        ! Not simvars%mu(1) = ms_grids%mu(1), because overwritten in simulations. Instead, calc mu0 from agg_grid.
        coeffs = Initialize(dir, n_coeffs,nz, estimate_from_simvars, ms_grids)
        grids  = MakeGrid(ms_grids,factor_k,factor_mu,cover_k, cover_mu, nk,nmu)
    endif
    output_path = construct_path(dir,calib_name)
    syserr = system('mkdir '//output_path)

    ! The following block is for debugging / understanding Krusell Smith
    if (save_all_iterations) then   ! format the file where intermediate coeffs will be saved
        open(unit=132, file=output_path//'/loms_it.txt', status = 'replace')
        write(132,*)
        write(132,'(t8,a91)') 'coeffs_k(1) , coeffs_k(2) , coeffs_k(3)     ---    coeffs_mu(1), coeffs_mu(2), coeffs_mu(3)'
        write(132,*)
        close(132)
    endif

    call solve_krusellsmith(grids, projectname, calib_name, output_path, it, coeffs, simvars, Phi, xgrid_ms, policies, value, lifecycles, err)
    if (err%not_converged) call err%print2stderr(dir)
    welfare = calc_average_welfare(simvars)

    ! Check distribution and save results
    print*, ' '
    call CheckPhi(Phi,output_path)
    if (.not. partial_equilibrium) call SaveUnformatted(grids, coeffs, simvars)
    call save_and_plot_results(dir, grids, err)
    print*,' **** Completed solution of calibration ', calib_name, ' **** '
    print*, ' '
    print*, ' '

    ! Solution of mean shock OLG model during transition ?
    ! Solution of model with stochastic transition ?

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - subroutine save_and_plot_results(dir, grids, err)
! - pure real(dp) function calc_average(simvars%welfare)
! - real(dp) function average_of_simulations()
!-------------------------------------------------------------------------------

    subroutine save_and_plot_results(dir, grids, err)
    ! Calc stats and save all results for each run, then plot and save to pdf
        use save_results_mod
	    character(len=*) ,intent(in)   :: dir
	    type(tAggGrids)  ,intent(in)   :: grids
        type(tErrors)    ,intent(in)   :: err
        real(dp) :: secs
        integer  :: end_time, count_rate

        print*, '- main: Saving results and plots to folder ./model_output '
        call system_clock(end_time,count_rate)
        secs= real(end_time-start_time,dp)/real(count_rate,dp) ! excludes time for saving and plotting results

        call save_results(Phi, simvars, coeffs, grids,lifecycles,&
                             policies, secs, it, projectname, calib_name, dir, err)

        call plot_results(output_path, 'plot_all')

    end subroutine save_and_plot_results

!-------------------------------------------------------------------------------
    pure real(dp) function calc_average_welfare(simvars)
	    use statistics        ,only: tStats
	    type(tSimvars), intent(in) :: simvars
	    type(tStats)      :: welfare_stats
	    call welfare_stats%calc_stats(simvars%welf, simvars%err_mu, simvars%err_K)
	    calc_average_welfare = welfare_stats%avg_exerr_() ! without the periods where simvars hit gridbounds
    end function calc_average_welfare

!-------------------------------------------------------------------------------
    real(dp) function inverted_average_mpk() !average_of_simulations()
        use statistics        ,only: tStats
        use params_mod, only: del_mean, de_ratio
        use income, only: f_net_mpk, alpha, delta, zeta
        type(tSimvars) :: s_temp
        type(tStats)   :: k_stats, mpk_stats
        real(dp), dimension(nt) :: mpk
        real(dp) :: rf, r
        integer :: tc

        call AllocateType(s_temp,nt)
        open(55,file='model_input/last_results/simvars_ge.unformatted',form='unformatted',action='read')
        read(55) s_temp%z, s_temp%K, s_temp%mu, s_temp%B, s_temp%C, s_temp%Phi_1, s_temp%Phi_nx, s_temp%err_aggr, &
        s_temp%err_income, s_temp%r, s_temp%rf, s_temp%wage, s_temp%pens, s_temp%tau, s_temp%welf, s_temp%bequests, s_temp%err_K, s_temp%err_mu
        close(55)

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

end subroutine run_model
end module run_model_mod
