module insurance_effect_mod
! In this module, the pure insurance effect is computed, excluding the mean effect.
! This module was added very late and implements experiments that are hard to insert into the existing code base.
! As a consequence, a lot of code is duplicated here.
    implicit none
contains

subroutine calc_insurance_effect(policies, value, agg_grid, simvars_in, Phi_in, calib_name, projectname, welfare)
    use ifport           ,only: system  ! Intel Fortran portability library
    use kinds            ,only: dp
    use classes_mod      ,only: tSimvars, tLifecycle, tPolicies, tAggGrids, tStats
    use lifecycles_class ,only: average
    use simvars_class    ,only: print_error
    use fun_locate
    use params_mod       ,only: exogenous_xgrid, nx_factor, stat_dist_z, stat_dist_eta, ccv, scale_IR, scale_AR, scale_IR_orig, etagrid, &
                                params_set, params_set_thisrun, CheckParams, construct_path
    use simulation_mod   ,only: simulate
    use interpolate_xgrid

    type(tPolicies)  ,intent(in)  :: policies
    real(dp)         ,intent(in)  :: value(:,:,:,:,:,:)
    type(tAggGrids)  ,intent(in)  :: agg_grid
    type(tSimvars)   ,intent(in)  :: simvars_in(:)    ! (zt, kt, mut, bt,...), first element contains starting values
    real(dp)         ,intent(in)  :: Phi_in(:,:,:) ! distribution. Returns: average Phi if (exogenous_xgrid), else Phi in nt
    character(len=*) ,intent(in)  :: calib_name, projectname
    real(dp)         ,intent(out) :: welfare(5)
    type(tLifecycle)              :: lifecycles, lifecycles_array(size(simvars_in))         ! lifecycle profiles
    type(tPolicies)               :: pol_fine, pol_minus_risk    ! if (exogenous_xgrid) then this will hold interpolated policies
    type(tStats)                  :: K, mu, netwage, pens, bequests, r, rf
    type(tSimvars) ,allocatable   :: simvars(:)
    real(dp)       ,allocatable   :: Phi(:,:,:)
    real(dp), allocatable     :: val_newx(:,:,:,:,:,:), val_minus_risk(:,:,:,:,:,:), xgrid_mean_old(:,:,:), xgrid_mean_new(:,:,:), Phi_spread(:,:,:,:)
    character(len=7)          :: runchar
    character(:),allocatable  :: output_path
    logical                   :: ccv_thisrun
    real(dp)                  :: scale_IR_thisrun, scale_AR_thisrun, w
    integer                   :: syserr, i

    output_path = construct_path(calib_name,'insurance')

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

    print *,'- insurance_effect: simulating'

    call InterpolateXgrid(nx_factor, policies, value, pol_fine, val_newx)
    xgrid_mean_old = sum(pol_fine%xgrid(:,:,1,:,:,1),4)/size(pol_fine%xgrid,5) ! only an approximation of the grid over which Phi is defined

    write(runchar,'(a6)') ',all'
    pol_minus_risk = policies
    val_minus_risk = value
    call sim_pe(welfare(1))

    call params_set('ccv', .false.)
    write(runchar,'(a6)') ',noCCV'
    pol_minus_risk = policies
    val_minus_risk = value
    ! not possible to interpolate policy functions, since ccv doesn't have its own dimension.
    ! So here only remove ccv from income in simulation, use the re-optimized policies and adjust them with aggregate consumption for behavioral response.
    call sim_pe(welfare(2))

    call params_set('scale_IR', -1.0_dp)
    pol_minus_risk = policies%mean(2,stat_dist_eta)
    val_minus_risk = value
    val_minus_risk = 0.0
    do i=1,size(stat_dist_eta)
        val_minus_risk(:,1,:,:,:,:)= val_minus_risk(:,1,:,:,:,:) + stat_dist_eta(i)*value(:,i,:,:,:,:)
    enddo
    val_minus_risk(:,2:,:,:,:,:) = spread(val_minus_risk(:,1,:,:,:,:),2, size(stat_dist_eta)-1)

    write(runchar,'(a5)') ',noIR'
    call sim_pe(welfare(3))

    call params_set('scale_AR', -1.0_dp)
    pol_minus_risk = policies%mean(3,stat_dist_z)
    pol_minus_risk = pol_minus_risk%interpolate(5,agg_grid%k , K%avg_exerr_())
    pol_minus_risk = pol_minus_risk%interpolate(6,agg_grid%mu,mu%avg_exerr_())
    val_minus_risk = value
    val_minus_risk = 0.0
    do i=1,size(stat_dist_z)
        val_minus_risk(:,:,1,:,:,:)= val_minus_risk(:,:,1,:,:,:) + stat_dist_z(i)*value(:,:,i,:,:,:)
    enddo
    val_minus_risk(:,:,2:,:,:,:) = spread(val_minus_risk(:,:,1,:,:,:),3, size(stat_dist_z)-1)

    i=f_locate(agg_grid%k,K%avg_exerr_())
    w=(K%avg_exerr_()-agg_grid%k(i))/(agg_grid%k(i+1)-agg_grid%k(i))
    val_minus_risk(:,:,:,:,:,:) = spread((1-w)*val_minus_risk(:,:,:,:,i,:) + w*val_minus_risk(:,:,:,:,i+1,:),5,size(agg_grid%k))

    i=f_locate(agg_grid%mu,mu%avg_exerr_())
    w=(mu%avg_exerr_()-agg_grid%mu(i))/(agg_grid%mu(i+1)-agg_grid%mu(i))
    val_minus_risk(:,:,:,:,:,:) = spread((1-w)*val_minus_risk(:,:,:,:,:,i) + w*val_minus_risk(:,:,:,:,:,i+1),6,size(agg_grid%mu))

    if (scale_IR_orig .ne. -1.0) call params_set('scale_IR', 0.0_dp)
    write(runchar,'(a5)') ',noAR'
    call sim_pe(welfare(4))

    call params_set('scale_AR', -1.0_dp)
    call params_set('scale_IR', -1.0_dp)
    pol_minus_risk = pol_minus_risk%mean(2,stat_dist_eta)
    do i=1,size(stat_dist_eta)
        val_minus_risk(:,1,:,:,:,:)= val_minus_risk(:,1,:,:,:,:) + stat_dist_eta(i)*value(:,i,:,:,:,:)
    enddo
    val_minus_risk(:,2:,:,:,:,:) = spread(val_minus_risk(:,1,:,:,:,:),2, size(stat_dist_eta)-1)

    write(runchar,'(a7)') ',norisk'
    call sim_pe(welfare(5))


    call params_set('ccv', ccv_thisrun)
    call params_set('scale_IR', scale_IR_thisrun)
    call params_set('scale_AR', scale_AR_thisrun)
    call params_set_thisrun

    contains

    subroutine sim_pe(welf)
        use distribution           ,only: TransitionPhi
        real(dp) ,intent(out) ::welf
        output_path=construct_path(calib_name,'insurance')//runchar
        syserr = system('mkdir '//output_path//' > /dev/null 2>&1') ! Creates directory for output files, suppresses error if dir exists
        Phi     = Phi_in    ! initial values, makes clear that Phi is not output
        simvars = simvars_in

        call params_set_thisrun
        call CheckParams
        if (runchar == ',noCCV' .or. runchar==',noIR' .or. runchar==',all') then
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

            call InterpolateXgrid(nx_factor, pol_minus_risk, val_minus_risk, pol_fine, val_newx)

            ! Alternatively, one could take Phi from the noCCV simulations, which is an average, but I think the following is more correct
            Phi = TransitionPhi(rf%avg_exerr_(),r%avg_exerr_(),netwage%avg_exerr_(),pens%avg_exerr_(),bequests%avg_exerr_(),pol_fine%xgrid(:,:,1,:,1,1),pol_fine%apgrid(:,:,1,:,1,1),pol_fine%stocks(:,:,1,:,1,1),etagrid(:,1))

            welf = sum(val_newx(:,:,1,1,1,1)*Phi(:,:,1))

        endif

        !call save_and_plot_results(runchar,agg_grid)
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

    subroutine save_and_plot_results(dir, grids)
    ! Calc stats and save all results for each run, then plot and save to pdf
        use coefficients_class ,only: tCoeffs
        use error_class        ,only: tErrors
        use save_results_mod
        character(len=*) ,intent(in)   :: dir
        type(tAggGrids)  ,intent(in)   :: grids
        type(tErrors)    :: err_temp
        type(tCoeffs) :: coeffs_temp
        real(dp)      :: secs_temp
        integer       :: it_temp

        print*, '- insurance_effects_mod: Saving results and plots to folder ./model_output '
        secs_temp = 0.0
        call coeffs_temp%allocate(size(pol_minus_risk%xgrid,5), size(pol_minus_risk%xgrid,6), size(pol_minus_risk%xgrid,3))
        call    err_temp%allocate(size(pol_minus_risk%xgrid,5), size(pol_minus_risk%xgrid,6))

        call save_results(Phi, simvars, coeffs_temp, grids,lifecycles,&
                             pol_minus_risk, secs_temp, it_temp, projectname, calib_name//',insurance', dir, err_temp)

        call plot_results(output_path, 'plot_all')

    end subroutine save_and_plot_results


end subroutine calc_insurance_effect

end module insurance_effect_mod
