module insurance_effect_mod
    implicit none
contains

subroutine calc_insurance_effect(policies, value, agg_grid, simvars, Phi, lc, welfare)
    use kinds            ,only: dp
    use classes_mod      ,only: tSimvars, tLifecycle, tPolicies, tAggGrids
    use lifecycles_class ,only: average
    use simvars_class    ,only: print_error
    use params_mod       ,only: exogenous_xgrid, nx_factor, stat_dist_z, stat_dist_eta
    use simulation_mod   ,only: simulate
    use interpolate_xgrid

    type(tPolicies)  ,intent(in)    :: policies
    real(dp)         ,intent(in)    :: value(:,:,:,:,:,:)
    type(tAggGrids)  ,intent(in)    :: agg_grid
    type(tSimvars)   ,intent(inout) :: simvars    ! (zt, kt, mut, bt,...), first element contains starting values
    real(dp)         ,intent(inout) :: Phi(:,:,:) ! distribution. Returns: average Phi if (exogenous_xgrid), else Phi in nt
    type(tLifecycle) ,intent(out)   :: lc         ! lifecycle profiles
    real(dp)         ,intent(out)   :: welfare
    type(tPolicies) :: pol_newx    ! if (exogenous_xgrid) then this will hold interpolated policies
    type(tLifecycle) :: lifecycles_array(size(simvars))
    real(dp), allocatable :: val_newx(:,:,:,:,:,:), xgrid_mean_old(:,:,:), Phi_spread(:,:,:,:)


    print *,'- insurance_effect: simulating'


    contains

    subroutine sim_pe()
        if (exogenous_xgrid) then
            ! This is the standard case which should always be used, because we make the xgrid much finer
            call InterpolateXgrid(nx_factor, policies, value, pol_newx, val_newx)
            ! We also want the initial Phi that we take from previous simulations to be defined over the new xgrid
            xgrid_mean_old = xgrid_mean_new ! Would be nicer to have a derived type Phi which carries its own xgrid.
            xgrid_mean_new = sum(pol_newx%xgrid(:,:,1,:,:,1),4)/size(pol_newx%xgrid,5) ! only an approximation of the grid over which Phi is defined
            call InterpolateXgrid(Phi, xgrid_mean_old, xgrid_mean_new)
            Phi_spread = spread(Phi,4,size(simvars))
            !$OMP  PARALLEL DO IF(size(simvars)>1)
            do i=1,size(simvars)
                call simulate(pol_newx, val_newx, grids, simvars(i), Phi_spread(:,:,:,i), lifecycles_array(i))
            enddo
            !$OMP END PARALLEL DO
        else
            ! In this case, we simultaneously solve for the rf(t+1) (actually the mu(t+1)) and the corresponding xgrid(since it depends on mu).
            ! That is very costly, so we do not refine xgrid. While more correct theoretically, the coarse xgrid makes the solution less precise.
            Phi_spread = spread(Phi,4,size(simvars))
            !$OMP  PARALLEL DO IF(size(simvars)>1)
            do i=1,size(simvars)
                call simulate(policies, value, grids, simvars(i), Phi_spread(:,:,:,i), lifecycles_array(i))
            enddo
            !$OMP END PARALLEL DO
        endif
        lifecycles=average(lifecycles_array)
        Phi = sum(Phi_spread,4)/size(Phi_spread,4)
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

end subroutine calc_insurance_effect

end module insurance_effect_mod
