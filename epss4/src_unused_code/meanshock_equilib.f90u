module meanshock_equilib

    use kinds
    use household_solution_mod ,only: tPolicies, olg_backwards_recursion
    use aggregate_grids        ,only: tAggGrids
    use laws_of_motion         ,only: tCoeffs
    use error_class            ,only: tErrors
    implicit none

    private

    real(dp), allocatable :: Phi(:,:,:)      ! distribution
    real(dp), allocatable :: v_fine(:,:,:,:,:,:)
    real(dp) ,dimension(:,:,:) ,allocatable  :: apgrid_ms, stocks_ms, xgrid_ms, kappa_ms, value_ms  ! policies /grids mean shock projection
    type(tCoeffs)   :: coeffs                !coeffs of loms
    type(tPolicies) :: policies, fine
    type(tAggGrids) :: grid
    type(tErrors)   :: errs
    real(dp)        :: mean_zeta, mean_delta ! mean shocks
    real(dp), allocatable :: w(:)                 ! weights (distance from the two means)
    real(dp), allocatable :: etagrid(:)        ! idio states in mean shock path

    public ms_equilib, ms_equilib_set, ms_equilib_get

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - function ms_equilib(msvars) result(distance)
! - subroutine ms_equilib_set(co, ag_gr, m_z, m_d, s_w, eta)
! - pure subroutine ms_equilib_get(P, pol, v, e_o)
!-------------------------------------------------------------------------------

function ms_equilib(msvars) result(distance)
! Solve for the 'mean shock equilibrium' (ms) to get initial values for k, mu, Phi
    use params_mod     ,only: L_N_ratio, n, g, stat_dist_z, de_ratio, nx_factor
    use income
    use policies_class ,only: tPolicies
    use distribution   ,only: TransitionPhi
    use interpolate_xgrid ,only: InterpolateXgrid

    implicit none
    real(dp) ,dimension(:),intent(in) :: msvars 			! k_ms, mu_ms
    real(dp) ,dimension(size(msvars)) :: distance           ! excess demands
    real(dp) ,dimension(:,:,:,:,:,:), allocatable :: value
    real(dp)						  :: kp_ms, agg_bond_demand
    real(dp)						  :: netwage_ms, pens_ms, r_ms, rf_ms
	integer 						  :: i

    grid%k  = msvars(1)
    grid%mu = msvars(2)

    call olg_backwards_recursion(policies,coeffs, grid, value, errs)

    call InterpolateXgrid(nx_factor, policies, value, fine, v_fine)

    if (.not. allocated(apgrid_ms)) &
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
        kappa_ms = stocks_ms/apgrid_ms  ! I don't need it right here, but in meanshock_wrapper
    elsewhere
        kappa_ms = 0.0
    end where

	! Prices in mean shock path
	netwage_ms	  = f_netwage (msvars(1), mean_zeta)
	pens_ms		  = f_pensions(msvars(1), mean_zeta)
    rf_ms         = f_riskfree_rate(msvars(1),msvars(2),stat_dist_z)
	r_ms		  = f_stock_return(msvars(1), mean_zeta, mean_delta, rf_ms)

	! Get distribution in mean shock path
	Phi	= TransitionPhi(rf_ms,r_ms,netwage_ms,pens_ms,xgrid_ms,apgrid_ms,stocks_ms,etagrid)

    ! Aggregate

	kp_ms           = sum(apgrid_ms *Phi)/(L_N_ratio*(1.0+n)*(1.0+g))
	agg_bond_demand = sum((apgrid_ms-stocks_ms)*Phi)
	! Exess demands
	distance(1)     = kp_ms - msvars(1)
	distance(2)     = kp_ms * de_ratio/(1.0 + de_ratio) - agg_bond_demand/(L_N_ratio*(1.0+n)*(1.0+g))

end function ms_equilib

!-------------------------------------------------------------------------------
subroutine ms_equilib_set(co, ag_gr, m_z, m_d, s_w, eta)
    type(tCoeffs), intent(in)         :: co
    type(tAggGrids), intent(in)       :: ag_gr
    real(dp), intent(in)              :: m_z, m_d
    real(dp),dimension(:), intent(in) :: s_w, eta
    coeffs     = co
    grid       = ag_gr
    mean_zeta  = m_z
    mean_delta = m_d
    w          = s_w
    etagrid    = eta
end subroutine ms_equilib_set

!-------------------------------------------------------------------------------
pure subroutine ms_equilib_get(P, pol, v_fi, er, xgr, ap, sto, kap, v_ms)
    real(dp), allocatable, intent(out)   :: P(:,:,:)
    real(dp), dimension(:,:,:,:,:,:) ,allocatable, intent(out):: v_fi
    real(dp), dimension(:,:,:)       ,allocatable, intent(out):: xgr, ap, sto, kap, v_ms
    type(tPolicies), intent(out)         :: pol
    type(tErrors), intent(inout)           :: er ! need inout because of er%not_converged
    P   = Phi
    pol = fine !policies
    v_fi = v_fine
    v_ms = value_ms
    xgr = xgrid_ms
    ap = apgrid_ms
    sto = stocks_ms
    kap = kappa_ms
    ! The trouble with err is that I do not want to overwrite er%not_converged
    er%asset   = errs%asset
    er%cons    = errs%cons
    er%kp      = errs%kp
    er%mup      = errs%mup
    er%rfp      = errs%rfp
end subroutine ms_equilib_get

end module meanshock_equilib
