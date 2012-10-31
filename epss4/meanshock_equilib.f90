module meanshock_equilib

    use kinds
    use params_mod      ,only: n_eta, nx, nj, nz
    use policyfunctions ,only: tPolicies
    use aggregate_grids ,only: tAggGrids
    use laws_of_motion  ,only: tCoeffs
    use error_class      ,only: tErrors
    implicit none

    private

    real(dp), allocatable :: Phi(:,:,:)      ! distribution
    real(dp), allocatable :: value(:,:,:,:,:,:)
    type(tCoeffs)   :: coeffs                !coeffs of loms
    type(tPolicies) :: policies
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
    use params_mod    ,only: L_N_ratio, n, g, stat_dist_z, de_ratio
    use income
    use policyfunctions
    use distribution ,only: TransitionPhi

    implicit none
    real(dp) ,dimension(:),intent(in) :: msvars 			! k_ms, mu_ms
    real(dp) ,dimension(size(msvars)) :: distance           ! excess demands
    real(dp) ,dimension(nx,n_eta,nj)  :: apgrid_ms, stocks_ms, xgrid_ms	! policies /grids mean shock projection
    real(dp)						  :: kp_ms, agg_bond_demand
    real(dp)						  :: netwage_ms, pens_ms, r_ms, rf_ms
	integer 						  :: i

    grid%k  = msvars(1)
    grid%mu = msvars(2)

    call policies%solve(coeffs, grid, value, errs)

	! Projection of policies / grids on mean shock
    xgrid_ms =0.0; apgrid_ms =0.0; stocks_ms =0.0
	do i=1,nz
		xgrid_ms  = xgrid_ms  + w(i)* policies%xgrid (:,:,i,:,1,1)
		apgrid_ms = apgrid_ms + w(i)* policies%apgrid(:,:,i,:,1,1)
		stocks_ms = stocks_ms + w(i)* policies%stocks(:,:,i,:,1,1)
	enddo

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
pure subroutine ms_equilib_get(P, pol, v, er)
    real(dp), allocatable, intent(out)   :: P(:,:,:), v(:,:,:,:,:,:)
    type(tPolicies), intent(out)         :: pol
    type(tErrors), intent(inout)           :: er ! need inout because of er%not_converged
    P   = Phi
    pol = policies
    v = value
    ! The trouble with err is that I do not want to overwrite er%not_converged
    er%asset   = errs%asset
    er%cons    = errs%cons
    er%kp      = errs%kp
    er%mup      = errs%mup
    er%rfp      = errs%rfp
end subroutine ms_equilib_get

end module meanshock_equilib
