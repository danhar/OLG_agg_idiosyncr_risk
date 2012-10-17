module fun_excessbonds
    use kinds           ,only: dp
    use aggregate_grids ,only: tAggGrids
    implicit none

    private
    type(tAggGrids)                  :: agg_grid
    real(dp)                         :: nwaget, penst, rt, rft
    real(dp) ,dimension(:,:,:,:), allocatable :: apgrid_zk, kappa_zk, xgrid_zk  ! policies for given z and K
    real(dp) ,dimension(:,:,:), allocatable   :: apgridm, kappam  ! policies for given z, K, and mu last period (minus)
    real(dp) ,dimension(:,:,:), allocatable   :: Phim             ! distribution
    real(dp) ,dimension(:), allocatable       :: etagridt

    public f_excessbonds, set_excessbondsvars

    interface set_excessbondsvars
        module procedure set_excessbondsvars_end_xgrid, set_excessbondsvars_ex_xgrid
    end interface

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure function f_excessbonds(mut)
! - subroutine set_excessbondsvars(wage,pens,r,rf,apgrid_minus, kappa_minus, apgrid, kappa ,xgrid, aggr ,Ph, etagri)
! - pure subroutine get_excessbondsvars(apgrid, kappa, aggr ,Ph)
!-------------------------------------------------------------------------------

	pure real(dp) function f_excessbonds(mut)
	! calculates excess bond demand
	    use params_mod    ,only: L_N_ratio, ap_numzero, exogenous_xgrid, de_ratio, nx, nj, n_eta
	    use fun_locate
	    use distribution ,only: TransitionPhi
    !    use bs3vl_int   ! IMSL Math.pdf, p. 754ff: evaluate a 3d tensor-product spline given B-spline-coeffs

	    real(dp) ,intent(in)             :: mut       ! expected equity premium
        real(dp) ,dimension(nx,n_eta,nj) :: Phi    ! distribution
        real(dp) ,dimension(nx,n_eta,nj) :: apgridt, stockst, xgridt ! for given z, K, AND mu
	    real(dp)                         :: w, bond_supply
	    integer                          :: i

        if (size(agg_grid%mu) >1) then
		    i        = f_locate(agg_grid%mu, mut)
		    w        = (mut - agg_grid%mu(i)) / (agg_grid%mu(i+1) - agg_grid%mu(i))

	        if (.not. exogenous_xgrid) then
		        xgridt   = (1-w)* xgrid_zk(:,:,:,i) + w* xgrid_zk(:,:,:,i+1)
		        Phi      = TransitionPhi(rft,rt,nwaget,penst,xgridt,apgridm,kappam,etagridt, Phim)
	        else
	            Phi = Phim
	        endif

		    apgridt  = (1-w)*apgrid_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)
		    stockst  = (1-w)*apgrid_zk(:,:,:,i)* kappa_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)* kappa_zk(:,:,:,i+1)

        else
            Phi = Phim
            apgridt  = apgrid_zk(:,:,:,1)
            stockst  = apgrid_zk(:,:,:,1)* kappa_zk(:,:,:,1)
        endif

        bond_supply   = de_ratio * sum(stockst*Phi)/L_N_ratio
	    f_excessbonds = bond_supply - sum((apgridt - stockst)*Phi)/L_N_ratio    ! analytically, L_N_ratio drops out

	end function f_excessbonds
!-------------------------------------------------------------------------------

    subroutine set_excessbondsvars_end_xgrid(wage,pens,r,rf,apgrid_minus, kappa_minus, apgrid, kappa ,xgrid, aggr ,Ph, etagri)
    ! this one is called if (.not. exogenous_xgrid)
        real(dp)                     ,intent(in) :: wage, pens, r, rf
        real(dp), dimension(:,:,:)   ,intent(in) :: apgrid_minus, kappa_minus ! policies for given z, K, and mu last period
        real(dp), dimension(:,:,:,:) ,intent(in) :: apgrid, kappa, xgrid      ! policies for given z and K
        type(tAggGrids)              ,intent(in) :: aggr
        real(dp), dimension(:,:,:)   ,intent(in) :: Ph
        real(dp), dimension(:)       ,intent(in) :: etagri

        nwaget    = wage
        penst     = pens
        rt        = r
        rft       = rf
        apgridm   = apgrid_minus
        kappam    = kappa_minus
        apgrid_zk = apgrid
        kappa_zk  = kappa
        xgrid_zk  = xgrid
        agg_grid  = aggr
        Phim      = Ph
        etagridt  = etagri

    end subroutine set_excessbondsvars_end_xgrid
!-------------------------------------------------------------------------------

    subroutine set_excessbondsvars_ex_xgrid(apgrid, kappa, aggr ,Ph)
    ! this one is called if (exogenous_xgrid)
        real(dp), dimension(:,:,:,:) ,intent(in) :: apgrid, kappa   ! policies for given z and K
        type(tAggGrids)              ,intent(in) :: aggr
        real(dp), dimension(:,:,:)   ,intent(in) :: Ph
        apgrid_zk = apgrid
        kappa_zk  = kappa
        agg_grid  = aggr
        Phim      = Ph
    end subroutine set_excessbondsvars_ex_xgrid
!-------------------------------------------------------------------------------

    pure subroutine get_excessbondsvars(apgrid, kappa, aggr ,Ph)
        ! This getter method might be superfluous.
        real(dp) ,dimension(:,:,:,:) ,intent(out) :: apgrid, kappa    ! policies for given z and K
        type(tAggGrids)              ,intent(out) :: aggr
        real(dp) ,dimension(:,:,:)   ,intent(out) :: Ph

        apgrid = apgrid_zk
        kappa  = kappa_zk
        aggr   = agg_grid
        Ph     = Phim

    end subroutine get_excessbondsvars

end module fun_excessbonds
