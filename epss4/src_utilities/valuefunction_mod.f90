module valuefunction_mod
    use kinds
    use classes_mod ,only: tAggGrids, tPolicies, tErrors, tCoeffs

    implicit none
    private
    public value

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - (pure) subroutine olg_backwards_recursion(p, coeffs, grids, value, err)
! - pure subroutine calc_vars_tomorrow(coeffs,grids,jc,zc,kc,muc,kp,mup,rp,rfp,yp, err_k, err_mu, err_rfp)
! - pure subroutine interp_policies_tomorrow(p,cons,value,kp, mup, grid, jc, consp, xgridp, vp, app_min)
! - pure function f_apgrid_j(rfp,yp, xgridp, app_min)
! - pure subroutine asset_allocation(xgridp, consp, vp, yp, rfp, rp, ap, pi_zp, pi_etap, xc, jc, kappa_out, error)
! - pure subroutine consumption(ap, kappa, xgridp, consp, vp, rfp,rp, yp, zc, xc, ec, betatildej, cons_out, evp, error)
!-------------------------------------------------------------------------------

pure function value(p, coeffs, grids)
! Get the policy functions for the entire state space, i.e. both individual and aggregate states
! This is the master subroutine, calling all module procedures below (which are appear in calling order)
! It it pure but for the OMP directives
    use params_mod      ,only: jr,surv, pi_z, pi_eta, trans_prob, cmin, g, beta, theta, gamm, apmax
    use makegrid_mod
    real(dp)            ,allocatable :: value(:,:,:,:,:,:)  ! could make optional
    type(tPolicies)                  ,intent(in) :: p      ! policies
    type(tCoeffs)                    ,intent(in)  :: coeffs ! coefficients for laws of motion
    type(tAggGrids)                  ,intent(in)  :: grids  ! grids for aggregate states k and mu

    real(dp) ,dimension(:,:,:,:,:,:) ,allocatable :: cons
    real(dp) ,dimension(:,:,:) ,allocatable             :: xgridp, consp, vp   ! xgrid, consumption, and value function tomorrow
    real(dp) ,dimension(:,:,:) ,allocatable                :: yp      ! income tomorrow, for every idiosyncr & aggr state
    real(dp) ,dimension(:), allocatable                       :: rp, mup ! risky return AND equity premium for every aggr state tomorrow
    real(dp)   :: kp, rfp, app_min                      ! tomorrow's capital, equity premium,risk-free rate, min aprime
    real(dp)   :: betatildej, evp_xc                            ! modified discount factor, expected value tomorrow
    integer    :: nx,n_eta,n_trans,nz, nj, nk, nmu, jc, muc, kc, zc, xc, ec

    nx = size(p%apgrid,1)
    n_eta = size(p%apgrid,2)
    n_trans=size(trans_prob)
    nz = size(p%apgrid,3)
    nj = size(p%apgrid,4)
    nk = size(p%apgrid,5)
    nmu= size(p%apgrid,6)
    allocate(yp(n_trans, n_eta,nz), rp(nz), mup(nz))
    allocate(xgridp(nx,n_eta,nz), consp(nx,n_eta,nz), vp(nx,n_eta,nz))
    allocate(cons(nx,n_eta,nz,nj,nk,nmu), value(nx,n_eta,nz,nj,nk,nmu))

    cons   = p%consumption()
    value(:,:,:,nj,:,:)   = cons(:,:,:,nj,:,:)

    do jc= nj-1,1,-1
        betatildej = beta*surv(jc)**(1.0/gamm)*(1.0+g)**((1.0-theta)/gamm)
        do kc=1,nk          ! Small performance notice: the outermost loop does not correspond to the rightmost state, because I interchanged loops for k and mu, so that OpenMP can work on k.
            do muc=1,nmu
                do zc=1,nz
                    call calc_vars_tomorrow(coeffs,grids,jc,zc,kc,muc,kp,mup,rp,rfp,yp)
                    call interp_policies_tomorrow(p, cons, value, kp, mup, grids, jc, consp, xgridp, vp, app_min)
                    do ec=1,n_eta
                        do xc=1,nx
                            evp_xc = evp(p%apgrid(xc,ec,zc,jc,kc,muc), p%kappa(xc,ec,zc,jc,kc,muc), xgridp, consp, vp, rfp,rp, yp, zc, xc, ec, betatildej)
                            value(xc,ec,zc,jc,kc,muc) = (cons(xc,ec,zc,jc,kc,muc)**((1.0-theta)/gamm) + betatildej*evp_xc**(1.0/gamm))**(gamm/(1.0-theta))
                        end do
                    enddo
                enddo
            enddo
        enddo
    enddo

end function value
!-------------------------------------------------------------------------------

pure subroutine calc_vars_tomorrow(coeffs,grids,jc,zc,kc,muc,kp,mup,rp,rfp,yp)
! Forecast kp, mup, and get corresonding prices  (might want to move into laws_of_motion)
    use params_mod     ,only: nz, pi_z, jr, ej, etagrid, trans_grid
    use laws_of_motion ,only: Forecast_k, Forecast_mu
    use income

    type(tCoeffs)   ,intent(in)  :: coeffs ! coefficients for laws of motion
    type(tAggGrids) ,intent(in)  :: grids  ! grids for aggregate states k and mu
    integer         ,intent(in)  :: jc, zc, kc, muc
    real(dp)        ,intent(out) :: kp, mup(:), rp(:), rfp, yp(:,:,:)
    integer :: zpc, epc, nk, nmu
    real(dp), parameter :: crit = 1.0e-10

    nk = size(grids%k)
    nmu= size(grids%mu)

    kp  = Forecast_k(coeffs%k(:,zc), grids%k(kc), grids%mu(muc))
    if (kp  > grids%k(nk)) then
        kp  = grids%k(nk)
    elseif (kp  < grids%k(1)) then
        kp  = grids%k(1)
    endif

    do zpc = 1,nz
        mup(zpc) = Forecast_mu(coeffs%mu(:,zpc), kp, grids%mu(muc))
    enddo
    where (mup > grids%mu(nmu)) mup = grids%mu(nmu) ! obsolete comment: This takes care of the wrong forecasts in the mean shock equilibrium
    where (mup < grids%mu(1))   mup = grids%mu(1)

    ! calculate tomorrow's risky returns and wage for given law of motion, for every zc today
    rfp = f_riskfree_rate(kp,grids%mu(muc),pi_z(zc,:))
    do zpc= 1,nz
        rp(zpc) = f_stock_return(kp, zeta(zpc), delta(zpc), rfp)
        if (jc+1>=jr) then
            yp(:,:,zpc) = f_pensions(kp, zeta(zpc))
        else
            do epc=1,size(etagrid,1)
                yp(:,epc,zpc) = ej(jc+1) * f_netwage(kp, zeta(zpc)) * etagrid(epc,zpc) * trans_grid + f_transfers(kp, zeta(zpc))
            enddo
        endif
    enddo
    if (rfp < rp(1)*(1.0 + sign(0.0001_dp,rp(1))) .and. grids%mu(muc)>0.0 ) then
        rfp = rp(1)*(1.0 + sign(0.0001_dp,rp(1)))
    endif

end subroutine calc_vars_tomorrow
!-------------------------------------------------------------------------------

pure subroutine interp_policies_tomorrow(p,cons,value,kp, mup, grid, jc, consp, xgridp, vp, app_min)
! Make a projection of tomorrow's policies on k prime and mu prime
    use fun_locate

    type(tPolicies) ,intent(in) :: p
    type(tAggGrids)  ,intent(in) :: grid
    real(dp), intent(in) :: kp, mup(:), cons(:,:,:,:,:,:), value(:,:,:,:,:,:)
    integer, intent(in) :: jc
    real(dp), dimension(:,:,:) ,intent(out) :: consp, xgridp, vp
    real(dp)                         ,intent(out) :: app_min
    integer :: iK, imu, zpc, nz
    real(dp) :: wK, wmu

    nz = size(cons,3)

    if(size(grid%K)>1 .and. size(grid%mu)>1) then
        ! Projection of policies on Kt
        iK        = f_locate(grid%K, kp)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        wK        = (kp - grid%k(iK)) / (grid%k(iK+1) - grid%k(iK))
        do zpc=1,nz
            imu        = f_locate(grid%mu, mup(zpc))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
            wmu        = (mup(zpc) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))

            ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds
            consp(:,:,zpc) = (1-wK)*(1-wmu)*   cons(:,:,zpc,jc+1,iK  ,imu  ) + &
                                wK *(1-wmu)*   cons(:,:,zpc,jc+1,iK+1,imu  ) + &
                             (1-wK)*   wmu *   cons(:,:,zpc,jc+1,iK  ,imu+1) + &
                                wK *   wmu *   cons(:,:,zpc,jc+1,iK+1,imu+1)

            xgridp(:,:,zpc)= (1-wK)*(1-wmu)*p%xgrid(:,:,zpc,jc+1,iK  ,imu  ) + &
                                wK *(1-wmu)*p%xgrid(:,:,zpc,jc+1,iK+1,imu  ) + &
                             (1-wK)*   wmu *p%xgrid(:,:,zpc,jc+1,iK  ,imu+1) + &
                                wK *   wmu *p%xgrid(:,:,zpc,jc+1,iK+1,imu+1)

            vp(:,:,zpc)    = (1-wK)*(1-wmu)*  value(:,:,zpc,jc+1,iK  ,imu  ) + &
                                wK *(1-wmu)*  value(:,:,zpc,jc+1,iK+1,imu  ) + &
                             (1-wK)*   wmu *  value(:,:,zpc,jc+1,iK  ,imu+1) + &
                                wK *   wmu *  value(:,:,zpc,jc+1,iK+1,imu+1)
        enddo

        imu        = f_locate(grid%mu, mup(1))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        wmu        = (mup(1) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))
        app_min= (1-wK)*(1-wmu)*p%apgrid(1,1,1,jc+1,iK  ,imu  ) + & ! This creates smallest aprime for forecasts
                    wK *(1-wmu)*p%apgrid(1,1,1,jc+1,iK+1,imu  ) + & ! I should move this into the loop and calc for all zpc
                 (1-wK)*   wmu *p%apgrid(1,1,1,jc+1,iK  ,imu+1) + & ! coz then better to understand.
                    wK *   wmu *p%apgrid(1,1,1,jc+1,iK+1,imu+1)

    elseif (size(grid%K)>1) then
        ! Projection of policies on Kt
        iK        = f_locate(grid%K, kp)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        wK        = (kp - grid%k(iK)) / (grid%k(iK+1) - grid%k(iK))
        consp  = (1.0 -wK)*cons    (:,:,:,jc+1,iK,1) + wK*cons    (:,:,:,jc+1,iK+1,1)
        xgridp = (1.0 -wK)*p%xgrid (:,:,:,jc+1,iK,1) + wK*p%xgrid (:,:,:,jc+1,iK+1,1)
        vp     = (1.0 -wK)*value   (:,:,:,jc+1,iK,1) + wK*value   (:,:,:,jc+1,iK+1,1)
        app_min= (1.0 -wK)*p%apgrid(1,1,1,jc+1,iK,1) + wK*p%apgrid(1,1,1,jc+1,iK+1,1)

    elseif (size(grid%mu)>1) then
        do zpc=1,nz
            imu        = f_locate(grid%mu, mup(zpc))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
            wmu        = (mup(zpc) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))

            ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds
            consp(:,:,zpc) = (1.0-wmu)*   cons(:,:,zpc,jc+1,1,imu) + wmu*   cons(:,:,zpc,jc+1,1,imu+1)
            xgridp(:,:,zpc)= (1.0-wmu)*p%xgrid(:,:,zpc,jc+1,1,imu) + wmu*p%xgrid(:,:,zpc,jc+1,1,imu+1)
            vp(:,:,zpc)    = (1.0-wmu)*  value(:,:,zpc,jc+1,1,imu) + wmu*  value(:,:,zpc,jc+1,1,imu+1)
        enddo

        imu        = f_locate(grid%mu, mup(1))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        wmu        = (mup(1) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))
        app_min= (1.0-wmu)*p%apgrid(1,1,1,jc+1,1,imu) + wmu*p%apgrid(1,1,1,jc+1,1,imu+1)
        ! This creates smallest aprime for forecasts. should move this into the loop and calc for all zpc coz then better to understand.

    else    ! Mean shock
        consp  = cons    (:,:,:,jc+1,1,1)
        xgridp = p%xgrid (:,:,:,jc+1,1,1)
        vp     = value   (:,:,:,jc+1,1,1)
        app_min= p%apgrid(1,1,1,jc+1,1,1)
    endif

end subroutine interp_policies_tomorrow
!-------------------------------------------------------------------------------

pure function evp(ap, kappa, xgridp, consp, vp, rfp,rp, yp, zc, xc, ec, betatildej)
    use fun_lininterp
    use params_mod, only : collateral_constraint, pi_z, pi_eta, trans_prob, nz, n_eta, n_trans, g, cmin, theta, gamm

    real(dp)                    :: evp  ! cons(xc,ec,zc,jc), evp
    real(dp)                   ,intent(in)  :: ap, kappa, rfp, betatildej, rp(:),yp(:,:,:), xgridp(:,:,:),consp(:,:,:), vp(:,:,:)
    integer                    ,intent(in)  :: zc, xc, ec
    real(dp)   ,dimension(n_eta)            :: cons_interp, vp_interp ! interpolated values of cons, vp
    real(dp)   ,dimension(nz)               :: evpz, rhs_temp ! temporary: cee_*: consumption euler eq.
    real(dp)                                :: rtildep, xp    ! rtilde and cash-at-hand tomorrow
    real(dp)                                :: cee_rhs, rhs_fac1, rhs_fac2 ! temporary: rhs of cee, factors 1 and 2
    integer                                 :: zpc, epc, tpc       ! counters for shocks tomorrow

    evpz     = 0.0
    rhs_temp = 0.0

    do zpc=1,nz
        if (pi_z(zc,zpc) == 0.0) cycle

        ! Get cash-at-hand tomorrow for all states
        rtildep= (1.0+rfp+kappa*(rp(zpc)-rfp))/(1.0+g)
        do tpc = 1,n_trans
            do epc = 1,n_eta
                xp = yp(tpc,epc,zpc)+rtildep*ap

                ! Interpolate consumption and value function
                cons_interp(epc) = f_lininterp(xgridp(:,epc,zpc),consp(:,epc,zpc),xp)
                vp_interp(epc)   = f_lininterp(xgridp(:,epc,zpc),vp(:,epc,zpc),xp)
                if (cons_interp(epc) <=cmin) then
                    cons_interp(epc) = cmin
                    xp               = cmin + ap
                    vp_interp(epc)   = f_lininterp(xgridp(:,epc,zpc),vp(:,epc,zpc),xp)
                endif
                if (vp_interp(epc) <=cmin) then
                    vp_interp(epc)   = cmin
                endif

            enddo
            ! Now different factors of consumption euler equation
            evpz(zpc)     = evpz(zpc) + trans_prob(tpc)* sum(pi_eta(ec,:)*vp_interp**(1.0-theta))
        enddo
    enddo

    evp=dot_product(pi_z(zc,:),evpz)    ! Expected V'

end function evp
!-------------------------------------------------------------------------------

end module valuefunction_mod
