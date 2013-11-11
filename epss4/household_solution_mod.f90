module household_solution_mod
    use kinds
    use classes_mod ,only: tAggGrids, tPolicies, tErrors, tCoeffs

    implicit none
    private
    public olg_backwards_recursion
    public interp_policies_tomorrow, consumption ! these are only needed for calculating the Euler errors in simulation_mod

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

subroutine olg_backwards_recursion(p, coeffs, grids, value, err)
! Get the policy functions for the entire state space, i.e. both individual and aggregate states
! This is the master subroutine, calling all module procedures below (which are appear in calling order)
! It it pure but for the OMP directives
    use params_mod      ,only: nj, nx, n_eta, nz, jr,surv, pi_z, pi_eta, cmin, g, beta, theta, gamm, apmax
    use makegrid_mod

    type(tPolicies)                  ,intent(out) :: p      ! policies
    type(tCoeffs)                    ,intent(in)  :: coeffs ! coefficients for laws of motion
    type(tAggGrids)                  ,intent(in)  :: grids  ! grids for aggregate states k and mu
    real(dp)            ,allocatable ,intent(out) :: value(:,:,:,:,:,:)  ! could make optional
    type(tErrors)                    ,intent(out) :: err
    real(dp) ,dimension(:,:,:,:,:,:) ,allocatable :: cons
    real(dp) ,dimension(nx,n_eta,nz)              :: xgridp, consp, vp   ! xgrid, consumption, and value function tomorrow
    real(dp) ,dimension(n_eta,nz)                 :: yp      ! income tomorrow, for every idiosyncr & aggr state
    real(dp) ,dimension(nz)                       :: rp, mup ! risky return AND equity premium for every aggr state tomorrow
    real(dp)   :: kp, rfp, app_min                      ! tomorrow's capital, equity premium,risk-free rate, min aprime
    real(dp)   :: betatildej, evp                            ! modified discount factor, expected value tomorrow
    integer    :: nk, nmu, jc, muc, kc, zc, xc, ec

    nk = size(grids%k)
    nmu= size(grids%mu)
    call p%allocate(nx,nz,nk,nmu)
    allocate(cons(nx,n_eta,nz,nj,nk,nmu), value(nx,n_eta,nz,nj,nk,nmu))
    call err%allocate(nk,nmu)

    !---------------------------------------------------------------------------
    ! Model solution, last generation
    !---------------------------------------------------------------------------
    do muc=1,nmu
        do kc=1,nk
            do zc=1,nz
                do ec=1,n_eta
                    p%xgrid(:,ec,zc,nj,kc,muc)=MakeGrid(cmin,apmax(ec,zc,nj),nx, 2.0_dp)
                enddo
            enddo
        enddo
    enddo
    ! Final period policy functions and value function
    p%apgrid(:,:,:,:,:,:) = 0.0
    p%kappa(:,:,:,:,:,:)  = 0.0
    p%stocks(:,:,:,:,:,:) = 0.0
    cons(:,:,:,nj,:,:)    = p%xgrid(:,:,:,nj,:,:)
    value(:,:,:,nj,:,:)   = cons(:,:,:,nj,:,:)

    !---------------------------------------------------------------------------
    ! Model solution, generations nj-1 to 1
    !---------------------------------------------------------------------------
!$OMP PARALLEL IF(nk>1) DEFAULT(NONE) &
!$OMP SHARED(p,value,cons,grids,coeffs,err,nmu,nk,nz,nj,n_eta,nx,beta,g,theta,gamm,surv, pi_eta, pi_z, apmax) &
!$OMP PRIVATE(jc,muc,kc,zc,betatildej,kp,mup,rp,rfp,yp,consp,xgridp,vp,app_min,evp)
jloop:do jc= nj-1,1,-1
        betatildej = beta*surv(jc)**(1.0/gamm)*(1.0+g)**((1.0-theta)/gamm)
!$OMP DO SCHEDULE(STATIC)
kloop:  do kc=1,nk          ! Small performance notice: the outermost loop does not correspond to the rightmost state, because I interchanged loops for k and mu, so that OpenMP can work on k.
muloop:     do muc=1,nmu
zloop:          do zc=1,nz

                    call calc_vars_tomorrow(coeffs,grids,jc,zc,kc,muc,kp,mup,rp,rfp,yp,err%kp (zc,kc,muc),err%mup(zc,kc,muc),err%rfp(zc,kc,muc))

                    call interp_policies_tomorrow(p, cons, value, kp, mup, grids, jc, consp, xgridp, vp, app_min)

etaloop:            do ec=1,n_eta
                        ! Create savings grid apgrid
                        p%apgrid(:,ec,zc,jc,kc,muc)= f_apgrid_j(rfp,yp, xgridp, app_min, apmax(ec,zc,jc), jc)

xloop:                  do xc=1,nx
                            call asset_allocation(xgridp, consp, vp, yp, rfp, rp, p%apgrid(xc,ec,zc,jc,kc,muc), pi_z(zc,:), pi_eta(ec,:), xc, jc, p%kappa(xc,ec,zc,jc,kc,muc), err%asset(xc,ec,zc,jc,kc,muc))
                            p%stocks(xc,ec,zc,jc,kc,muc) = p%apgrid(xc,ec,zc,jc,kc,muc) * p%kappa(xc,ec,zc,jc,kc,muc)

                            call consumption(p%apgrid(xc,ec,zc,jc,kc,muc), p%kappa(xc,ec,zc,jc,kc,muc), xgridp, consp, vp, rfp,rp, yp, zc, xc, ec, betatildej, cons(xc,ec,zc,jc,kc,muc), evp, err%cons(:,xc,ec,zc,jc,kc,muc))

                            value(xc,ec,zc,jc,kc,muc) = (cons(xc,ec,zc,jc,kc,muc)**((1.0-theta)/gamm) + betatildej*evp**(1.0/gamm))**(gamm/(1.0-theta))

                            ! create new grid for cash at hand (xgrid)
                            p%xgrid(xc,ec,zc,jc,kc,muc)=p%apgrid(xc,ec,zc,jc,kc,muc)+cons(xc,ec,zc,jc,kc,muc)

                        end do xloop
                    enddo etaloop
                enddo zloop
            enddo muloop
        enddo kloop
!$OMP END DO
    enddo jloop
!$OMP END PARALLEL

end subroutine olg_backwards_recursion
!-------------------------------------------------------------------------------

pure subroutine calc_vars_tomorrow(coeffs,grids,jc,zc,kc,muc,kp,mup,rp,rfp,yp, err_k, err_mu, err_rfp)
! Forecast kp, mup, and get corresonding prices  (might want to move into laws_of_motion)
    use params_mod     ,only: nz, pi_z, jr, ej, etagrid
    use laws_of_motion ,only: Forecast_k, Forecast_mu
    use income

    type(tCoeffs)   ,intent(in)  :: coeffs ! coefficients for laws of motion
    type(tAggGrids) ,intent(in)  :: grids  ! grids for aggregate states k and mu
    integer         ,intent(in)  :: jc, zc, kc, muc
    real(dp)        ,intent(out) :: kp, mup(:), rp(:), rfp, yp(:,:)
    logical(1)      ,intent(out) :: err_k, err_mu, err_rfp
    integer :: zpc, nk, nmu
    real(dp), parameter :: crit = 1.0e-10

    err_k   = .false.
    err_mu  = .false.
    err_rfp = .false.
    nk = size(grids%k)
    nmu= size(grids%mu)

    kp  = Forecast_k(coeffs%k(:,zc), grids%k(kc), grids%mu(muc))
    if (kp  > grids%k(nk)) then
        if (kp - grids%k(nk) > crit) err_k = .true.
        kp  = grids%k(nk)
    elseif (kp  < grids%k(1)) then
        if (grids%k(1) - kp  > crit) err_k = .true.
        kp  = grids%k(1)
    endif

    do zpc = 1,nz
        mup(zpc) = Forecast_mu(coeffs%mu(:,zpc), kp, grids%mu(muc))
    enddo
    if (any(mup - grids%mu(nmu) > crit) .or. any(grids%mu(1) - mup > crit)) err_mu = .true.
    where (mup > grids%mu(nmu)) mup = grids%mu(nmu) ! obsolete comment: This takes care of the wrong forecasts in the mean shock equilibrium
    where (mup < grids%mu(1))   mup = grids%mu(1)

    ! calculate tomorrow's risky returns and wage for given law of motion, for every zc today
    rfp = f_riskfree_rate(kp,grids%mu(muc),pi_z(zc,:))
    do zpc= 1,nz
        rp(zpc) = f_stock_return(kp, zeta(zpc), delta(zpc), rfp)
        if (jc+1>=jr) then
            yp(:,zpc) = f_pensions(kp, zeta(zpc))
        else
            yp(:,zpc) = ej(jc+1) * f_netwage(kp, zeta(zpc)) * etagrid(:,zpc)
        endif
    enddo
    if (rfp < rp(1)*(1.0 + sign(0.0001_dp,rp(1))) .and. grids%mu(muc)>0.0 ) then
        rfp = rp(1)*(1.0 + sign(0.0001_dp,rp(1)))
        err_rfp = .true.
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

pure function f_apgrid_j(rfp,yp, xgridp, app_min, apmax, jc)
! Create savings grid for generation j, apgrid(:,zc,jc,kc,muc)
    use params_mod , only: nap, collateral_constraint, cmin, g
    use makegrid_mod

    real(dp) ,dimension(nap) :: f_apgrid_j
    real(dp) ,intent(in)     :: rfp, yp(:,:), xgridp(:,:,:), app_min, apmax
    integer  ,intent(in)     :: jc
    real(dp)                 :: rtildemax_debt, apmin, xp_min

CC: if (collateral_constraint) then
        f_apgrid_j(1) = 0.0 ! a'=0, but I can still borrow and invest in stock
        f_apgrid_j(2:nap) = MakeGrid(0.0_dp, apmax,nap-1,1.5_dp)

    else    ! here comes the difficult part: need to determine maximum borrowing (natural borrowing limit)

        ! need to check whether zpc=1 or zpc=nz gives the smaller apmin (in absolute terms)
        rtildemax_debt  = (1.0+rfp)/(1.0+g)
        xp_min          = min(cmin, cmin + app_min, xgridp(1,1,1))
        apmin           = (xp_min-yp(1,1))/rtildemax_debt
        if (jc < 6) then
            apmin = xp_min*.8_dp
        endif

        f_apgrid_j= MakeGrid(apmin,apmax,nap,1.5_dp) !2.0 !,'chebyshev'

    endif CC

end function f_apgrid_j
!-------------------------------------------------------------------------------

pure subroutine asset_allocation(xgridp, consp, vp, yp, rfp, rp, ap, pi_zp, pi_etap, xc, jc, kappa_out, error)
    ! This subroutine will pass the asset euler equation as a function argument to a root finder.
    ! In this version, the function argument is an internal procedure, which is a thread-safe Fortran 2008 feature implemented
    ! in the Intel Fortran Compiler >= 11.0 and in gfortran >= 4.5
    use params_mod ,only: tol_asset_eul, opt_zbrak, kappa_in_01, nonnegative_stock, scale_AR, de_ratio, g
    use fun_zbrent     ! NR: Brent method
    use sub_zbrac      ! NR: outwards bracketing
    use sub_zbrak      ! NR: inwards bracketing

    real(dp) ,dimension(:,:,:) ,intent(in) :: xgridp, consp, vp
    real(dp)   ,intent(in)  :: ap, rfp, yp(:,:), rp(:), pi_zp(:), pi_etap(:)        ! apgrid(xc,ec,zc,jc)
    integer   ,intent(in)   :: xc, jc
    real(dp)   ,intent(out) :: kappa_out ! kappa(xc,ec,zc,jc)
    logical(1) ,intent(out) :: error
    real(dp)                :: kappa1, kappa2,kappa_brack1, kappa_brack2  ! temporary variables
    real(dp) ,dimension(1)  :: kappa_result
    real(dp) ,dimension(5)  :: f_test, kappa_test
    real(dp) ,dimension(:) ,allocatable :: kappal, kappau !sub_zbrak: lower and upper bounds of segment containing a root
    logical                 :: bracket_found

    error = .false.

    if (scale_AR == -1.0) then
        kappa_out = 1.0/(1.0 + de_ratio) ! So that market clearing obtains.
        return
    endif
    if (ap==0.0) then       ! this happens if (collateral_constraint), because then ap(1) == 0.0 and ap(2) == 0.0
        kappa_out = 0.0     ! it can also catch the (coincidental, unintended) case where apgrid contains zero
        return
    endif

    ! The following may lead to a lot of errs%cons, but that should not be a worry then
    if (xc == 1 .and. jc>=6 ) then  ! For jc<6, xc=1 does not correspond to (natural) borrowing limit, see f_apgrid_j above
        kappa_out = 0.0     ! At the (natural) borrowing limit, agents cannot invest in risky stock
        return
    endif

    if (kappa_in_01) then   ! since ap <0.0 is possible (if .not. collateral_constraint) this is not the same as no_short_selling
        if (nonnegative_stock .and. ap <=0.0) then
            kappa_out = 0.0
            return
        endif

        kappa1 = 0.0
        kappa2 = 1.0
        if (asseteuler_f(kappa1) <= 0.0) then
            kappa_out = kappa1
            return
        elseif(asseteuler_f(kappa2) >= 0.0 ) then
            kappa_out = kappa2
            return
        endif
    else
        ! Determine maximum leverage without risk of not paying back (i.e. dropping below xgrid(1,1,jc+1))
        ! only need to look at the worst state tomorrow
        kappa1=((xgridp(1,1,1)-yp(1,1))*(1.0+g)/ap-(1.0+rfp))/(rp(1)-rfp)
        ! It would be easier and just as justifiable to calculate the kappa1 that implies a return of -100%, i.e. zero assets left, i.e. Rtilde = -1.0
        ! This kappa1 would be somewhat lower, because in the above calculation, return can be even worse than 100% as long as agent receives some income.
        ! That is at the moment he has to pay up, there is no limited liability. However, during solution this never happens, but it can happen in simulations.

        kappa2=-1.0*sign(1.0,ap)
        !kappa2=((xgrid(1,1,jc+1)-yp(nz))*(1.0+g)/ap-(1.0+rfp(zc)))/(rp(zc,nz)-rfp(zc))
    endif

    ! Now find root
    bracket_found=.false.
    if (asseteuler_f(kappa1)*asseteuler_f(kappa2)<0.0) then
        ! Our maximum leverage bounds bracket a root
        kappa_brack1    = kappa1
        kappa_brack2    = kappa2
        bracket_found=.true.
    else
        if (opt_zbrak) then
            ! Look 'inwards', maybe multiple roots bracketed by maximum leverage bounds
            call s_zbrak(asseteuler_f,kappa1,kappa2,1000,kappal,kappau)
            if (size(kappal) > 0) then
                bracket_found=.true.
                kappa_brack1    = kappal(1)
                kappa_brack2    = kappau(1)
                deallocate(kappal,kappau)
            endif
        endif
        if (.not. bracket_found) then
            if (kappa_in_01) then
                kappa_brack1    = kappa2/3.0
                kappa_brack2    = kappa2*2.0/3.0
            else
                ! Look 'outwards', i.e. extend the bounds
                kappa_brack1 = -20.0_dp
                kappa_brack2 = 20.0_dp
                call s_zbrac(asseteuler_f,kappa_brack1,kappa_brack2,bracket_found)
            endif
        endif
    endif

    if (bracket_found) then
        ! call d_zbren(asseteuler_f,kappa_brack1,kappa_brack2,errabs=tol_asset_eul,errrel=tol_asset_eul)
        ! kappa(xc,zc,jc)=kappa2
        kappa_out= f_zbrent(asseteuler_f,kappa_brack1,kappa_brack2,tol_asset_eul)
    else
        kappa_test(1) = kappa1
        f_test(1) = asseteuler_f(kappa1)
        kappa_test(2) = kappa2
        f_test(2) = asseteuler_f(kappa2)
        kappa_test(3) = kappa_brack1
        f_test(3) = asseteuler_f(kappa_brack1)
        kappa_test(4) = kappa_brack2
        f_test(4) = asseteuler_f(kappa_brack2)
        !kappa_test(5) = p%kappa(xc,ec,zc,jc+1,kc,muc)
        !f_test(5) = asseteuler_f(p%kappa(xc,ec,zc,jc+1,kc,muc))

        kappa_result = minloc(abs(f_test))
        kappa_out= kappa_test(kappa_result(1))
        error=.true.
    endif

    if (nonnegative_stock .and. (ap * kappa_out < 0.0)) kappa_out = 0.0

contains
!---------------------------------------------------------------------------
! Euler equation
!---------------------------------------------------------------------------
    pure function asseteuler_f(kappa)
    use params_mod         ,only: theta, gamm, g, cmin
    use fun_lininterp

    real(dp)                  :: asseteuler_f
    real(dp) ,intent(in)      :: kappa
    real(dp) ,dimension(:) ,allocatable :: aeez      ! asset euler equation for each z
    real(dp) ,dimension(1)    :: cons_interp, vp_interp, xp, aeetemp
    real(dp) :: rtildep     ! rtilde prime, aprime
    integer                   :: zpc, epc, nz, n_eta

    nz = size(pi_zp)
    n_eta = size(pi_etap)
    allocate(aeez(nz))
    aeez = 0.0
    do zpc=1,nz
        rtildep     = (1.0+rfp+kappa*(rp(zpc)-rfp))/(1.0+g)
        do epc=1,n_eta
            xp          = yp(epc,zpc)+rtildep*ap
            cons_interp = f_lininterp(xgridp(:,epc,zpc), consp(:,epc,zpc), xp)
            vp_interp   = f_lininterp(xgridp(:,epc,zpc), vp(:,epc,zpc), xp)


            if (cons_interp(1) <=cmin) then
                cons_interp = cmin
                xp               = cmin + ap
                vp_interp   = f_lininterp(xgridp(:,epc,zpc),vp(:,epc,zpc),xp)
            endif
            if (vp_interp(1) <=cmin) then
                vp_interp(1)   = cmin
            endif
            aeetemp = vp_interp**((1.0-theta)*(gamm-1.0)/gamm)*cons_interp**((1.0-theta-gamm)/gamm)
!
!           if (xp(1) < xgridp(1,epc,zpc)) then
!               aeetemp = taylor_expansion(cons_interp, vp_interp,epc,zpc)  ! only works for theta\=1
!           else
!               aeetemp = vp_interp**((1.0-theta)*(gamm-1.0)/gamm)*cons_interp**((1.0-theta-gamm)/gamm)
!           endif

            aeez(zpc) = aeez(zpc) + pi_etap(epc)*aeetemp(1)
        enddo
    enddo

    asseteuler_f    = dot_product(pi_zp,aeez*(rp-rfp))

    end function asseteuler_f

    !-------------------------------------------------------------------------------
    ! Taylor expansion of Euler euqation (need to outcomment above)
    !-------------------------------------------------------------------------------
    pure function taylor_expansion(cons_in, vp_in,epc,zpc)
    ! only works for theta\=1
        use params_mod, only: theta, gamm
        real(dp), dimension(:), intent(in) :: cons_in, vp_in
        integer, intent(in)                :: epc, zpc
        real(dp), dimension(size(vp_in))   :: taylor_expansion

    taylor_expansion = vp(1,epc,zpc)**((1.0-theta)*(gamm-1.0)/gamm)*consp(1,epc,zpc)**((1.0-theta-gamm)/gamm) &
    ! First order Taylor expansion
        + (vp_in-vp(1,epc,zpc))*((1.0-theta)*(gamm-1.0)/gamm)*vp(1,epc,zpc)**((1.0-theta)*(gamm-1.0)/gamm -1.0) &
          *consp(1,epc,zpc)**((1.0-theta-gamm)/gamm) &
        + vp(1,epc,zpc)**((1.0-theta)*(gamm-1.0)/gamm) &
          *(cons_in-consp(1,epc,zpc))*((1.0-theta-gamm)/gamm)*consp(1,epc,zpc)**((1.0-theta-gamm)/gamm-1.0) &
    ! Second order Taylor expansion
        + 0.5_dp*(vp_in-vp(1,epc,zpc))**2.0*((1.0-theta)*(gamm-1.0)/gamm -1.0)*((1.0-theta)*(gamm-1.0)/gamm)*vp(1,epc,zpc)**((1.0-theta)*(gamm-1.0)/gamm-2.0) &
          *consp(1,epc,zpc)**((1.0-theta-gamm)/gamm) &
        + vp(1,epc,zpc)**((1.0-theta)*(gamm-1.0)/gamm) &
          *0.5_dp*(cons_in-consp(1,epc,zpc))**2.0*((1.0-theta-gamm)/gamm)*((1.0-theta-gamm)/gamm-1.0)*consp(1,epc,zpc)**((1.0-theta-gamm)/gamm-2.0) &
        + (vp_in-vp(1,epc,zpc))*((1.0-theta)*(gamm-1.0)/gamm)*vp(1,epc,zpc)**((1.0-theta)*(gamm-1.0)/gamm -1.0) &
          *(cons_in-consp(1,epc,zpc))*((1.0-theta-gamm)/gamm)*consp(1,epc,zpc)**((1.0-theta-gamm)/gamm-1.0)
    end function taylor_expansion

end subroutine asset_allocation
!-------------------------------------------------------------------------------

pure subroutine consumption(ap, kappa, xgridp, consp, vp, rfp,rp, yp, zc, xc, ec, betatildej, cons_out, evp, error)
    use fun_lininterp
    use params_mod, only : collateral_constraint, pi_z, pi_eta, nz, n_eta, g, cmin, theta, gamm

    real(dp)                   ,intent(in)  :: ap, kappa, rfp, betatildej, rp(:),yp(:,:), xgridp(:,:,:),consp(:,:,:), vp(:,:,:)
    integer                    ,intent(in)  :: zc, xc, ec
    real(dp)                   ,intent(out) :: cons_out, evp  ! cons(xc,ec,zc,jc), evp
    logical(1) ,dimension(nz)  ,intent(out) :: error
    real(dp)   ,dimension(n_eta)            :: cons_interp, vp_interp ! interpolated values of cons, vp
    real(dp)   ,dimension(nz)               :: evpz, rhs_temp ! temporary: cee_*: consumption euler eq.
    real(dp)                                :: rtildep, xp    ! rtilde and cash-at-hand tomorrow
    real(dp)                                :: cee_rhs, rhs_fac1, rhs_fac2 ! temporary: rhs of cee, factors 1 and 2
    integer                                 :: zpc, epc       ! counters for shocks tomorrow

    error = .false.

    do zpc=1,nz
        if (pi_z(zc,zpc) == 0.0) then
            evpz(zpc)     = 0.0
            rhs_temp(zpc) = 0.0
            cycle
        endif

        ! Get cash-at-hand tomorrow for all states
        rtildep= (1.0+rfp+kappa*(rp(zpc)-rfp))/(1.0+g)
        ! As mentioned in the subroutine asset_allocation above, one could limit kappa so that rtildep >= 0.0.
        ! In that case, this should be checked here also, because kappa can take different (very high) values during the simulations (when ap is close to zero).

        do epc = 1,n_eta
            xp = yp(epc,zpc)+rtildep*ap

            ! Interpolate consumption and value function
            cons_interp(epc) = f_lininterp(xgridp(:,epc,zpc),consp(:,epc,zpc),xp)
            vp_interp(epc)   = f_lininterp(xgridp(:,epc,zpc),vp(:,epc,zpc),xp)
            if (cons_interp(epc) <=cmin) then
                error(zpc)       = .true.
                cons_interp(epc) = cmin
                xp               = cmin + ap
                vp_interp(epc)   = f_lininterp(xgridp(:,epc,zpc),vp(:,epc,zpc),xp)
            endif
            if (vp_interp(epc) <=cmin) then
                error(zpc)       = .true.
                vp_interp(epc)   = cmin
            endif

        enddo
        ! Now different factors of consumption euler equation
        evpz(zpc)     = sum(pi_eta(ec,:)*vp_interp**(1.0-theta))
        rhs_temp(zpc) = sum(pi_eta(ec,:)*vp_interp**((1.0-theta)*&
                            (gamm-1.0)/gamm)*cons_interp**((1.0-theta-gamm)/gamm))
    enddo

    evp=dot_product(pi_z(zc,:),evpz)    ! Expected V'

    if (collateral_constraint .and. xc ==1) then
        cons_out = cmin/100.0
    else
        rhs_fac1=betatildej*evp**((1.0-gamm)/gamm)
        rhs_fac2=(1.0+rfp)/(1.0+g)*dot_product(pi_z(zc,:),rhs_temp)

        cee_rhs=rhs_fac1*rhs_fac2

        cons_out=max(cee_rhs**(gamm/(1.0-theta-gamm)), cmin)
    endif

end subroutine consumption
!-------------------------------------------------------------------------------

end module household_solution_mod
