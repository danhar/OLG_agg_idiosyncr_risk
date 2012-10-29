module policyfunctions
    use kinds
    implicit none
    private

    public tPolicies

    type tPolicies
        real(dp), allocatable, dimension(:,:,:,:,:,:) :: apgrid, kappa, stocks, xgrid ! policies /grids
    contains
        procedure :: solve => solve_policyfunctions
        procedure :: allocate => allocate_policies
        procedure :: deallocate => deallocate_policies
        procedure :: calc_kappa
        procedure :: interpolate    ! At the moment, I am not using this anywhere
        procedure :: mean           ! At the moment, I am not using this anywhere
    end type tPolicies

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure subroutine solve_policyfunctions(coeffs, grids, p, value, err_o)
! - pure subroutine allocate_policies(p,nk,nmu)
! - pure subroutine deallocate_policies(p)
! - pure function mean(this,dimension_o,weight_o) result(mean_policy)
! - pure function interpolate(this,dim_x, gridx, x) result(pol_int)
! - pure subroutine calc_kappa(p)
!-------------------------------------------------------------------------------

pure subroutine solve_policyfunctions(p, coeffs, grids, value, err_o)
! Get the policy functions for the entire state space, i.e. both individual and aggregate states
    use params_mod      ,only: nj, nx, n_eta, nz, jr,surv, pi_z, pi_eta, cmin, g, beta, theta, gamm, apmax
    use aggregate_grids ,only: tAggGrids
    use laws_of_motion  ,only: tCoeffs
    use error_class
    use makegrid_mod

    class(tPolicies)                 ,intent(out) :: p      ! policies
    type(tCoeffs)                    ,intent(in)  :: coeffs ! coefficients for laws of motion
    type(tAggGrids)                  ,intent(in)  :: grids  ! grids for aggregate states k and mu
    real(dp)            ,allocatable ,intent(out) :: value(:,:,:,:,:,:)  ! could make optional
    type(tErrors)        ,optional   ,intent(out) :: err_o
    real(dp) ,dimension(:,:,:,:,:,:) ,allocatable :: cons
    real(dp) ,dimension(nx,n_eta,nz)              :: xgridp, consp, vp   ! xgrid, consumption, and value function tomorrow
    real(dp) ,dimension(n_eta,nz)                 :: yp      ! income tomorrow, for every idiosyncr & aggr state
    real(dp) ,dimension(nz)                       :: rp, mup ! risky return AND equity premium for every aggr state tomorrow
    real(dp)   :: kp, rfp, app_min                      ! tomorrow's capital, equity premium,risk-free rate, min aprime
    real(dp)   :: betatildej, evp                            ! modified discount factor, expected value tomorrow
    logical(1) :: err_kp, err_mup, err_rfp, err_cons(nz), err_asset
	integer    :: nk, nmu, jc, muc, kc, zc, xc, ec

    nk = size(grids%k)
    nmu= size(grids%mu)
    call p%allocate(nz,nk,nmu)
    allocate(cons(nx,n_eta,nz,nj,nk,nmu), value(nx,n_eta,nz,nj,nk,nmu))
    if (present(err_o)) call err_o%allocate(nk,nmu)

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
    associate(apgrid=> p%apgrid, kappa=> p%kappa, stocks => p%stocks, xgrid => p%xgrid)
    associate(gridk=> grids%k, gridmu=> grids%mu, coeffs_k => coeffs%k, coeffs_mu=> coeffs%mu)
    if (present(err_o)) associate(err_o_asset=> err_o%asset, err_o_cons=> err_o%cons, err_o_kp=> err_o%kp, err_o_mup=> err_o%mup, err_o_rfp=> err_o%rfp)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(apgrid,kappa,stocks,xgrid,gridk,gridmu,coeffs_k,coeffs_mu,value,cons,err_o_asset,err_o_cons,err_o_kp,err_o_mup,err_o_rfp,nmu,nk,nz,nj,n_eta,nx,beta,g,theta,gamm,surv) &
!$OMP PRIVATE(jc,muc,kc,zc,betatildej,kp,mup,rp,rfp,yp,err_kp,err_mup,err_rfp,consp,xgridp,vp,app_min,err_asset,evp,err_cons)
jloop:do jc= nj-1,1,-1
        betatildej = beta*surv(jc)*(1.0+g)**((1.0-theta)/gamm)
!$OMP DO SCHEDULE(STATIC)
muloop: do muc=1,nmu
kloop:      do kc=1,nk
zloop: 			do zc=1,nz

                    call calc_vars_tomorrow(kp,mup,rp,rfp,yp, err_kp, err_mup, err_rfp)
                    if (present(err_o_kp) .and. (err_kp .or. err_mup .or. err_rfp)) then
                        err_o_kp (zc,kc,muc) = err_kp
                        err_o_mup(zc,kc,muc) = err_mup
                        err_o_rfp(zc,kc,muc) = err_rfp
                    endif

                    call interp_policies(kp, mup, grids, p, consp, xgridp, vp, app_min)

etaloop:            do ec=1,n_eta
	                    ! Create savings grid apgrid
	                    apgrid(:,ec,zc,jc,kc,muc)= f_apgrid_j(rfp,yp, xgridp, app_min)

xloop:                  do xc=1,nx
					        ! Asset allocation problem, see internal subroutine below
							call asset_allocation(apgrid(xc,ec,zc,jc,kc,muc), kappa(xc,ec,zc,jc,kc,muc), err_asset)
							stocks(xc,ec,zc,jc,kc,muc) = apgrid(xc,ec,zc,jc,kc,muc) * kappa(xc,ec,zc,jc,kc,muc)

							! Consumption problem, see internal subroutine below. Also returns evp.
							call consumption(apgrid(xc,ec,zc,jc,kc,muc), kappa(xc,ec,zc,jc,kc,muc), cons(xc,ec,zc,jc,kc,muc), evp, err_cons)

							! calculate new optimal value
							value(xc,ec,zc,jc,kc,muc) = (cons(xc,ec,zc,jc,kc,muc)**((1.0-theta)/gamm) + betatildej*evp**(1.0/gamm))**(gamm/(1.0-theta))

							! create new grid for cash at hand (xgrid)
							xgrid(xc,ec,zc,jc,kc,muc)=apgrid(xc,ec,zc,jc,kc,muc)+cons(xc,ec,zc,jc,kc,muc)

							if(present(err_o) .and. (err_asset .or. any(err_cons))) call error_handling(err_o, err_asset, err_cons, err_kp, err_mup)

	                    end do xloop
                    enddo etaloop
				enddo zloop
			enddo kloop
		enddo muloop
!$OMP END DO
	enddo jloop
!$OMP END PARALLEL
if (present(err_o)) end associate
end associate
end associate

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - pure subroutine calc_vars_tomorrow(kp,mup,rp,rfp,yp)
! - pure function f_apgrid_j(kp, mup, rfp,yp, xgridp)
! - pure subroutine interp_policies(kp, mup, grid, p, consp, xgridp, vp)
! - pure subroutine asset_allocation(ap, kappa_out, error)
! - pure function asseteuler_f(kappa)
! - pure function taylor_expansion(cons_in, vp_in)
! - pure subroutine consumption(ap, kappa, cons_out, evp, error)
! - (pure) subroutine error_handling(err, err_asset, err_cons)
!-------------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! Forecast kp, mup, and get corresonding prices
    !---------------------------------------------------------------------------
    pure subroutine calc_vars_tomorrow(kp,mup,rp,rfp,yp, err_k, err_mu, err_rfp)
        use params_mod     ,only: ej, etagrid
        use laws_of_motion ,only: Forecast
        use income
        real(dp), intent(out) :: kp, mup(nz), rp(nz), rfp, yp(n_eta, nz)
        logical(1), intent(out) :: err_k, err_mu, err_rfp
        integer :: zpc
        real(dp), parameter :: crit = 1.0e-10

        err_k   = .false.
        err_mu  = .false.
        err_rfp = .false.

        kp  = Forecast(coeffs_k(:,zc), gridk(kc))
        if (kp  > gridk(nk)) then
            if (kp - gridk(nk) > crit) err_k = .true.
            kp  = gridk(nk)
        elseif (kp  < gridk(1)) then
            if (gridk(1) - kp  > crit) err_k = .true.
            kp  = gridk(1)
        endif

        do zpc = 1,nz
            mup(zpc) = Forecast(coeffs_mu(:,zpc), kp, gridmu(muc)) !gridmu(muc) is hackish to distinguish for ms and STY
        enddo
        if (any(mup - gridmu(nmu) > crit) .or. any(gridmu(1) - mup > crit)) err_mu = .true.
        where (mup > gridmu(nmu)) mup = gridmu(nmu) ! obsolete comment: This takes care of the wrong forecasts in the mean shock equilibrium
        where (mup < gridmu(1))   mup = gridmu(1)

        ! calculate tomorrow's risky returns and wage for given law of motion, for every zc today
	    rfp = f_riskfree_rate(kp,gridmu(muc),pi_z(zc,:))
	    do zpc= 1,nz
	        rp(zpc) = f_stock_return(kp, zeta(zpc), delta(zpc), rfp)
            if (jc+1>=jr) then
                yp(:,zpc) = f_pensions(kp, zeta(zpc))
            else
                yp(:,zpc) = ej(jc+1) * f_netwage(kp, zeta(zpc)) * etagrid(:,zpc)
            endif
	    enddo
        if (rfp < rp(1)*(1.0 + sign(0.0001_dp,rp(1))) .and. gridmu(muc)>0.0 ) then
            rfp = rp(1)*(1.0 + sign(0.0001_dp,rp(1)))
            err_rfp = .true.
        endif

    end subroutine calc_vars_tomorrow

    !---------------------------------------------------------------------------
    ! Create savings grid for generation j, apgrid(:,zc,jc,kc,muc)
    !---------------------------------------------------------------------------
    pure function f_apgrid_j(rfp,yp, xgridp, app_min)
        use params_mod , only: nap, collateral_constraint

        real(dp), dimension(nap) :: f_apgrid_j
        real(dp), intent(in)     :: rfp, yp(n_eta,nz), xgridp(nx,n_eta,nz), app_min
        real(dp)                 :: rtildemax_debt, apmin, xp_min

CC:     if (collateral_constraint) then
            f_apgrid_j(1) = 0.0 ! a'=0, but I can still borrow and invest in stock
            f_apgrid_j(2:nap) = MakeGrid(0.0_dp, apmax(ec,zc,jc),nap-1,1.5_dp)

        else    ! here comes the difficult part: need to determine maximum borrowing (natural borrowing limit)

	        ! need to check whether zpc=1 or zpc=nz gives the smaller apmin (in absolute terms)
	        rtildemax_debt  = (1.0+rfp)/(1.0+g)
	        xp_min          = min(cmin, cmin + app_min, xgridp(1,1,1))
	        apmin           = (xp_min-yp(1,1))/rtildemax_debt
	        if (jc < 6) then
	            apmin = xp_min*.8_dp
	        endif

            f_apgrid_j= MakeGrid(apmin,apmax(ec,zc,jc),nap,1.5_dp) !2.0 !,'chebyshev'

        endif CC

    end function f_apgrid_j

    !---------------------------------------------------------------------------
    ! Make a projection of tomorrow's policies on k prime and mu prime
    !---------------------------------------------------------------------------
    pure subroutine interp_policies(kp, mup, grid, p, consp, xgridp, vp, app_min)
        use fun_locate

	    real(dp), intent(in) :: kp, mup(nz)
	    type(tAggGrids), intent(in) :: grid
	    type(tPolicies), intent(in) :: p
	    real(dp), dimension(nx,n_eta,nz) ,intent(out) :: consp, xgridp, vp
	    real(dp)                         ,intent(out) :: app_min
	    integer :: iK, imu, zpc
	    real(dp) :: wK, wmu

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

		        xgridp(:,:,zpc)= (1-wK)*(1-wmu)*xgrid(:,:,zpc,jc+1,iK  ,imu  ) + &
				                    wK *(1-wmu)*xgrid(:,:,zpc,jc+1,iK+1,imu  ) + &
				                 (1-wK)*   wmu *xgrid(:,:,zpc,jc+1,iK  ,imu+1) + &
				                    wK *   wmu *xgrid(:,:,zpc,jc+1,iK+1,imu+1)

		        vp(:,:,zpc)    = (1-wK)*(1-wmu)*  value(:,:,zpc,jc+1,iK  ,imu  ) + &
				                    wK *(1-wmu)*  value(:,:,zpc,jc+1,iK+1,imu  ) + &
				                 (1-wK)*   wmu *  value(:,:,zpc,jc+1,iK  ,imu+1) + &
				                    wK *   wmu *  value(:,:,zpc,jc+1,iK+1,imu+1)
            enddo

            imu        = f_locate(grid%mu, mup(1))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
            wmu        = (mup(1) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))
            app_min= (1-wK)*(1-wmu)*apgrid(1,1,1,jc+1,iK  ,imu  ) + & ! This creates smallest aprime for forecasts
                        wK *(1-wmu)*apgrid(1,1,1,jc+1,iK+1,imu  ) + & ! I should move this into the loop and calc for all zpc
                     (1-wK)*   wmu *apgrid(1,1,1,jc+1,iK  ,imu+1) + & ! coz then better to understand.
                        wK *   wmu *apgrid(1,1,1,jc+1,iK+1,imu+1)

        elseif (size(grid%K)>1) then
            ! Projection of policies on Kt
            iK        = f_locate(grid%K, kp)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
            wK        = (kp - grid%k(iK)) / (grid%k(iK+1) - grid%k(iK))
            consp  = (1.0 -wK)*cons    (:,:,:,jc+1,iK,1) + wK*cons    (:,:,:,jc+1,iK+1,1)
            xgridp = (1.0 -wK)*xgrid (:,:,:,jc+1,iK,1) + wK*xgrid (:,:,:,jc+1,iK+1,1)
            vp     = (1.0 -wK)*value   (:,:,:,jc+1,iK,1) + wK*value   (:,:,:,jc+1,iK+1,1)
            app_min= (1.0 -wK)*apgrid(1,1,1,jc+1,iK,1) + wK*apgrid(1,1,1,jc+1,iK+1,1)

        elseif (size(grid%mu)>1) then
            do zpc=1,nz
                imu        = f_locate(grid%mu, mup(zpc))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
                wmu        = (mup(zpc) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))

                ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds
                consp(:,:,zpc) = (1.0-wmu)*   cons(:,:,zpc,jc+1,1,imu) + wmu*   cons(:,:,zpc,jc+1,1,imu+1)
                xgridp(:,:,zpc)= (1.0-wmu)*xgrid(:,:,zpc,jc+1,1,imu) + wmu*xgrid(:,:,zpc,jc+1,1,imu+1)
                vp(:,:,zpc)    = (1.0-wmu)*  value(:,:,zpc,jc+1,1,imu) + wmu*  value(:,:,zpc,jc+1,1,imu+1)
            enddo

            imu        = f_locate(grid%mu, mup(1))   ! In 'default', returns iu-1 if x>xgrid(iu-1)
            wmu        = (mup(1) - grid%mu(imu)) / (grid%mu(imu+1) - grid%mu(imu))
            app_min= (1.0-wmu)*apgrid(1,1,1,jc+1,1,imu) + wmu*apgrid(1,1,1,jc+1,1,imu+1)
            ! This creates smallest aprime for forecasts. should move this into the loop and calc for all zpc coz then better to understand.

        else    ! Mean shock
            consp  = cons    (:,:,:,jc+1,1,1)
            xgridp = xgrid (:,:,:,jc+1,1,1)
            vp     = value   (:,:,:,jc+1,1,1)
            app_min= apgrid(1,1,1,jc+1,1,1)
        endif

    end subroutine interp_policies

    !---------------------------------------------------------------------------
    ! Asset Allocation Problem
    !---------------------------------------------------------------------------
    pure subroutine asset_allocation(ap, kappa_out, error)
        use params_mod ,only: opt_zbren, tol_asset_eul, opt_zbrak, kappa_in_01, scale_AR, de_ratio
!       use zreal_int      ! IMSL Math.pdf, p. 1195f: Mullers Method to find roots (like secant but quadratic).
                           ! expects array as input (i.e. define scalar as dimension(1)).
!       use zbren_int      ! IMSL Math.pdf, p. 1192f: Brent's method to find roots
        use fun_zbrent     ! NR: Brent method
        use sub_zbrac      ! NR: outwards bracketing
        use sub_zbrak      ! NR: inwards bracketing

        real(dp)   ,intent(in)  :: ap        ! apgrid(xc,ec,zc,jc)
        real(dp)   ,intent(out) :: kappa_out ! kappa(xc,ec,zc,jc)
        logical(1) ,intent(out) :: error
        real(dp)                :: kappa1, kappa2,kappa_brack1, kappa_brack2  ! temporary variables
        real(dp) ,dimension(1)  :: kappa_start, kappa_result ! zreal
        real(dp) ,dimension(5)  :: f_test, kappa_test
        integer  ,dimension(1)  :: zreal_its ! zreal: no of iterations to convergence. need to test convergence
        real(dp) ,dimension(:) ,allocatable :: kappal, kappau !sub_zbrak: lower and upper bounds of segment containing a root
        logical                 :: bracket_found
        !real(dp), external     :: asseteuler_f ! function with asset euler equation, needs to be external for IMSL solvers (or can it be internal?)

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
	        kappa2=-1.0*sign(1.0,ap)
	        !kappa2=((xgrid(1,1,jc+1)-yp(nz))*(1.0+g)/ap-(1.0+rfp(zc)))/(rp(zc,nz)-rfp(zc))
        endif

        ! Now find root, depending on root finder
        if (opt_zbren) then ! Brent's Method
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
                if (.not. bracket_found .and. .not. kappa_in_01) then
                    ! Look 'outwards', i.e. extend the bounds
                    kappa_brack1 = -20.0_dp
                    kappa_brack2 = 20.0_dp
                    call s_zbrac(asseteuler_f,kappa_brack1,kappa_brack2,bracket_found)
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
                kappa_test(5) = kappa(xc,ec,zc,jc+1,kc,muc)
                f_test(5) = asseteuler_f(kappa(xc,ec,zc,jc+1,kc,muc))

                kappa_result = minloc(abs(f_test))
                kappa_out= kappa_test(kappa_result(1))
                error=.true.
            endif

        else    ! Mueller's Method
            kappa_start=(kappa1+kappa2)/2
            ! maybe I should scale as mentioned in the IMSL documentation, since some kappa<1
            !call d_zreal(asseteuler_f,kappa_result,xguess=kappa_start,itmax=1000,info=zreal_its,errabs=tol_asset_eul,errrel=tol_asset_eul)
            if (zreal_its(1)>1000) then
                error=.true.
                kappa_out=kappa(xc-1,ec,zc,jc,kc,muc)
            else
                kappa_out=kappa_result(1)
                !kappa(xc,zc,jc)=min(max(kappa_result(1),kappa1),kappa2)
                !kappa(xc,zc,jc)=min(max(kappa_result(1),-kappamax),kappamax)
            endif
        endif
    end subroutine asset_allocation

    !---------------------------------------------------------------------------
    ! Euler equation
    !---------------------------------------------------------------------------

	pure function asseteuler_f(kappa)
    use params_mod         ,only: theta, gamm, g, cmin
    use fun_lininterp

    real(dp)                  :: asseteuler_f
    real(dp) ,intent(in)      :: kappa
    real(dp) :: aeez(nz)                 ! asset euler equation for each z
    real(dp) ,dimension(1)    :: cons_interp, vp_interp, xp, aeetemp
    real(dp) :: rtildep, ap     ! rtilde prime, aprime
    integer                   :: zpc, epc

    ap = apgrid(xc,ec,zc,jc,kc,muc)
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

            aeez(zpc) = aeez(zpc) + pi_eta(ec,epc)*aeetemp(1)
        enddo
    enddo

    asseteuler_f    = dot_product(pi_z(zc,:),aeez*(rp-rfp))

    end function asseteuler_f

    !-------------------------------------------------------------------------------
    ! Taylor expansion of Euler euqation (optional)
    !-------------------------------------------------------------------------------

    pure function taylor_expansion(cons_in, vp_in,epc,zpc)
    ! only works for theta\=1
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


    !---------------------------------------------------------------------------
    ! Consumption problem
    !---------------------------------------------------------------------------
    pure subroutine consumption(ap, kappa, cons_out, evp, error)
        use fun_lininterp
        use params_mod, only : collateral_constraint

        real(dp)                   ,intent(in)  :: ap, kappa      ! apgrid(xc,ec,zc,jc), kappa(xc,ec,zc,jc)
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
                                (gamm-1.0)/gamm)*cons_interp**((1.0-theta-gamm)/gamm)*rtildep)
        enddo

        evp=dot_product(pi_z(zc,:),evpz)    ! Expected V'

        if (collateral_constraint .and. xc ==1) then
            cons_out = cmin/100.0
        else
	        rhs_fac1=betatildej*evp**((1.0-gamm)/gamm)
	        rhs_fac2=dot_product(pi_z(zc,:),rhs_temp)
	        cee_rhs=rhs_fac1*rhs_fac2

	        cons_out=max(cee_rhs**(gamm/(1.0-theta-gamm)), cmin)
        endif

    end subroutine consumption

    !---------------------------------------------------------------------------
    ! Print errors
    !---------------------------------------------------------------------------
    pure subroutine error_handling(err, err_asset, err_cons, err_kp, err_mup)
        use params_mod ,only: detailed_euler_errs
        type(tErrors)  ,intent(inout) :: err
        logical(1)     ,intent(in)    :: err_asset, err_cons(:), err_kp, err_mup
        err%asset(xc,ec,zc,jc,kc,muc)  = err_asset
	    err%cons(:,xc,ec,zc,jc,kc,muc) = err_cons
	    if (detailed_euler_errs) then ! the following is only for serious debugging. Outcommenting destroys pure!
!            print '(a54,6i3)','ERROR: in policyfunctions at (xc,ec,zc,jc,kc,muc)=', xc,ec,zc,jc, kc, muc
!            print '(t8,a9,l1,a12,l1,a14,l1,a14,<nz>(l1,x))', 'err_kp = ', err_kp, ' ,err_mup = ', err_mup, ' ,err_asset = ',err_asset,',  err_cons = ', err_cons
        endif
    end subroutine error_handling
end subroutine solve_policyfunctions
!-------------------------------------------------------------------------------

pure subroutine allocate_policies(this,nz,nk,nmu)
    use params_mod, only: nx, n_eta, nj
    class(tPolicies), intent(inout)  :: this
    integer,    intent(in)      :: nz,nk,nmu
    call deallocate_policies(this)
    allocate(this%apgrid(nx,n_eta,nz,nj,nk,nmu),this%kappa(nx,n_eta,nz,nj,nk,nmu),this%stocks(nx,n_eta,nz,nj,nk,nmu),this%xgrid(nx,n_eta,nz,nj,nk,nmu))
end subroutine allocate_policies
!-------------------------------------------------------------------------------

pure subroutine deallocate_policies(this)
    class(tPolicies), intent(inout)  :: this
    ! deallocating in reverse order to allocation for memory purposes
    if (allocated(this%xgrid)) deallocate(this%xgrid)
    if (allocated(this%stocks)) deallocate(this%stocks)
    if (allocated(this%kappa)) deallocate(this%kappa)
    if (allocated(this%apgrid)) deallocate(this%apgrid)
end subroutine deallocate_policies
!-------------------------------------------------------------------------------

pure function mean(this,dimension_o,weight_o) result(mean_policy)
    ! At the moment, I am not using this anywhere, because typically I want the result to have one dimension less, which is not the case here
    class(tPolicies), intent(in)  :: this
    type(tPolicies)               :: mean_policy
    integer, intent(in), optional:: dimension_o
    real(dp), intent(in), optional:: weight_o
    integer :: dimension, nd, dc
    real(dp) :: weight

    if (.not. present(dimension_o)) then
        dimension = 3 ! corresponds to aggregate shocks nz
    else
        dimension = dimension_o
    endif

    nd = size(this%apgrid,dimension)

    if (.not. present(weight_o)) then
        weight = 1.0/nd
    else
        weight = weight_o
    endif

    select case (dimension)
    case (3)
        call mean_policy%allocate(1,size(this%apgrid,5),size(this%apgrid,6)) ! policies for given z, K, and mu
        mean_policy%apgrid = 0.0
        mean_policy%stocks = 0.0
        mean_policy%xgrid  = 0.0
        do dc = 1,nd
            mean_policy%apgrid(:,:,1,:,:,:) = mean_policy%apgrid(:,:,1,:,:,:) + weight*this%apgrid(:,:,dc,:,:,:)
            mean_policy%stocks(:,:,1,:,:,:) = mean_policy%stocks(:,:,1,:,:,:) + weight*this%stocks(:,:,dc,:,:,:)
            mean_policy%xgrid (:,:,1,:,:,:) = mean_policy%xgrid (:,:,1,:,:,:) + weight*this%xgrid (:,:,dc,:,:,:)
        enddo
    case default ! same as case 3
        call mean_policy%allocate(1,size(this%apgrid,5),size(this%apgrid,6)) ! policies for given z, K, and mu
        mean_policy%apgrid = 0.0
        mean_policy%stocks = 0.0
        mean_policy%xgrid  = 0.0
        do dc = 1,nd
            mean_policy%apgrid(:,:,1,:,:,:) = mean_policy%apgrid(:,:,1,:,:,:) + weight*this%apgrid(:,:,dc,:,:,:)
            mean_policy%stocks(:,:,1,:,:,:) = mean_policy%stocks(:,:,1,:,:,:) + weight*this%stocks(:,:,dc,:,:,:)
            mean_policy%xgrid (:,:,1,:,:,:) = mean_policy%xgrid (:,:,1,:,:,:)  + weight*this%xgrid(:,:,dc,:,:,:)
        enddo
    end select

    call mean_policy%calc_kappa
end function mean
!-------------------------------------------------------------------------------

pure function interpolate(this,dim_x, gridx, x) result(pol_int)
    ! linear interpolation for one value in one dimension
    ! At the moment, I am not using this anywhere, because typically I want the result to have one dimension less, which is not the case here
    ! For this, would need another type tPolicies5 (or Policies_d5) and
    ! class(tPolicies),allocatable  :: pol_int ! allocate later
    ! or
    ! type(Policies_d5)  :: pol_int
    use fun_locate      ,only: f_locate

    class(tPolicies), intent(in)  :: this
    class(tPolicies),allocatable  :: pol_int
    integer, intent(in)           :: dim_x
    real(dp), intent(in) :: x, gridx(:)
    integer :: i
    real(dp) :: w

    i        = f_locate(gridx, x)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
    w        = (x - gridx(i)) / (gridx(i+1) - gridx(i))
    ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds

    select case (dim_x)
    case (5)
        call pol_int%allocate(size(this%apgrid,3),1,size(this%apgrid,6)) ! policies for given z, K, and mu
        pol_int%apgrid(:,:,:,:,1,:)= (1-w)*this%apgrid(:,:,:,:,i,:) +w*this%apgrid(:,:,:,:,i+1,:)
        pol_int%stocks(:,:,:,:,1,:)= (1-w)*this%stocks(:,:,:,:,i,:) +w*this%stocks(:,:,:,:,i+1,:)
        pol_int%xgrid (:,:,:,:,1,:) = (1-w)*this%xgrid (:,:,:,:,i,:) +w*this%xgrid (:,:,:,:,i+1,:)
    case default ! same as case 5
        call pol_int%allocate(size(this%apgrid,3),1,size(this%apgrid,6)) ! policies for given z, K, and mu
        pol_int%apgrid(:,:,:,:,1,:)= (1-w)*this%apgrid(:,:,:,:,i,:) +w*this%apgrid(:,:,:,:,i+1,:)
        pol_int%stocks(:,:,:,:,1,:)= (1-w)*this%stocks(:,:,:,:,i,:) +w*this%stocks(:,:,:,:,i+1,:)
        pol_int%xgrid (:,:,:,:,1,:) = (1-w)*this%xgrid (:,:,:,:,i,:) +w*this%xgrid (:,:,:,:,i+1,:)
    end select

    call pol_int%calc_kappa
end function interpolate
!-------------------------------------------------------------------------------

pure subroutine calc_kappa(this)
    class(tPolicies), intent(inout)  :: this
    where (this%apgrid .ne. 0.0)
        this%kappa = this%stocks/this%apgrid
    elsewhere
        this%kappa = 0.0
    end where
end subroutine calc_kappa
!-------------------------------------------------------------------------------

end module policyfunctions
