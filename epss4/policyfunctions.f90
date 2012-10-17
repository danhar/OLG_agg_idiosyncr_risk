module policyfunctions
    implicit none
contains

subroutine calc_policyfunctions(coeffs, grids, p, value, err_o)
! Get the policy functions for the entire state space, i.e. both individual and aggregate states
    use kinds
    use params_mod      ,only: nj, nx, n_eta, nz, jr,surv, pi_z, pi_eta, cmin, g, beta, theta, gamm, apmax
    use types           ,only: tPolicies, AllocateType
    use aggregate_grids ,only: tAggGrids
    use laws_of_motion  ,only: tCoeffs
    use error_class
    use makegrid_mod
    use asseteuler     ,only: asseteuler_set

    type(tCoeffs)                    ,intent(in)  :: coeffs ! coefficients for laws of motion
    type(tAggGrids)                  ,intent(in)  :: grids  ! grids for aggregate states k and mu
    type(tPolicies)                  ,intent(out) :: p      ! policies
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
    call AllocateType(p,nk,nmu)
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
    cons(:,:,:,nj,:,:)    = p%xgrid(:,:,:,nj,:,:)
    value(:,:,:,nj,:,:)   = cons(:,:,:,nj,:,:)

    !---------------------------------------------------------------------------
    ! Model solution, generations nj-1 to 1
    !---------------------------------------------------------------------------
jloop:do jc= nj-1,1,-1
        betatildej = beta*surv(jc)*(1.0+g)**((1.0-theta)/gamm)

muloop: do muc=1,nmu
kloop:      do kc=1,nk
zloop: 			do zc=1,nz

                    call calc_vars_tomorrow(kp,mup,rp,rfp,yp, err_kp, err_mup, err_rfp)
                    if (present(err_o) .and. (err_kp .or. err_mup .or. err_rfp)) then
                        err_o%kp (zc,kc,muc) = err_kp
                        err_o%mup(zc,kc,muc) = err_mup
                        err_o%rfp(zc,kc,muc) = err_rfp
                    endif

                    call interp_policies(kp, mup, grids, p, consp, xgridp, vp, app_min)

etaloop:            do ec=1,n_eta

	                    ! Create savings grid apgrid
	                    p%apgrid(:,ec,zc,jc,kc,muc)= f_apgrid_j(rfp,yp, xgridp, app_min)

						call asseteuler_set(yp, rp, rfp, pi_z(zc,:), pi_eta(ec,:), consp, xgridp, vp)

xloop:                  do xc=1,nx

					        ! Set last euler variable
					        call asseteuler_set(p%apgrid(xc,ec,zc,jc,kc,muc))
					        ! Asset allocation problem, see internal subroutine below
							call asset_allocation(p%apgrid(xc,ec,zc,jc,kc,muc), p%kappa(xc,ec,zc,jc,kc,muc), err_asset)

							! Consumption problem, see internal subroutine below. Also returns evp.
							call consumption(p%apgrid(xc,ec,zc,jc,kc,muc), p%kappa(xc,ec,zc,jc,kc,muc), cons(xc,ec,zc,jc,kc,muc), evp, err_cons)

							! calculate new optimal value
							value(xc,ec,zc,jc,kc,muc) = (cons(xc,ec,zc,jc,kc,muc)**((1.0-theta)/gamm) + betatildej*evp**(1.0/gamm))**(gamm/(1.0-theta))

							! create new grid for cash at hand (xgrid)
							p%xgrid(xc,ec,zc,jc,kc,muc)=p%apgrid(xc,ec,zc,jc,kc,muc)+cons(xc,ec,zc,jc,kc,muc)

							if(present(err_o) .and. (err_asset .or. any(err_cons))) call error_handling(err_o, err_asset, err_cons, err_kp, err_mup)

	                    end do xloop
                    enddo etaloop
				enddo zloop
			enddo kloop
		enddo muloop
	enddo jloop

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - pure subroutine calc_vars_tomorrow(kp,mup,rp,rfp,yp)
! - pure function f_apgrid_j(kp, mup, rfp,yp, xgridp)
! - pure subroutine interp_policies(kp, mup, grid, p, consp, xgridp, vp)
! - pure subroutine asset_allocation(ap, kappa_out, error)
! - pure subroutine consumption(ap, kappa, cons_out, evp, error)
! - subroutine error_handling(err, err_asset, err_cons)
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

        kp  = Forecast(coeffs%k(:,zc), grids%k(kc))
        if (kp  > grids%k(nk)) then
            if (kp - grids%k(nk) > crit) err_k = .true.
            kp  = grids%k(nk)
        elseif (kp  < grids%k(1)) then
            if (grids%k(1) - kp  > crit) err_k = .true.
            kp  = grids%k(1)
        endif

        do zpc = 1,nz
            mup(zpc) = Forecast(coeffs%mu(:,zpc), kp, grids%mu(muc)) !grids%mu(muc) is hackish to distinguish for ms and STY
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

    !---------------------------------------------------------------------------
    ! Create savings grid for generation j, apgrid(:,zc,jc,kc,muc)
    !---------------------------------------------------------------------------
    pure function f_apgrid_j(rfp,yp, xgridp, app_min)
        use params_mod , only: nap, ap_numzero, collateral_constraint

        real(dp), dimension(nap) :: f_apgrid_j
        real(dp), intent(in)     :: rfp, yp(n_eta,nz), xgridp(nx,n_eta,nz), app_min
        real(dp)                 :: rtildemax_debt, apmin, xp_min
        integer                  :: nap1    ! number of of points below zero.

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

	        if (apmin <= -ap_numzero((jc+1)/jr+1)/2.0) then
	            nap1    = ceiling(abs(apmin)/(apmax(ec,zc,jc)-apmin)*real(nap,dp)) !+1
	            if (nap-nap1 <2) nap1=nap-2
	            nap1 =nap/2 !nap*2/3 !
	        else
	            nap1    = 0
	        endif

	        ! Create apgrid so that more points around ap_numzero
	        ! At the moment asymmetric around zero, but could make symmetric:
	        ! Closest point to zero from left is half the distance as the one from right
	        if (nap1==1) then
	            f_apgrid_j(1)=apmin !min(apmin, -ap_numzero((jc+1)/jr+1)))
	        elseif (nap1>1) then
	            f_apgrid_j(1:nap1)= MakeGrid(apmin,-ap_numzero((jc+1)/jr+1)/2.0,nap1,1.5_dp) !2.0 !,'chebyshev'
	        endif
	        if (nap-nap1<1) then
	            nap1 = nap1
	        endif
	        f_apgrid_j(nap1+1:nap)   =  MakeGrid(ap_numzero((jc+1)/jr+1), apmax(ec,zc,jc),nap-nap1,1.5_dp) !,'chebyshev'
	        !f_apgrid_j =  MakeGrid(apmin, apmax(jc),nap,1.0_dp) !,'chebyshev'
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

    end subroutine interp_policies

    !---------------------------------------------------------------------------
    ! Asset Allocation Problem
    !---------------------------------------------------------------------------
    pure subroutine asset_allocation(ap, kappa_out, error)
        use params_mod ,only: opt_zbren, tol_asset_eul, opt_zbrak, kappa_in_01, scale_AR, de_ratio
        use asseteuler ,only: asseteuler_f
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
        !real(dp), external     :: asseteuler_f ! function with asset euler equation, needs to be external for IMSL solvers

        error = .false.

        if (scale_AR == -1.0) then
            kappa_out = 1.0/(1.0 + de_ratio) ! So that market clearing obtains.
            return
        endif
        if (ap==0.0) then       ! this happens if (collateral_constraint), because then ap(1) == 0.0 and ap(2) == 0.0
            kappa_out = 0.0     ! it can also catch the (coincidental, unintended) case where apgrid contains zero (should not happen with ap_numzero>0)
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
                kappa_test(5) = p%kappa(xc,ec,zc,jc+1,kc,muc)
                f_test(5) = asseteuler_f(p%kappa(xc,ec,zc,jc+1,kc,muc))

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
                kappa_out=p%kappa(xc-1,ec,zc,jc,kc,muc)
            else
                kappa_out=kappa_result(1)
                !kappa(xc,zc,jc)=min(max(kappa_result(1),kappa1),kappa2)
                !kappa(xc,zc,jc)=min(max(kappa_result(1),-kappamax),kappamax)
            endif
        endif
    end subroutine asset_allocation

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
    subroutine error_handling(err, err_asset, err_cons, err_kp, err_mup)
        use params_mod ,only: detailed_euler_errs
        type(tErrors)  ,intent(inout) :: err
        logical(1)     ,intent(in)    :: err_asset, err_cons(:), err_kp, err_mup
        err%asset(xc,ec,zc,jc,kc,muc)  = err_asset
	    err%cons(:,xc,ec,zc,jc,kc,muc) = err_cons
	    if (detailed_euler_errs) then
            print '(a54,6i3)','ERROR: in policyfunctions at (xc,ec,zc,jc,kc,muc)=', xc,ec,zc,jc, kc, muc
            print '(t8,a9,l1,a12,l1,a14,l1,a14,<nz>(l1,x))', 'err_kp = ', err_kp, ' ,err_mup = ', err_mup, ' ,err_asset = ',err_asset,',  err_cons = ', err_cons
        endif
    end subroutine error_handling

!-------------------------------------------------------------------------------
! End of internal procedures and of main subroutine
!-------------------------------------------------------------------------------

end subroutine calc_policyfunctions
end module policyfunctions
