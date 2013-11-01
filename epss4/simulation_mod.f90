module simulation_mod
    use kinds           ,only: dp
    use classes_mod     ,only: tSimvars, tLifecycle, tPolicies, tAggGrids, tCoeffs

    implicit none
    private
    public simulate, calc_inequality_measures, f_euler_errors
contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure subroutine simulate(policies, value, agg_grid, simvars, Phi, lc)
! - pure subroutine calc_inequality_measures(simvars, xgridt, apgridt, stockst, Phi, etagridt, penst, netwaget, tc)
! - pure real(dp) function f_euler_errors()
!-------------------------------------------------------------------------------

pure subroutine simulate(policies, value, agg_grid, coeffs, calc_euler_errors, simvars, Phi, lc)
! Performs the Krusell-Smith simulation step and records lifecycle statistics
    use params_mod      ,only: n,g,L_N_ratio,pi_z,etagrid,t_scrap,exogenous_xgrid, partial_equilibrium, zeta, delta, alpha, tol_mut=> tol_simulation_marketclearing
    use income          ,only: f_netwage, f_pensions, f_stock_return, f_riskfree_rate, f_tau
    use fun_locate      ,only: f_locate
    use distribution    ,only: TransitionPhi, CheckPhi
    use fun_zbrent
    use fun_aggregate_diff
    use sorting_partial_mod ,only : valnth     ! for Median

    type(tPolicies)  ,intent(in)    :: policies
    real(dp)         ,intent(in)    :: value(:,:,:,:,:,:)
    type(tAggGrids)  ,intent(in)    :: agg_grid
    type(tCoeffs)    ,intent(in)    :: coeffs
    logical          ,intent(in)    :: calc_euler_errors
    type(tSimvars)   ,intent(inout) :: simvars    ! (zt, kt, mut, bt,...), first element contains starting values
    real(dp)         ,intent(inout) :: Phi(:,:,:) ! distribution. Returns: average Phi if (exogenous_xgrid), else Phi in nt
    type(tLifecycle) ,intent(out)   :: lc         ! lifecycle profiles
    real(dp) ,dimension(:,:,:,:) ,allocatable :: apgrid_zk, kappa_zk, xgrid_zk, stocks_zk ! policies for given z and K
    real(dp) ,dimension(:,:,:)   ,allocatable :: apgridt, kappat, xgridt, stockst         ! policies for given z, K, and mu
    real(dp) ,dimension(:,:,:)   ,allocatable :: Phi_avg, r_pf, val_j1_zk ! portfolio return
    real(dp) ,dimension(:,:)     ,allocatable :: val_j1_t ! value of j=1 for given z and K, and mu
    real(dp) ,dimension(:)       ,allocatable :: ap_lct, stocks_lct, cons_lct, cons_var_lct, return_lct, return_var_lct
    real(dp) ,dimension(:)       ,allocatable :: Knew       ! partial equilibrium: save aggregate stock in t
    real(dp)  :: Kt, mut, rt, netwaget, penst, w, eul_err_temp(2) ! variables in period t
    integer   :: tc, i, zt, jc, nmu, nx, n_eta, nj, nt, nk

    ! Intel Fortran Compiler XE 13.0 Update 1 (and previous) has a bug on realloc on assignment. If that is corrected, I think I can remove this whole allocation block
    nmu = size(agg_grid%mu); nk= size(agg_grid%k); nx=size(value,1); n_eta=size(value,2); nj=size(value,4); nt=size(simvars%z)
    allocate(apgrid_zk(nx,n_eta,nj,nmu), kappa_zk(nx,n_eta,nj,nmu), xgrid_zk(nx,n_eta,nj,nmu), stocks_zk(nx,n_eta,nj,nmu))
    allocate(apgridt(nx,n_eta,nj), kappat(nx,n_eta,nj), xgridt(nx,n_eta,nj), stockst(nx,n_eta,nj), Phi_avg(nx,n_eta,nj), r_pf(nx,n_eta,nj))
    allocate(val_j1_zk(nx,n_eta,nmu), val_j1_t(nx,n_eta))
    allocate(ap_lct(nj), stocks_lct(nj), cons_lct(nj), cons_var_lct(nj), return_lct(nj), return_var_lct(nj))
    allocate(Knew(nt+1))

    Phi_avg         = 0.0
    call lc%set_number(0.0_dp)
    simvars%err_K   = .false.
    simvars%err_mu  = .false.
    Knew(1)   = simvars%K(1)    ! This is only interesting for PE

    do tc=1,nt
        zt = simvars%z(tc)
        Kt = simvars%K(tc)
        if (Kt <= agg_grid%K(1) .or. Kt >= agg_grid%K(nk)) simvars%err_K(tc) = .true.

        ! Prices
        netwaget = f_netwage(Kt, zeta(zt))
        penst    = f_pensions(Kt, zeta(zt))
        rt       = f_stock_return(Kt, zeta(zt), delta(zt), simvars%rf(tc))

        simvars%bequests(tc) = f_bequests(simvars%rf(tc), rt, stockst, apgridt, Phi)

        ! Projection of policies on Kt
        i        = f_locate(agg_grid%K, Kt)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        w        = (Kt - agg_grid%K(i)) / (agg_grid%K(i+1) - agg_grid%K(i))
        ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds
        !! Once the Intel Fortran bug on reallocation is fixed, I can simply write xgrid_zk= ...
        xgrid_zk (:,:,:,:)= (1-w)*policies%xgrid (:,:,zt,:,i,:) +w*policies%xgrid (:,:,zt,:,i+1,:)
        apgrid_zk(:,:,:,:)= (1-w)*policies%apgrid(:,:,zt,:,i,:) +w*policies%apgrid(:,:,zt,:,i+1,:)
        stocks_zk(:,:,:,:)= (1-w)*policies%stocks(:,:,zt,:,i,:) +w*policies%stocks(:,:,zt,:,i+1,:)
        val_j1_zk(:,:,:)  = (1-w)*value(:,:,zt,1,i,:) +w*value(:,:,zt,1,i+1,:)
        where (apgrid_zk .ne. 0.0)
            kappa_zk = stocks_zk/apgrid_zk
        elsewhere
            kappa_zk = 0.0
        end where

tc1:    if (tc == 1) then
            if (nmu > 1) then
                i        = f_locate(agg_grid%mu, simvars%mu(1))
                w        = (simvars%mu(1) - agg_grid%mu(i)) / (agg_grid%mu(i+1) - agg_grid%mu(i))
                apgridt  = (1-w)*apgrid_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)
                stockst  = (1-w)*stocks_zk(:,:,:,i) + w*stocks_zk(:,:,:,i+1)
            else
                apgridt  = apgrid_zk(:,:,:,1)
                stockst  = stocks_zk(:,:,:,1)
            endif
        endif tc1

ex:     if (exogenous_xgrid .or. (nmu ==1) ) then
            xgridt   = xgrid_zk(:,:,:,1)
            ! To calc distribution, we need xgridt , along with netwaget etc, of this period (tc), but policy projections of last
            Phi      = TransitionPhi(simvars%rf(tc),rt,netwaget,penst,simvars%bequests(tc),xgridt,apgridt,stockst,etagrid(:,zt), Phi)
        endif ex

mu:     if (partial_equilibrium) then
            mut      = simvars%mu(tc)
            if ((agg_grid%mu(1)-mut)>tol_mut .or. (mut-agg_grid%mu(nmu))>tol_mut) simvars%err_mu(tc) = .true.

        else   ! Find mut that clears bondmarket
            if (f_excessbonds(agg_grid%mu(1)) * f_excessbonds(agg_grid%mu(nmu)) < 0.0 ) then
                mut  = f_zbrent(f_excessbonds,agg_grid%mu(1),agg_grid%mu(nmu),tol_mut)
            else
                simvars%err_mu(tc) = .true.
                if (abs(f_excessbonds(agg_grid%mu(1))) < abs(f_excessbonds(agg_grid%mu(nmu)))) then
                    mut = agg_grid%mu(1)
                else
                    mut = agg_grid%mu(nmu)
                endif
            endif
        endif mu

        ! Projection of policies on mut
        if (nmu > 1) then
	        i        = f_locate(agg_grid%mu, mut)
	        w        = (mut - agg_grid%mu(i)) / (agg_grid%mu(i+1) - agg_grid%mu(i))
	        if (.not. exogenous_xgrid) then
		        xgridt   = (1-w)* xgrid_zk(:,:,:,i) + w* xgrid_zk(:,:,:,i+1)
		        ! To calc distribution, we need xgridt , along with netwaget etc, of this period (tc), but policy projections of last
		        Phi      = TransitionPhi(simvars%rf(tc),rt,netwaget,penst,simvars%bequests(tc),xgridt,apgridt,stockst,etagrid(:,zt), Phi)
	        endif

	        apgridt  = (1-w)*apgrid_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)
	        stockst  = (1-w)*stocks_zk(:,:,:,i) + w*stocks_zk(:,:,:,i+1)
	        val_j1_t = (1-w)*val_j1_zk(:,:,i)   + w*val_j1_zk(:,:,i+1)
	        where (apgridt .ne. 0.0)
	            kappat = stockst/ apgridt
	        elsewhere ! includes apgrid(:,nj) =0
	            kappat = 0.0
	        end where
        else
            apgridt  = apgrid_zk(:,:,:,1)
            stockst  = stocks_zk(:,:,:,1)
            val_j1_t = val_j1_zk(:,:,1)
            kappat   = kappa_zk (:,:,:,1)
        endif

        call CheckPhi(Phi, simvars%Phi_1(tc), simvars%Phi_nx(tc)) ! Do here because Phi now was transitioned for all cases

        ! Aggregate
        Knew(tc+1) = sum(apgridt*Phi) /(L_N_ratio*(1.0+n)*(1.0+g))  ! in units of efficient labor
        ! Could allow some deviation out of grid - recorded as error in next period. Also in PE to keep comparable!
        if (Knew(tc+1)     < agg_grid%K(1 )) then
            Knew(tc+1)     = agg_grid%K(1 )
        elseif (Knew(tc+1) > agg_grid%K(nk)) then
            Knew(tc+1)     = agg_grid%K(nk)
        endif

        if (.not. partial_equilibrium) simvars%K(tc+1) = Knew(tc+1)

        simvars%output(tc)= zeta(zt)*simvars%K(tc)**alpha
        simvars%stock(tc) = sum(stockst*Phi)/L_N_ratio ! different from K(t+1) since it is in today's per capita terms! Thus no (1+g)(1+n) in denominator.
        simvars%bonds(tc) = sum((apgridt-stockst)*Phi)/L_N_ratio
        simvars%invest(tc)   = simvars%stock(tc) + simvars%bonds(tc) -simvars%K(tc)*(1.0-delta(zt))
        simvars%C(tc)     = sum((xgridt-apgridt)*Phi) / L_N_ratio ! in units of efficient labor
        simvars%r(tc)     = rt
        simvars%wage(tc)  = netwaget
        simvars%pens(tc)  = penst
        simvars%tau (tc)  = f_tau(Kt, zeta(zt))
        simvars%welf(tc)  = sum(val_j1_t*Phi(:,:,1))
        simvars%err_aggr(tc)   = f_aggregate_diff(simvars%output(tc),  simvars%invest(tc), simvars%C(tc),simvars%bequests(tc))
        simvars%B(tc)          = f_excessbonds(mut)
        simvars%err_income(tc) = f_income_diff(simvars%K(tc), zeta(zt), simvars%r(tc), simvars%rf(tc), delta(zt))

        if (.not. partial_equilibrium) then
            simvars%mu(tc)   = mut
            simvars%rf(tc+1) = f_riskfree_rate(Knew(tc+1),mut,pi_z(zt,:))
        elseif (partial_equilibrium .and.  mut == 0.0) then ! no aggregate risk case
            simvars%rf(tc+1) = f_riskfree_rate(simvars%K(tc+1),mut,pi_z(zt,:))
        endif
        r_pf = sign(1.0,apgridt)*(simvars%rf(tc+1) + kappat*simvars%mu(tc))/(1.0+g)
        simvars%r_pf_median(tc) = valnth(pack(r_pf, Phi/=0.0), ceiling(size(pack(r_pf, Phi/=0.0))/2.0))

        ! The next calculation is neglecting sign(1.0,apgridt), but that would become unnecessarily tedious
        simvars%r_pf_kappa_med(tc)=(simvars%rf(tc+1) + valnth(pack(kappat,Phi/=0.0), ceiling(size(pack(kappat, Phi/=0.0))/2.0)) *simvars%mu(tc))/(1.0+g)

        if (.not. calc_euler_errors) then
            simvars%eul_err_max(tc)=0.0
            simvars%eul_err_avg(tc)=0.0
        else
            if (tc <= t_scrap) then
            ! do not calculate Euler equation errors, because that is very costly and will be thrown away
                simvars%eul_err_max(tc)=0.0
                simvars%eul_err_avg(tc)=0.0
            else
                eul_err_temp = f_euler_errors(zt, simvars%rf(tc+1), mut,simvars%K(tc+1),coeffs, agg_grid, policies, value, xgridt, apgridt, kappat, Phi)
                simvars%eul_err_max(tc)=eul_err_temp(1)
                simvars%eul_err_avg(tc)=eul_err_temp(2)
            endif
        endif

        call calc_inequality_measures(simvars, xgridt, apgridt, stockst, Phi, etagrid(:,zt), penst, netwaget, tc)

        ! Average life cycle profiles and average Phi
        if (tc > t_scrap) then ! 'Throw away' first t_scrap
            ap_lct      = sum(sum(apgridt * Phi,1),1)
            stocks_lct  = sum(sum(stockst * Phi,1),1)
            cons_lct    = sum(sum((xgridt-apgridt) * Phi,1),1)
            return_lct  = sum(sum(Phi*r_pf,1),1)
            lc%ap       = lc%ap    + ap_lct    /(nt-t_scrap)
            lc%cons     = lc%cons  + cons_lct  /(nt-t_scrap)
            lc%stock    = lc%stock + stocks_lct/(nt-t_scrap)
            do jc=1,nj
                cons_var_lct(jc)   = sum(((xgridt(:,:,jc)-apgridt(:,:,jc)) - cons_lct(jc))**2 * Phi(:,:,jc))
                return_var_lct(jc)= sum((sign(1.0,apgridt(:,:,jc))*(simvars%rf(tc+1) + kappat(:,:,jc)*simvars%mu(tc))/(1.0+g) - return_lct(jc))**2 * Phi(:,:,jc))
            enddo
            lc%cons_var = lc%cons_var   + cons_var_lct  /(nt-t_scrap)
            lc%return   = lc%return     + return_lct    /(nt-t_scrap)
            lc%return_var=lc%return_var + return_var_lct/(nt-t_scrap)
            Phi_avg     = Phi_avg + Phi
        endif

    enddo
    where (lc%ap .ne. 0.0)
        lc%kappa = lc%stock/lc%ap
    elsewhere
        lc%kappa = 0.0
    end where

    if (partial_equilibrium) simvars%K = Knew
    if (exogenous_xgrid) Phi = Phi_avg/(nt-t_scrap)    ! else averaging does not make sense

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - pure real(dp) function f_excessbonds(mut)
!-------------------------------------------------------------------------------

    pure real(dp) function f_excessbonds(mut)
    ! calculates excess bond demand
        use params_mod    ,only: de_ratio
    !    use bs3vl_int   ! IMSL Math.pdf, p. 754ff: evaluate a 3d tensor-product spline given B-spline-coeffs

        real(dp) ,intent(in) :: mut       ! expected equity premium
        real(dp) ,dimension(:,:,:) ,allocatable :: Phit, apgridtt, stockstt, xgridtt    ! temporary distribution, policies for given z, K, AND mu
        real(dp)                         :: w, bond_supply
        integer                          :: i

        if (size(agg_grid%mu) >1) then
            i        = f_locate(agg_grid%mu, mut)
            w        = (mut - agg_grid%mu(i)) / (agg_grid%mu(i+1) - agg_grid%mu(i))

            if (.not. exogenous_xgrid) then
                xgridtt   = (1-w)* xgrid_zk(:,:,:,i) + w* xgrid_zk(:,:,:,i+1)
                Phit      = TransitionPhi(simvars%rf(tc),rt,netwaget,penst,simvars%bequests(tc),xgridtt,apgridt,stockst,etagrid(:,zt), Phi)
            else
                Phit = Phi
            endif

            apgridtt  = (1-w)*apgrid_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)
            stockstt  = (1-w)*stocks_zk(:,:,:,i) + w*stocks_zk(:,:,:,i+1)

        else
            Phit = Phi
            apgridtt  = apgrid_zk(:,:,:,1)
            stockstt  = stocks_zk(:,:,:,1)
        endif

        bond_supply   = de_ratio * sum(stockstt*Phit)/L_N_ratio
        f_excessbonds = bond_supply - sum((apgridtt - stockstt)*Phit)/L_N_ratio    ! analytically, L_N_ratio drops out

    end function f_excessbonds

end subroutine simulate
!-------------------------------------------------------------------------------

pure subroutine calc_inequality_measures(simvars, xgridt, apgridt, stockst, Phi, etagridt, penst, netwaget, tc)
    use statistics  ,only: lorenz_calc
    use params_mod  ,only: t_scrap, jr, ej
    type(tSimvars)             ,intent(inout) :: simvars
    real(dp) ,dimension(:,:,:) ,intent(in)    :: xgridt, apgridt, stockst, Phi
    real(dp) ,dimension(:)     ,intent(in)    :: etagridt
    real(dp)                   ,intent(in)    :: penst, netwaget
    integer                    ,intent(in)    :: tc
    real(dp) ,dimension(:), allocatable       :: lorenz_x, lorenz_y
    real(dp) ,dimension(:,:,:) ,allocatable    :: const
    real(dp) ,dimension(size(Phi,2),size(Phi,3)) :: income_dist
    real(dp)                                  :: mean, std
    integer                                   :: error, jc
    logical                                   :: ms_equilibrium

    if (size(simvars%rf) < t_scrap) then
        ms_equilibrium = .true.
    else
        ms_equilibrium = .false.
    endif

    if (tc <= t_scrap .and. .not. ms_equilibrium) then
    ! do not calculate Ginis errors, because that is very costly and will be thrown away
        simvars%gini_income(tc)      = 0.0
        simvars%gini_assets(tc)      = 0.0
        simvars%gini_stocks(tc)      = 0.0
        simvars%gini_consumption(tc) = 0.0
        simvars%cv_income(tc)        = 0.0
        simvars%cv_assets(tc)        = 0.0
        simvars%cv_stocks(tc)        = 0.0
        simvars%cv_consumption(tc)   = 0.0
    else

        const = xgridt -apgridt

        do jc=1,size(income_dist,2)
            if (jc>=jr) then
                income_dist(:,jc) = penst
            else
                income_dist(:,jc) = netwaget*ej(jc)*etagridt
            endif
        enddo

        lorenz_x = pack(income_dist, .true.)
        lorenz_y = lorenz_x
        call lorenz_calc(pack(sum(Phi,1), .true.), lorenz_x, lorenz_y, simvars%gini_income(tc), error)
        lorenz_x = pack(apgridt, .true.)
        lorenz_y = lorenz_x
        call lorenz_calc(pack(Phi, .true.), lorenz_x, lorenz_y, simvars%gini_assets(tc), error)
        lorenz_x = pack(stockst, .true.)
        call lorenz_calc(pack(Phi, .true.), lorenz_x, lorenz_y, simvars%gini_stocks(tc), error)
        lorenz_x = pack(const, .true.)
        call lorenz_calc(pack(Phi, .true.), lorenz_x, lorenz_y, simvars%gini_consumption(tc), error)

        mean = sum(sum(Phi,1)*income_dist)
        if (mean == 0.0) mean = 1.0
        std  = sqrt(sum(sum(Phi,1)*(income_dist - mean)**2))
        simvars%cv_income(tc)        = std/mean

        mean = sum(Phi*apgridt)
        if (mean == 0.0) mean = 1.0
        std  = sqrt(sum(Phi*(apgridt - mean)**2))
        simvars%cv_assets(tc)        = std/mean

        mean = sum(Phi*stockst)
        if (mean == 0.0) mean = 1.0
        std  = sqrt(sum(Phi*(stockst - mean)**2))
        simvars%cv_stocks(tc)        = std/mean

        mean = sum(Phi*const)
        if (mean == 0.0) mean = 1.0
        std  = sqrt(sum(Phi*(const - mean)**2))
        simvars%cv_consumption(tc)        = std/mean

    endif

end subroutine calc_inequality_measures
!-------------------------------------------------------------------------------

pure function f_euler_errors(zt, rfp, mut,kp,coeffs, grids, policies, value, xgridt, apgridt, kappat, Phi)
    use params_mod ,only: beta, gamm, g, theta, jr, surv, ej, etagrid, cmin
    use household_solution_mod ,only: interp_policies_tomorrow, consumption
    use laws_of_motion ,only: Forecast_mu
    use income ,only: f_stock_return, f_pensions, f_netwage, zeta, delta

    real(dp), dimension(2) :: f_euler_errors
    integer                          ,intent(in) :: zt
    real(dp)                         ,intent(in) :: rfp, mut, kp
    type(tCoeffs)                    ,intent(in) :: coeffs ! coefficients for laws of motion
    type(tAggGrids)                  ,intent(in) :: grids  ! grids for aggregate states k and mu
    type(tPolicies)                  ,intent(in) :: policies
    real(dp) ,dimension(:,:,:,:,:,:) ,intent(in) :: value
    real(dp) ,dimension(:,:,:)       ,intent(in) :: xgridt, apgridt, kappat, Phi

    real(dp) :: betatildej, app_min, evp
    real(dp) ,allocatable :: mup(:), rp(:), yp(:,:)
    real(dp) ,allocatable ,dimension(:,:,:) :: consp, xgridp, vp, cons_opt, cons_t, eul_err
    integer :: zpc, jc, ec, xc, nj, nz, n_eta, nx, nmu
    logical(1) ,dimension(size(coeffs%mu,2)) :: error

    nz = size(coeffs%mu,2)
    nj = size(policies%apgrid,4)
    n_eta = size(etagrid,1)
    nx = size(policies%apgrid,1)
    nmu = size(grids%mu)

    allocate(mup(nz), rp(nz), yp(nz,n_eta))
    allocate(consp(nx,n_eta,nz), xgridp(nx,n_eta,nz), vp(nx,n_eta,nz))
    allocate(cons_opt(nx,n_eta,nj), cons_t(nx,n_eta,nj), eul_err(nx,n_eta,nj))

    do zpc = 1,nz
        mup(zpc) = Forecast_mu(coeffs%mu(:,zpc), kp, mut)
        rp(zpc) = f_stock_return(kp, zeta(zpc), delta(zpc), rfp)
    enddo
    where (mup > grids%mu(nmu)) mup = grids%mu(nmu) ! because I do this in household_solution_mod:calc_vars_tomorrow
    where (mup < grids%mu(1))   mup = grids%mu(1)

    cons_t=xgridt -apgridt

    do jc=1,nj-1
        betatildej = beta*surv(jc)**(1.0/gamm)*(1.0+g)**((1.0-theta)/gamm)

        do zpc = 1,nz
            if (jc+1>=jr) then
                yp(:,zpc) = f_pensions(kp, zeta(zpc))
            else
                yp(:,zpc) = ej(jc+1) * f_netwage(kp, zeta(zpc)) * etagrid(:,zpc)
            endif
        enddo

        call interp_policies_tomorrow(policies, policies%consumption(), value, kp, mup, grids, jc, consp, xgridp, vp, app_min)
        do ec=1,n_eta
            do xc=1,nx
                if (Phi(xc,ec,jc)==0.0) then
                    cons_opt(xc,ec,jc)=cons_t(xc,ec,jc)
                else
                    call consumption(apgridt(xc,ec,jc), kappat(xc,ec,jc), xgridp, consp, vp, rfp,rp, yp, zt, xc, ec, betatildej, cons_opt(xc,ec,jc), evp, error)
                    if ((cons_opt(xc,ec,jc) < cmin*10.0_dp) .or. any(error)) cons_opt(xc,ec,jc) = cons_t(xc,ec,jc)
                endif
            enddo
        enddo
    enddo
    cons_opt(:,:,nj)=cons_t(:,:,nj)

    eul_err = abs(1.0-cons_opt/cons_t)

    f_euler_errors(1) = maxval(eul_err)
    !f_euler_errors(2) = sum(eul_err)/(size(eul_err)-zero_mass)
    f_euler_errors(2) = sum(Phi*eul_err) ! average
end function f_euler_errors
!-------------------------------------------------------------------------------

end module simulation_mod
