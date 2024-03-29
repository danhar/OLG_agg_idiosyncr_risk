!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

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
    use params_mod      ,only: n,g,L_N_ratio,pop_frac,pi_z,etagrid,t_scrap,exogenous_xgrid, &
                               partial_equilibrium, zeta, delta, alpha, check_dynamic_efficiency, &
                               tol_mut=> tol_simulation_marketclearing, nx_factor
    use income          ,only: f_netwage, f_transfers, f_pensions, f_stock_return, f_riskfree_rate, f_tau, f_net_mpk
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
    real(dp) ,dimension(:,:,:,:) ,allocatable :: apgrid_zk, kappa_zk, xgrid_zk, stocks_zk, value_zk    ! policies for given z and K
    real(dp) ,dimension(:,:,:)   ,allocatable :: apgridt, kappat, xgridt, stockst, valuet, const, exp_value_t, weight ! policies for given z, K, and mu
    real(dp) ,dimension(:,:,:)   ,allocatable :: Phi_avg, r_pf ! portfolio return
    real(dp) ,dimension(:)       ,allocatable :: ap_lct, stocks_lct, cons_lct, cons_var_lct, return_lct, return_var_lct, log_cons_lct, var_log_cons_lct
    real(dp) ,dimension(:)       ,allocatable :: Knew       ! partial equilibrium: save aggregate stock in t
    real(dp)  :: Kt, mut, rt, netwaget, transfert, penst, w, eul_err_temp(2) ! variables in period t
    integer   :: tc, i, zt, jc, nmu, nx, n_eta, nj, nt, nk, dyn_eff_b_counter

    ! Intel Fortran Compiler XE 13.0 Update 1 (and previous) has a bug on realloc on assignment. If that is corrected, I think I can remove this whole allocation block.
    ! Intel Compiler 2017.0.098 works if I remove the block, but performance degrades substantially.
    nmu = size(agg_grid%mu); nk= size(agg_grid%k); nx=size(value,1); n_eta=size(value,2); nj=size(value,4); nt=size(simvars%z)
    allocate(apgrid_zk(nx,n_eta,nj,nmu), kappa_zk(nx,n_eta,nj,nmu), xgrid_zk(nx,n_eta,nj,nmu), stocks_zk(nx,n_eta,nj,nmu), value_zk(nx,n_eta,nj,nmu))
    allocate(apgridt(nx,n_eta,nj), kappat(nx,n_eta,nj), xgridt(nx,n_eta,nj), const(nx,n_eta,nj), stockst(nx,n_eta,nj), valuet(nx,n_eta,nj))
    allocate(exp_value_t(nx,n_eta,nj), weight(nx,n_eta,nj), Phi_avg(nx,n_eta,nj), r_pf(nx,n_eta,nj))
    allocate(ap_lct(nj), stocks_lct(nj), cons_lct(nj), cons_var_lct(nj), return_lct(nj), return_var_lct(nj), log_cons_lct(nj), var_log_cons_lct(nj))
    allocate(Knew(nt+1))

    Phi_avg         = 0.0
    call lc%allocate(nj,nx)
    call lc%set_number(0.0_dp)
    simvars%err_K   = .false.
    simvars%err_mu  = .false.
    dyn_eff_b_counter = 0
    Knew(1)   = simvars%K(1)    ! This is only interesting for PE

    do tc=1,nt
        zt = simvars%z(tc)
        Kt = simvars%K(tc)
        if (Kt <= agg_grid%K(1) .or. Kt >= agg_grid%K(nk)) simvars%err_K(tc) = .true.

        ! Prices
        netwaget = f_netwage(Kt, zeta(zt))
        transfert= f_transfers(Kt, zeta(zt))
        penst    = f_pensions(Kt, zeta(zt))
        rt       = f_stock_return(Kt, zeta(zt), delta(zt), simvars%rf(tc))

        if (tc == 1) then
            simvars%bequests(tc) = 0.0
        else
            simvars%bequests(tc) = f_bequests(simvars%rf(tc), rt, stockst, apgridt, Phi)
        endif

        ! Projection of policies on Kt
        i        = f_locate(agg_grid%K, Kt)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        w        = (Kt - agg_grid%K(i)) / (agg_grid%K(i+1) - agg_grid%K(i))
        ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds
        !! Once the Intel Fortran bug on reallocation is fixed, I can simply write xgrid_zk= ...
        xgrid_zk (:,:,:,:)= (1-w)*policies%xgrid (:,:,zt,:,i,:) +w*policies%xgrid (:,:,zt,:,i+1,:)
        apgrid_zk(:,:,:,:)= (1-w)*policies%apgrid(:,:,zt,:,i,:) +w*policies%apgrid(:,:,zt,:,i+1,:)
        stocks_zk(:,:,:,:)= (1-w)*policies%stocks(:,:,zt,:,i,:) +w*policies%stocks(:,:,zt,:,i+1,:)
        value_zk (:,:,:,:)= (1-w)*value          (:,:,zt,:,i,:) +w*value          (:,:,zt,:,i+1,:)
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
            Phi      = TransitionPhi(simvars%rf(tc),rt,netwaget,transfert,penst,simvars%bequests(tc),xgridt,apgridt,stockst,etagrid(:,zt), Phi)
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
		        Phi      = TransitionPhi(simvars%rf(tc),rt,netwaget,transfert,penst,simvars%bequests(tc),xgridt,apgridt,stockst,etagrid(:,zt), Phi)
	        endif

	        apgridt  = (1-w)*apgrid_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)
	        stockst  = (1-w)*stocks_zk(:,:,:,i) + w*stocks_zk(:,:,:,i+1)
	        valuet   = (1-w)*value_zk (:,:,:,i) + w*value_zk (:,:,:,i+1)
	        where (apgridt .ne. 0.0)
	            kappat = stockst/ apgridt
	        elsewhere ! includes apgrid(:,nj) =0
	            kappat = 0.0
	        end where
        else
            apgridt  = apgrid_zk(:,:,:,1)
            stockst  = stocks_zk(:,:,:,1)
            valuet   = value_zk (:,:,:,1)
            kappat   = kappa_zk (:,:,:,1)
        endif
        const        = xgridt - apgridt
        exp_value_t  = valuet*Phi(:,:,:)

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

        simvars%output(tc)     = zeta(zt)*simvars%K(tc)**alpha
        simvars%stock(tc)      = sum(stockst*Phi)/L_N_ratio ! different from K(t+1) since it is in today's per capita terms! Thus no (1+g)(1+n) in denominator.
        simvars%bonds(tc)      = sum((apgridt-stockst)*Phi)/L_N_ratio
        simvars%invest(tc)     = simvars%stock(tc) + simvars%bonds(tc) -simvars%K(tc)*(1.0-delta(zt))
        simvars%C(tc)          = sum((xgridt-apgridt)*Phi) / L_N_ratio ! in units of efficient labor
        simvars%net_mpk(tc)    = f_net_mpk(Kt, zeta(zt), delta(zt))
        simvars%r(tc)     = rt
        simvars%wage(tc)  = netwaget
        simvars%trans(tc) = transfert
        simvars%pens(tc)  = penst
        simvars%tau (tc)  = f_tau(Kt, zeta(zt))
        simvars%welf(tc)  = sum(exp_value_t(:,:,1))
        simvars%err_aggr(tc)   = f_aggregate_diff(simvars%output(tc),  simvars%invest(tc), simvars%C(tc),simvars%bequests(tc))
        simvars%B(tc)          = f_excessbonds(mut)
        simvars%err_income(tc) = f_income_diff(simvars%K(tc), zeta(zt), simvars%r(tc), simvars%rf(tc), delta(zt))

        if (.not. partial_equilibrium) then
            simvars%mu(tc)   = mut
            simvars%rf(tc+1) = f_riskfree_rate(Knew(tc+1),mut,pi_z(zt,:))
        elseif (partial_equilibrium .and.  mut == 0.0) then ! no aggregate risk case
            simvars%rf(tc+1) = f_riskfree_rate(simvars%K(tc+1),mut,pi_z(zt,:))
        endif
        ! r_pf = sign(1.0,apgridt)*(simvars%rf(tc+1) + kappat*simvars%mu(tc))/(1.0+g)
        r_pf = sign(1.0,apgridt)*(simvars%rf(tc+1) + kappat*simvars%mu(tc))
        simvars%r_pf(tc) = sum(r_pf*Phi)

        if (nx_factor < 5) then
            ! The following two medians for portfolio holdings may stop program execution for nx_factor > 8 due to memory limits.
            simvars%r_pf_median(tc) = valnth(pack(r_pf, Phi/=0.0), ceiling(size(pack(r_pf, Phi/=0.0))/2.0))
            ! The next calculation is neglecting sign(1.0,apgridt), but that would become unnecessarily tedious
            simvars%r_pf_kappa_med(tc)=(simvars%rf(tc+1) + valnth(pack(kappat,Phi/=0.0), ceiling(size(pack(kappat, Phi/=0.0))/2.0)) *simvars%mu(tc))
        else
            simvars%r_pf_median(tc) = 0.0
            simvars%r_pf_kappa_med(tc)= 0.0
        endif

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

        call calc_inequality_measures(simvars, xgridt, apgridt, stockst, Phi, etagrid(:,zt), penst, netwaget, transfert, tc)

        if (check_dynamic_efficiency .and. .not. partial_equilibrium) then
            dyn_eff_b_counter = dyn_eff_b_counter + 1
            simvars%dyn_eff_b(tc) = dyn_eff_b_counter
            if ((1.0 + simvars%rf(tc+1))>(1+n)*(1+g)) then
                ! dyn_eff_a defined as violation of condition A.
                simvars%dyn_eff_a(tc) = dyn_eff_a(simvars%rf(tc+1), simvars%K(tc+1), stockst, apgridt, policies, agg_grid, Phi)
                dyn_eff_b_counter = 0
            else
                simvars%dyn_eff_a(tc) = .false. ! no violation of condition A.
            endif
        endif

        ! Average life cycle profiles and average Phi
        if (tc > t_scrap) then ! 'Throw away' first t_scrap
            do jc=1,nj
                weight(:,:,jc) = Phi(:,:,jc)/pop_frac(jc)
            enddo

            ap_lct      = sum(sum(apgridt    * weight,1),1)
            stocks_lct  = sum(sum(stockst    * weight,1),1)
            cons_lct    = sum(sum(const      * weight,1),1)
            log_cons_lct= sum(sum(log(const) * weight,1),1)
            return_lct  = sum(sum(r_pf       * weight,1),1)
            do jc=1,nj
                cons_var_lct(jc)     = sum((const(:,:,jc) - cons_lct(jc))**2 * weight(:,:,jc))
                ! According to Wikipedia article "Algorithms for calculating variance" this alternative is more unstable:
                !cons_var_lct(jc)    = sum(const(:,:,jc)**2 * weight(:,:,jc)) - cons_lct(jc)**2
                var_log_cons_lct(jc) = sum((log(const(:,:,jc)) - log_cons_lct(jc))**2 * weight(:,:,jc))
                return_var_lct(jc)   = sum((sign(1.0,apgridt(:,:,jc))*(simvars%rf(tc+1) + kappat(:,:,jc)*simvars%mu(tc)) - return_lct(jc))**2 * weight(:,:,jc))
            enddo
            lc%ap           = lc%ap           + ap_lct          /(nt-t_scrap)
            lc%cons         = lc%cons         + cons_lct        /(nt-t_scrap)
            lc%stock        = lc%stock        + stocks_lct      /(nt-t_scrap)
            lc%cons_var     = lc%cons_var     + cons_var_lct    /(nt-t_scrap)
            lc%return       = lc%return       + return_lct      /(nt-t_scrap)
            lc%return_var   = lc%return_var   + return_var_lct  /(nt-t_scrap)
            lc%log_cons     = lc%log_cons     + log_cons_lct    /(nt-t_scrap)
            lc%var_log_cons = lc%var_log_cons + var_log_cons_lct/(nt-t_scrap)
            lc%exp_value    = lc%exp_value    + exp_value_t     /(nt-t_scrap)
            lc%xgrid        = lc%xgrid        + xgridt          /(nt-t_scrap)
            Phi_avg         = Phi_avg + Phi
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
                Phit      = TransitionPhi(simvars%rf(tc),rt,netwaget,transfert,penst,simvars%bequests(tc),xgridtt,apgridt,stockst,etagrid(:,zt), Phi)
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

pure subroutine calc_inequality_measures(simvars, xgridt, apgridt, stockst, Phi, etagridt, penst, netwaget, transfert, tc)
    use statistics  ,only: lorenz_calc
    use params_mod  ,only: t_scrap, jr, ej, trans_grid, trans_prob
    type(tSimvars)             ,intent(inout) :: simvars
    real(dp) ,dimension(:,:,:) ,intent(in)    :: xgridt, apgridt, stockst, Phi
    real(dp) ,dimension(:)     ,intent(in)    :: etagridt
    real(dp)                   ,intent(in)    :: penst, netwaget, transfert
    integer                    ,intent(in)    :: tc
    real(dp) ,dimension(:), allocatable       :: lorenz_x, lorenz_y
    real(dp) ,dimension(:,:,:) ,allocatable   :: const, assets_t, Phi_trans, income_dist
    real(dp) ,dimension(:,:)   ,allocatable   :: income_dist_a ! just for calc assets
    real(dp)                                  :: mean, std
    integer                                   :: error, jc, ec, n_trans, n_eta, n_j
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
        n_trans = size(trans_grid)
        n_eta = size(Phi,2)
        n_j   = size(Phi,3)
        allocate(income_dist(n_trans, n_eta, n_j))
        allocate(income_dist_a(n_eta, n_j))
        income_dist   = 0.0
        Phi_trans     = income_dist
        income_dist_a = 0.0
        const = xgridt -apgridt

        do jc=1, n_j
            do ec=1,n_eta
                if (jc>=jr) then
                    income_dist_a(ec,jc) = 0.0 !penst
                    income_dist(:,ec,jc) = 0.0 !penst
                    Phi_trans(:,ec,jc) = sum(Phi(:,ec,jc))/real(n_trans,dp)
                else
                    income_dist_a(ec,jc) = netwaget*ej(jc)*etagridt(ec) + transfert
                    income_dist(:,ec,jc) = netwaget*ej(jc)*etagridt(ec)*trans_grid + transfert
                    Phi_trans(:,ec,jc) = sum(Phi(:,ec,jc))*trans_prob
                endif
            enddo
        enddo

        ! Only here is income_dist_a used.
        assets_t = xgridt - spread(income_dist_a, 1, size(xgridt,1))

        lorenz_x = pack(income_dist, .true.)
        lorenz_y = lorenz_x
        call lorenz_calc(pack(Phi_trans, .true.), lorenz_x, lorenz_y, simvars%gini_income(tc), error)
        lorenz_x = pack(assets_t, .true.)
        lorenz_y = lorenz_x
        call lorenz_calc(pack(Phi, .true.), lorenz_x, lorenz_y, simvars%gini_assets(tc), error)
        lorenz_x = pack(stockst, .true.)
        call lorenz_calc(pack(Phi, .true.), lorenz_x, lorenz_y, simvars%gini_stocks(tc), error)
        lorenz_x = pack(const, .true.)
        call lorenz_calc(pack(Phi, .true.), lorenz_x, lorenz_y, simvars%gini_consumption(tc), error)

        mean = sum(Phi_trans*income_dist)
        if (mean == 0.0) mean = 1.0
        std  = sqrt(sum(Phi_trans*(income_dist - mean)**2))
        simvars%cv_income(tc)        = std/mean

        mean = sum(Phi*assets_t)
        if (mean == 0.0) mean = 1.0
        std  = sqrt(sum(Phi*(assets_t - mean)**2))
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
    use params_mod ,only: beta, gamm, g, theta, jr, surv, ej, etagrid, trans_grid, cmin
    use household_solution_mod ,only: interp_policies_tomorrow, consumption
    use laws_of_motion ,only: Forecast_mu
    use income ,only: f_stock_return, f_pensions, f_netwage, f_transfers, zeta, delta

    real(dp), dimension(2) :: f_euler_errors
    integer                          ,intent(in) :: zt
    real(dp)                         ,intent(in) :: rfp, mut, kp
    type(tCoeffs)                    ,intent(in) :: coeffs ! coefficients for laws of motion
    type(tAggGrids)                  ,intent(in) :: grids  ! grids for aggregate states k and mu
    type(tPolicies)                  ,intent(in) :: policies
    real(dp) ,dimension(:,:,:,:,:,:) ,intent(in) :: value
    real(dp) ,dimension(:,:,:)       ,intent(in) :: xgridt, apgridt, kappat, Phi

    real(dp) :: betatildej, app_min, evp
    real(dp) ,allocatable :: mup(:), rp(:), yp(:,:,:)
    real(dp) ,allocatable ,dimension(:,:,:) :: consp, xgridp, vp, cons_opt, cons_t, eul_err
    integer :: zpc, jc, ec, xc, nj, nz, n_eta, n_trans, nx, nmu
    logical(1) ,dimension(size(coeffs%mu,2)) :: error

    nz = size(coeffs%mu,2)
    nj = size(policies%apgrid,4)
    n_eta = size(etagrid,1)
    n_trans = size(trans_grid)
    nx = size(policies%apgrid,1)
    nmu = size(grids%mu)

    allocate(mup(nz), rp(nz), yp(n_trans, n_eta,nz))
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
                yp(:,:,zpc) = f_pensions(kp, zeta(zpc))
            else
                do ec=1,n_eta
                    yp(:,ec,zpc) = ej(jc+1) * f_netwage(kp, zeta(zpc)) * etagrid(ec,zpc) * trans_grid + f_transfers(kp, zeta(zpc))
                enddo
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

pure logical function dyn_eff_a(rf_t, Kt, stockst, apgridt, policies, agg_grid, Phi)
    ! Checking the Demange 2002 criterion for dynamic efficiency. Specifically, here we check
    ! condition (a) of Prop. 1 of Kubler and Krueger AER 2006.
    ! Not implemented for partial equilibrium (since that is efficient if GE is).
    ! See my note in /home/elessar/work/Research/EPSocSec/notes/170515-Dynamic_efficiency.
    ! Important: Output .true. indicates a violation of the condition!

    use params_mod      ,only: n,g,L_N_ratio,pi_z,etagrid,exogenous_xgrid, zeta, delta, tol_mut=> tol_simulation_marketclearing
    use income          ,only: f_netwage, f_transfers, f_pensions, f_stock_return, f_riskfree_rate
    use fun_locate      ,only: f_locate
    use distribution    ,only: TransitionPhi, CheckPhi
    use fun_aggregate_diff
    use fun_zbrent

    real(dp)         ,intent(in)    :: rf_t, Kt ! these are actually t+1, but for notational purposes just t.
    real(dp) ,dimension(:,:,:) ,intent(in) :: apgridt, stockst, Phi ! policies for given z, K, and mu
    type(tPolicies)  ,intent(in)    :: policies
    type(tAggGrids)  ,intent(in)    :: agg_grid
    real(dp) ,dimension(:,:,:,:) ,allocatable :: apgrid_zk, xgrid_zk, stocks_zk    ! policies for given z and K
    real(dp) ,dimension(:,:,:)   ,allocatable :: apgrid_new, xgridt, Phi_new
    real(dp)  :: Knew, rf_new, mut, rt, bequestst, netwaget, transfert, penst, w ! variables in period t
    integer   :: i, zt, nmu, nx, n_eta, nj, nk, nz
    logical   :: r_low, r_high

    nmu = size(agg_grid%mu); nk= size(agg_grid%k); nx=size(policies%apgrid,1); n_eta=size(policies%apgrid,2); nj=size(policies%apgrid,4); nz=size(zeta)
    allocate(apgrid_zk(nx,n_eta,nj,nmu), xgrid_zk(nx,n_eta,nj,nmu), stocks_zk(nx,n_eta,nj,nmu))
    allocate(xgridt, apgrid_new, Phi_new, source=apgridt)
    r_low = .false.
    r_high= .false.

    do zt = 1, nz ! for all shocks tomorrow
        ! Prices
        netwaget = f_netwage(Kt, zeta(zt))
        transfert= f_transfers(Kt, zeta(zt))
        penst    = f_pensions(Kt, zeta(zt))
        rt       = f_stock_return(Kt, zeta(zt), delta(zt), rf_t)
        bequestst = f_bequests(rf_t, rt, stockst, apgridt, Phi)

        ! Projection of policies on Kt
        i        = f_locate(agg_grid%K, Kt)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
        w        = (Kt - agg_grid%K(i)) / (agg_grid%K(i+1) - agg_grid%K(i))
        ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds
        xgrid_zk (:,:,:,:)= (1-w)*policies%xgrid (:,:,zt,:,i,:) +w*policies%xgrid (:,:,zt,:,i+1,:)
        apgrid_zk(:,:,:,:)= (1-w)*policies%apgrid(:,:,zt,:,i,:) +w*policies%apgrid(:,:,zt,:,i+1,:)
        stocks_zk(:,:,:,:)= (1-w)*policies%stocks(:,:,zt,:,i,:) +w*policies%stocks(:,:,zt,:,i+1,:)

ex:     if (exogenous_xgrid .or. (nmu ==1) ) then
            xgridt   = xgrid_zk(:,:,:,1)
            ! To calc distribution, we need xgridt , along with netwaget etc, of this period (tc), but policy projections of last
            Phi_new  = TransitionPhi(rf_t,rt,netwaget,transfert,penst,bequestst,xgridt,apgridt,stockst,etagrid(:,zt), Phi)
        endif ex

        ! Find mut that clears bondmarket
        if (f_excessbonds(agg_grid%mu(1)) * f_excessbonds(agg_grid%mu(nmu)) < 0.0 ) then
            mut  = f_zbrent(f_excessbonds,agg_grid%mu(1),agg_grid%mu(nmu),tol_mut)
        else
            ! simvars%err_mu(tc) = .true.
            if (abs(f_excessbonds(agg_grid%mu(1))) < abs(f_excessbonds(agg_grid%mu(nmu)))) then
                mut = agg_grid%mu(1)
            else
                mut = agg_grid%mu(nmu)
            endif
        endif

        ! Projection of policies on mut
        if (nmu > 1) then
            i        = f_locate(agg_grid%mu, mut)
            w        = (mut - agg_grid%mu(i)) / (agg_grid%mu(i+1) - agg_grid%mu(i))
            if (.not. exogenous_xgrid) then
                xgridt   = (1-w)* xgrid_zk(:,:,:,i) + w* xgrid_zk(:,:,:,i+1)
                ! To calc distribution, we need xgridt , along with netwaget etc, of this period (tc), but policy projections of last
                Phi_new = TransitionPhi(rf_t,rt,netwaget,transfert,penst,bequestst,xgridt,apgridt,stockst,etagrid(:,zt), Phi)
            endif
            apgrid_new  = (1-w)*apgrid_zk(:,:,:,i) + w*apgrid_zk(:,:,:,i+1)
        else
            apgrid_new  = apgrid_zk(:,:,:,1)
        endif

        ! call CheckPhi(Phi, simvars%Phi_1(tc), simvars%Phi_nx(tc)) ! Do here because Phi now was transitioned for all cases

        ! Aggregate
        Knew = sum(apgrid_new*Phi_new) /(L_N_ratio*(1.0+n)*(1.0+g))  ! in units of efficient labor
        ! Could allow some deviation out of grid - recorded as error in next period. Also in PE to keep comparable!
!        if (Knew     < agg_grid%K(1 )) then
!            Knew     = agg_grid%K(1 )
!        elseif (Knew > agg_grid%K(nk)) then
!            Knew     = agg_grid%K(nk)
!        endif

        rf_new = f_riskfree_rate(Knew,mut,pi_z(zt,:))

        if ((1.0 + rf_new) > (1.0+n)*(1.0+g)) then
            if (rt > rf_new) then
                r_high = .true.
            else
                r_low  = .true.
            endif
        endif

        if (r_low .and. r_high) exit
    enddo

dyn_eff_a = .not.(r_low .and. r_high) ! report as a violation

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
                Phit      = TransitionPhi(rf_t,rt,netwaget,transfert,penst,bequestst,xgridtt,apgridt,stockst,etagrid(:,zt), Phi)
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

end function dyn_eff_a
!-------------------------------------------------------------------------------

end module simulation_mod
