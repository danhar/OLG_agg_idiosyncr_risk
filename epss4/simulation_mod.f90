module simulation_mod
    implicit none
contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure subroutine simulate(policies, agg_grid, simvars, Phi, lc)
! - subroutine print_error_msg_sim(simvars)
!-------------------------------------------------------------------------------

pure subroutine simulate(policies, value, agg_grid, simvars, Phi, lc)
! Performs the Krusell-Smith simulation step and records lifecycle statistics
    use kinds           ,only: dp
    use types           ,only: tSimvars, tLifecycle, AllocateType, set_number
    use policies_class  ,only: tPolicies
    use aggregate_grids ,only: tAggGrids
    use params_mod      ,only: n,g,L_N_ratio,pi_z,etagrid,t_scrap,exogenous_xgrid, partial_equilibrium, zeta, delta, alpha
    use income          ,only: f_netwage, f_pensions, f_stock_return, f_riskfree_rate, f_tau
    use fun_locate      ,only: f_locate
    use distribution    ,only: TransitionPhi, CheckPhi
    use fun_zbrent
    use fun_aggregate_diff
    use partial_sorting     ! function valnth

    type(tPolicies)  ,intent(in)    :: policies
    real(dp)         ,intent(in)    :: value(:,:,:,:,:,:)
    type(tAggGrids)  ,intent(in)    :: agg_grid
    type(tSimvars)   ,intent(inout) :: simvars    ! (zt, kt, mut, bt,...), first element contains starting values
    real(dp)         ,intent(inout) :: Phi(:,:,:) ! distribution. Returns: average Phi if (exogenous_xgrid), else Phi in nt
    type(tLifecycle) ,intent(out)   :: lc         ! lifecycle profiles
    real(dp) ,dimension(:,:,:,:) ,allocatable :: apgrid_zk, kappa_zk, xgrid_zk, stocks_zk ! policies for given z and K
    real(dp) ,dimension(:,:,:)   ,allocatable :: apgridt, kappat, xgridt, stockst         ! policies for given z, K, and mu
    real(dp) ,dimension(:,:,:)   ,allocatable :: Phi_avg, r_pf, val_j1_zk ! portfolio return
    real(dp) ,dimension(:,:)     ,allocatable :: val_j1_t ! value of j=1 for given z and K, and mu
    real(dp) ,dimension(:)       ,allocatable :: ap_lct, stocks_lct, cons_lct, cons_var_lct, return_lct, return_var_lct
    real(dp) ,dimension(:)       ,allocatable :: Knew       ! partial equilibrium: save aggregate stock in t
    real(dp) ,parameter :: tol_mut = 1e-8_dp
    real(dp)  :: Kt, mut, rt, netwaget, penst, w, brackl, bracku ! variables in period t
    integer   :: tc, i, zt, jc, nmu, nx, n_eta, nj, nt, nk
    logical   :: brack_found

    ! Intel Fortran Compiler XE 13.0 Update 1 (and previous) has a bug on realloc on assignment. If that is corrected, I think I can remove this whole allocation block
    nmu = size(agg_grid%mu); nk= size(agg_grid%k); nx=size(value,1); n_eta=size(value,2); nj=size(value,4); nt=size(simvars%z)
    allocate(apgrid_zk(nx,n_eta,nj,nmu), kappa_zk(nx,n_eta,nj,nmu), xgrid_zk(nx,n_eta,nj,nmu), stocks_zk(nx,n_eta,nj,nmu))
    allocate(apgridt(nx,n_eta,nj), kappat(nx,n_eta,nj), xgridt(nx,n_eta,nj), stockst(nx,n_eta,nj), Phi_avg(nx,n_eta,nj), r_pf(nx,n_eta,nj))
    allocate(val_j1_zk(nx,n_eta,nmu), val_j1_t(nx,n_eta))
    allocate(ap_lct(nj), stocks_lct(nj), cons_lct(nj), cons_var_lct(nj), return_lct(nj), return_var_lct(nj))
    allocate(Knew(nt+1))

    Phi_avg         = 0.0
    call set_number(lc, 0.0_dp)
    simvars%err_K   = .false.
    simvars%err_mu  = .false.
    Knew(1)   = simvars%K(1)    ! This is only interesting for PE

    ! Initialize rf and policies for 'previous' period (i.e. t=0) by setting them to mean-shock value
    call get_initial_values(apgridt, stockst, simvars%rf(1))

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

ex:     if (exogenous_xgrid .or. (nmu ==1) ) then
            xgridt   = xgrid_zk(:,:,:,1)
            ! To calc distribution, we need xgridt , along with netwaget etc, of this period (tc), but policy projections of last
            Phi      = TransitionPhi(simvars%rf(tc),rt,netwaget,penst,xgridt,apgridt,stockst,etagrid(:,zt), Phi)
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
		        Phi      = TransitionPhi(simvars%rf(tc),rt,netwaget,penst,xgridt,apgridt,stockst,etagrid(:,zt), Phi)
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
! - pure subroutine get_initial_values(apgridt, stockst, rf1)
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
                Phit      = TransitionPhi(simvars%rf(tc),rt,netwaget,penst,xgridtt,apgridt,stockst,etagrid(:,zt), Phi)
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
!-------------------------------------------------------------------------------

    pure subroutine get_initial_values(apgridt, stockst, rf1)
    ! Calculate policies for mean shock, which correspond to initial (=mean shock) distribution, K, mu

        use params_mod ,only: stat_dist_z
	    real(dp) ,dimension(:,:,:) ,intent(inout) :: apgridt, stockst !intent(in) to get allocation status and size
	    real(dp)                   ,intent(out)   :: rf1
	    real(dp) :: K0, mu0, wK, wmu, w_ms
	    integer  :: iK, imu, zc, jc, nz

        nz = size(policies%apgrid,3)
        w_ms    = 0.25_dp ! corresponds to MeanShockGE variable w
        apgridt = 0.0
        stockst = 0.0

	    K0      = simvars%K(1)
        iK      = f_locate(agg_grid%K, K0)     ! In 'default', returns ju-1 if x>xgrid(ju-1) !
        wK      = (K0 - agg_grid%K(iK)) / (agg_grid%K(iK+1) - agg_grid%K(iK))
        mu0     = (agg_grid%mu(1) + agg_grid%mu(nmu))/2  ! very close to mean shock value
        imu     = f_locate(agg_grid%mu, mu0)

        if (nmu > 1) then
            wmu     = (mu0 - agg_grid%mu(imu)) / (agg_grid%mu(imu+1) - agg_grid%mu(imu))
		    do zc=1,nz
	            apgridt = apgridt + w_ms*( &
	                              (1-wK)*(1-wmu)*policies%apgrid(:,:,zc,:,iK  ,imu  ) + &
	                                 wK *(1-wmu)*policies%apgrid(:,:,zc,:,iK+1,imu  ) + &
	                              (1-wK)*   wmu *policies%apgrid(:,:,zc,:,iK  ,imu+1) + &
	                                 wK *   wmu *policies%apgrid(:,:,zc,:,iK+1,imu+1) )

	            stockst = stockst + w_ms*( &
	                              (1-wK)*(1-wmu)*policies%stocks(:,:,zc,:,iK  ,imu  ) + &
	                                 wK *(1-wmu)*policies%stocks(:,:,zc,:,iK+1,imu  ) + &
	                              (1-wK)*   wmu *policies%stocks(:,:,zc,:,iK  ,imu+1) + &
	                                 wK *   wmu *policies%stocks(:,:,zc,:,iK+1,imu+1) )
			enddo
        else
            do zc=1,nz
                apgridt = apgridt + w_ms*((1-wK)*policies%apgrid(:,:,zc,:,iK,1) + wK*policies%apgrid(:,:,zc,:,iK+1,1))
                stockst = stockst + w_ms*((1-wK)*policies%stocks(:,:,zc,:,iK,1) + wK*policies%stocks(:,:,zc,:,iK+1,1))
            enddo
        endif

        rf1     = f_riskfree_rate(simvars%K(1),mu0,stat_dist_z)

    end subroutine get_initial_values

end subroutine simulate
!-------------------------------------------------------------------------------

subroutine print_error_msg(simvars)
    use types      ,only: tSimvars
    use kinds      ,only: dp
    use params_mod ,only: nt, t_scrap
    type(tSimvars) ,intent(in) :: simvars(:)
    integer                    :: count_err, i, n
    real(dp)                   :: perc_err

    n = size(simvars)
    if (any([(simvars(i)%err_K ,i=1,n)])) then
        count_err = count([(simvars(i)%err_K(t_scrap+1:) ,i=1,n)])
        perc_err  = count_err/real((nt-t_scrap)*n,dp)*100_dp
        print 214, ' Warning: simulate_economy: # K  not in grid =', count_err,'  (',perc_err,'%)'
    endif

    if (any([(simvars(i)%err_mu ,i=1,n)])) then
        count_err = count([(simvars(i)%err_mu(t_scrap+1:) ,i=1,n)])
        perc_err  = real(count_err,dp)/real((nt-t_scrap)*n,dp)*100_dp
        print 214, ' Warning: simulate_economy: # mu not in grid =', count_err,'  (',perc_err,'%)'
    endif

214 format((a,i6,a3,f5.1,a2))

end subroutine print_error_msg

end module simulation_mod
