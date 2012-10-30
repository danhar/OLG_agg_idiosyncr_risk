module policyfunctions
    use kinds
    implicit none
    private

    public tPolicies

    type tPolicies
        real(dp), allocatable, dimension(:,:,:,:,:,:) :: apgrid, kappa, stocks, xgrid ! policies /grids
    contains
        procedure :: solve => solve_policyfunctions
        procedure :: interpolate_tomorrow => interp_policies
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
! - pure subroutine interp_policies(kp, mup, grid, p, consp, xgridp, vp, app_min)
! - pure subroutine allocate_policies(p,nk,nmu)
! - pure subroutine deallocate_policies(p)
! - pure function mean(this,dimension_o,weight_o) result(mean_policy)
! - pure function interpolate(this,dim_x, gridx, x) result(pol_int)
! - pure subroutine calc_kappa(p)
!-------------------------------------------------------------------------------

subroutine solve_policyfunctions(p, coeffs, grids, value, err_o)
! Get the policy functions for the entire state space, i.e. both individual and aggregate states
    use params_mod      ,only: nj, nx, n_eta, nz, jr,surv, pi_z, pi_eta, cmin, g, beta, theta, gamm, apmax
    use household_solution_mod
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
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(p,value,cons,grids,coeffs,err_o,nmu,nk,nz,nj,n_eta,nx,beta,g,theta,gamm,surv, pi_eta, pi_z) &
!$OMP PRIVATE(jc,muc,kc,zc,betatildej,kp,mup,rp,rfp,yp,err_kp,err_mup,err_rfp,consp,xgridp,vp,app_min,err_asset,evp,err_cons)
jloop:do jc= nj-1,1,-1
        betatildej = beta*surv(jc)*(1.0+g)**((1.0-theta)/gamm)
!$OMP DO SCHEDULE(STATIC)
muloop: do muc=1,nmu
kloop:      do kc=1,nk
zloop: 			do zc=1,nz

                    call calc_vars_tomorrow(coeffs,grids,jc,zc,kc,muc,kp,mup,rp,rfp,yp,err_kp,err_mup,err_rfp)
                    if (present(err_o) .and. (err_kp .or. err_mup .or. err_rfp)) then
                        err_o%kp (zc,kc,muc) = err_kp
                        err_o%mup(zc,kc,muc) = err_mup
                        err_o%rfp(zc,kc,muc) = err_rfp
                    endif

                    call p%interpolate_tomorrow(cons, value, kp, mup, grids, jc, consp, xgridp, vp, app_min)

etaloop:            do ec=1,n_eta
	                    ! Create savings grid apgrid
	                    p%apgrid(:,ec,zc,jc,kc,muc)= f_apgrid_j(rfp,yp, xgridp, app_min)

xloop:                  do xc=1,nx
					        ! Asset allocation problem, see internal subroutine below
							call asset_allocation(xgridp, consp, vp, yp, rfp, rp, p%apgrid(xc,ec,zc,jc,kc,muc), pi_z(zc,:), pi_eta(ec,:), xc, jc, p%kappa(xc,ec,zc,jc,kc,muc), err_asset)
							p%stocks(xc,ec,zc,jc,kc,muc) = p%apgrid(xc,ec,zc,jc,kc,muc) * p%kappa(xc,ec,zc,jc,kc,muc)

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
!$OMP END DO
	enddo jloop
!$OMP END PARALLEL

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - pure function f_apgrid_j(kp, mup, rfp,yp, xgridp)
! - pure subroutine consumption(ap, kappa, cons_out, evp, error)
! - (pure) subroutine error_handling(err, err_asset, err_cons)
!-------------------------------------------------------------------------------

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

pure subroutine interp_policies(p,cons,value,kp, mup, grid, jc, consp, xgridp, vp, app_min)
! Make a projection of tomorrow's policies on k prime and mu prime
    use laws_of_motion  ,only: tCoeffs
    use aggregate_grids ,only: tAggGrids
    use fun_locate

    class(tPolicies) ,intent(in) :: p
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

end subroutine interp_policies
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
