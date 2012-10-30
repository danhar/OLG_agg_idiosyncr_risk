module asset_allocation_mod
    implicit none

contains

    !---------------------------------------------------------------------------
    ! Asset Allocation Problem
    !---------------------------------------------------------------------------
    pure subroutine asset_allocation(xgridp, consp, vp, yp, rfp, rp, ap, pi_zp, pi_etap, xc, jc, kappa_out, error)
        use kinds
        use params_mod ,only: opt_zbren, tol_asset_eul, opt_zbrak, kappa_in_01, scale_AR, de_ratio, g
!       use zreal_int      ! IMSL Math.pdf, p. 1195f: Mullers Method to find roots (like secant but quadratic).
                           ! expects array as input (i.e. define scalar as dimension(1)).
!       use zbren_int      ! IMSL Math.pdf, p. 1192f: Brent's method to find roots
        use fun_zbrent     ! NR: Brent method
        use sub_zbrac      ! NR: outwards bracketing
        use sub_zbrak      ! NR: inwards bracketing

        real(dp) ,dimension(:,:,:) ,intent(in) :: xgridp, consp, vp
        real(dp)   ,intent(in)  :: ap, rfp, yp(:,:), rp(:), pi_zp(:), pi_etap(:)        ! apgrid(xc,ec,zc,jc)
        integer   ,intent(in)   :: xc, jc
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
                !kappa_test(5) = p%kappa(xc,ec,zc,jc+1,kc,muc)
                !f_test(5) = asseteuler_f(p%kappa(xc,ec,zc,jc+1,kc,muc))

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
                kappa_out=kappa_start(1) !p%kappa(xc-1,ec,zc,jc,kc,muc)
            else
                kappa_out=kappa_result(1)
                !kappa(xc,zc,jc)=min(max(kappa_result(1),kappa1),kappa2)
                !kappa(xc,zc,jc)=min(max(kappa_result(1),-kappamax),kappamax)
            endif
        endif

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
        ! Taylor expansion of Euler euqation (optional)
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

end module asset_allocation_mod
