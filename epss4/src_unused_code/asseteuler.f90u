module asseteuler
    use kinds
    use params_mod, only: nx, n_eta, nz
    implicit none

    private
    public asseteuler_f, asseteuler_set
    real(dp) ,dimension(:,:,:), allocatable :: xgridp, consp, vp  ! xgridp, cons, vp, for j+1
    real(dp), allocatable :: yp(:,:), rp(:)                       ! Tomorrow's income, risky return
    real(dp), allocatable :: pi_zp(:), pi_etap(:)                  ! Shock probabilities for given z, eta
    real(dp) :: ap, rfp                                    ! savings tomorrow, risk-free rate prime

    interface asseteuler_set
        module procedure asseteuler_set_most, asseteuler_set_last
    end interface asseteuler_set
contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure function asseteuler_f(kappa)
! - subroutine asseteuler_set_most(ypi, rpi, rfpi, pi_zi, pi_etai, conspi, xgridpi, vpi)
! - subroutine asseteuler_set_last(api)
!-------------------------------------------------------------------------------

pure function asseteuler_f(kappa)
    use params_mod         ,only: theta, gamm, g, cmin
	use fun_lininterp

    real(dp)                  :: asseteuler_f
    real(dp) ,intent(in)	  :: kappa
    real(dp) :: aeez(nz)                 ! asset euler equation for each z
    real(dp) ,dimension(1)    :: cons_interp, vp_interp, xp, aeetemp
    real(dp) :: rtildep          ! rtilde prime
    integer			          :: zpc, epc

    aeez = 0.0
    do zpc=1,nz
	    rtildep		= (1.0+rfp+kappa*(rp(zpc)-rfp))/(1.0+g)
	    do epc=1,n_eta
		    xp			= yp(epc,zpc)+rtildep*ap
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
!			if (xp(1) < xgridp(1,epc,zpc)) then
!				aeetemp = taylor_expansion(cons_interp, vp_interp)	! only works for theta\=1
!			else
!				aeetemp = vp_interp**((1.0-theta)*(gamm-1.0)/gamm)*cons_interp**((1.0-theta-gamm)/gamm)
!			endif

			aeez(zpc) = aeez(zpc) + pi_etap(epc)*aeetemp(1)
        enddo
    enddo

    asseteuler_f	= dot_product(pi_zp,aeez*(rp-rfp))

	contains
	!-------------------------------------------------------------------------------
	! Internal procedures in order:
	! - pure function taylor_expansion(cons_in, vp_in)
	!-------------------------------------------------------------------------------

		pure function taylor_expansion(cons_in, vp_in)
		! only works for theta\=1
			real(dp), dimension(:), intent(in) :: cons_in, vp_in
			real(dp), dimension(size(vp_in))	:: taylor_expansion

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

end function asseteuler_f

subroutine asseteuler_set_most(ypi, rpi, rfpi, pi_zi, pi_etai, conspi, xgridpi, vpi)
    real(dp) ,intent(in) :: ypi(:,:), rpi(:), rfpi, pi_zi(:), pi_etai(:)
    real(dp) ,dimension(nx,n_eta,nz) ,intent(in) :: xgridpi, conspi, vpi
    yp     = ypi
    rp     = rpi
    rfp    = rfpi
    pi_zp  = pi_zi
    pi_etap= pi_etai
    xgridp = xgridpi
    consp  = conspi
    vp     = vpi
end subroutine asseteuler_set_most

subroutine asseteuler_set_last(api)
    real(dp), intent(in) :: api
    ap    = api
end subroutine asseteuler_set_last

end module asseteuler
