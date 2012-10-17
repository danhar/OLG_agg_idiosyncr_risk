module distribution
    use kinds
    use params_mod, only: nj, nx, n_eta, n, jr, g, pi_eta, ej, pop_frac, surv, stat_dist_eta
    implicit none

    private check_Phi_full, check_Phi_small, check_Phi_pure
    interface CheckPhi
        module procedure check_Phi_full, check_Phi_small, check_Phi_pure
    end interface

contains
!-------------------------------------------------------------------------------
! - pure function TransitionPhi(rf,r,netwage,pens,xgrid,apgrid,kappa,etagrid, Phitm_o) result(Phi)
! - pure subroutine check_Phi_pure(Phi_full, bound_Phi_1, bound_Phi_2)
! - subroutine check_Phi_full(Phi_full,dir)
! - subroutine check_Phi_small(Phi,dir)
!-------------------------------------------------------------------------------

pure function TransitionPhi(rf,r,netwage,pens,xgrid,apgrid,kappa,etagrid, Phitm_o) result(Phi)
! Calculate the transition of Phi
    use fun_locate

    real(dp) ,dimension(nx,n_eta,nj) ,target     :: Phi                  ! new distribution
    real(dp) ,dimension(nx,n_eta,nj) ,intent(in) ,optional:: Phitm_o     ! distribution from t-1 (t minus)
    real(dp)                         ,intent(in) :: rf, r, netwage, pens ! today's risk-free rate, risky return, net wage, pensions
    real(dp) ,dimension(nx,n_eta,nj) ,intent(in) :: apgrid, kappa, xgrid ! optimal policies/ grids at today's aggregate state
    real(dp) ,dimension(n_eta)       ,intent(in) :: etagrid              ! idiosyncratic income shocks today
    real(dp) ,dimension(:,:,:)       ,pointer    :: Phitm                ! distribution previous period/ generation (Phi 'T M'inus one)
    real(dp) ,dimension(size(etagrid)) :: y ! income
    real(dp) :: wx, x, rtilde ! weight, cash at hand, weighted return
    integer  :: ix, jc, ec, xmc, emc   ! xmc,emc: x/eta previous period or generation

    if (present(Phitm_o)) then
        allocate(Phitm(nx,n_eta,nj))
        Phitm = Phitm_o
    else
        Phitm => Phi
    endif

    Phi         = 0.0
    ! Generation jc=1
    y           = netwage*ej(1)*etagrid
    do ec = 1,n_eta ! this loop is necessary because several ix could have same value, or +/- 1, so that Phi would be overwritten
	    x           = y(ec)             ! agents born with zero assets
	    ix          = f_locate(xgrid(:,ec,1),x)
	    wx          = (x-xgrid(ix,ec,1))/(xgrid(ix+1,ec,1)-xgrid(ix,ec,1))
	    wx          = max(min(wx,1.0),0.0)
        Phi(ix  ,ec,1) = Phi(ix  ,ec,1) + stat_dist_eta(ec) * (1.0 - wx) * pop_frac(1)
        Phi(ix+1,ec,1) = Phi(ix+1,ec,1) + stat_dist_eta(ec) *        wx  * pop_frac(1)
    enddo

    ! Generations jc=2 to nj
    do jc=2,nj
        if (jc>=jr) then
            y=pens
        else
            y = netwage*ej(jc)*etagrid
        endif
        do emc= 1,n_eta
            do ec = 1,n_eta
                do xmc = 1, nx ! this loop is necessary because several ix could have same value, or +/- 1, so that Phi would be overwritten
		            if (Phitm(xmc,emc,jc-1)==0.0) cycle
		            rtilde      = (1.0+rf+kappa(xmc,emc,jc-1)*(r-rf))/(1.0+g)
		            x           = y(ec)+rtilde*apgrid(xmc,emc,jc-1)
		            ix          = f_locate(xgrid(:,ec, jc),x)
		            wx          = (x-xgrid(ix,ec,jc))/(xgrid(ix+1,ec,jc)-xgrid(ix,ec,jc))
		            wx          = max(min(wx,1.0),0.0)
	                Phi(ix  ,ec,jc) = Phi(ix  ,ec,jc) + pi_eta(emc, ec) * (1.0 - wx)*Phitm(xmc,emc,jc-1)/(1.0 + n) * surv(jc)
	                Phi(ix+1,ec,jc) = Phi(ix+1,ec,jc) + pi_eta(emc, ec) *        wx *Phitm(xmc,emc,jc-1)/(1.0 + n) * surv(jc)
	            enddo
            enddo
        enddo
    enddo

    if (present(Phitm_o)) then
        deallocate(Phitm)
    else
        nullify(Phitm)
    endif

end function TransitionPhi
!-------------------------------------------------------------------------------

pure subroutine check_Phi_pure(Phi_full, bound_Phi_1, bound_Phi_2)
! Phi including eta
    real(dp), dimension(nx,n_eta,nj), intent(in) :: Phi_full
    real(dp), intent(out)       :: bound_Phi_1, bound_Phi_2
    real(dp), dimension(nx,nj)  :: Phi

    Phi = sum(Phi_full,2)
    bound_Phi_1=sum(Phi(1,:))
    bound_Phi_2=sum(Phi(nx,:))

end subroutine check_Phi_pure
!-------------------------------------------------------------------------------

subroutine check_Phi_full(Phi_full,path)
! Phi including eta
    real(dp), dimension(nx,n_eta,nj), intent(in) :: Phi_full
    character(len=*), intent(in)        :: path  ! directory to write to
    real(dp), dimension(nx,nj)  :: Phi

    Phi = sum(Phi_full,2)

    call check_Phi_small(Phi,path)

end subroutine check_Phi_full
!-------------------------------------------------------------------------------

subroutine check_Phi_small(Phi,path)
! Phi excluding eta
    real(dp), dimension(nx,nj), intent(in) :: Phi
    character(len=*), intent(in)        :: path  ! directory to write to
    real(dp)      :: marg_Phi(nj), bounds_Phi(2), mass_Phi, crit

    crit= 1e-8_dp
    marg_Phi=sum(Phi,dim=1)
    bounds_Phi(1)=sum(Phi(1,:))
    bounds_Phi(2)=sum(Phi(nx,:))
    mass_Phi = sum(Phi)

    if (any(abs(marg_Phi- pop_frac)>crit) .or. any(bounds_Phi>crit) .or. abs(mass_Phi-1.0)>crit) then
        print*, 'Warning: check_Phi: marginals too large, see err_distribution.txt'
        open(unit=301, file=path//'/err_distribution.txt', status = 'replace')
        write (301,*) 'min(marg_Phi-pop_frac) = ', minval(marg_Phi-pop_frac)
        write (301,*) 'max(marg_Phi-pop_frac) = ', maxval(marg_Phi-pop_frac)
        write (301,*) 'sum(marg_Phi)        = ', sum(marg_Phi)
        write (301,*) '-------------------------'
        write (301,*) 'sum(Phi(1,:)) = ', bounds_Phi(1)
        write (301,*) 'sum(Phi(nx,:))= ', bounds_Phi(2)
        write (301,*) 'sum(Phi)      = ', mass_Phi
        close(301)
    endif
end subroutine check_Phi_small

end module distribution
