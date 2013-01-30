module income
! contains functions for net wages, net returns, and pensions
    use kinds
    use params_mod,     only: alpha, tau, zeta, delta
    implicit none

contains

    pure real(dp) function f_net_mpk(k,zeta, delta)
        real(dp),intent(in) :: k, zeta, delta     ! k is capital per capita
        f_net_mpk= zeta*alpha*k**(alpha-1.0) - delta
    end function f_net_mpk

    pure real(dp) function f_expected_return(k,cond_prob) result(Erp)
        real(dp),intent(in) :: k, cond_prob(:)
        integer :: zpc
        Erp = 0.0
        do zpc = 1,size(cond_prob)
            Erp = Erp + cond_prob(zpc)* f_net_mpk(k,zeta(zpc),delta(zpc))
        enddo
    end function f_expected_return

    pure real(dp) function f_riskfree_rate(k,mu,cond_prob) result(rf)
        use params_mod, only: de_ratio
        real(dp),intent(in) :: k, mu, cond_prob(:)
        real(dp) :: Erp
        Erp = f_expected_return(k,cond_prob)
        rf  = Erp - mu/(1.0 + de_ratio)
    end function f_riskfree_rate

    pure real(dp) function f_stock_return(k,zeta, delta,rf) result(r)
        use params_mod, only: de_ratio
        real(dp),intent(in) :: k, zeta, delta, rf     ! k is capital per capita
        r = f_net_mpk(k,zeta, delta)*(1.0+de_ratio) - de_ratio * rf
    end function f_stock_return

    pure function f_netwage(k,zeta)
        real(dp),intent(in) :: k, zeta   ! k is capital per capita
        real(dp)            :: f_netwage
        f_netwage = f_grosswage(k,zeta)*(1.0-f_tau(k,zeta))
    end function f_netwage

    pure function f_grosswage(k,zeta)
        real(dp),intent(in) :: k, zeta
        real(dp)            :: f_grosswage
        f_grosswage = zeta*(1.0-alpha)*k**(alpha)
    end function f_grosswage

    pure function f_pensions(k,zeta)
        use params_mod, only : P_L_ratio
        real(dp),intent(in) :: k, zeta
        real(dp)            :: f_pensions
        f_pensions=f_grosswage(k,zeta)*f_tau(k,zeta)/P_L_ratio
    end function f_pensions

    pure function f_tau(k,zeta)
        use params_mod, only : def_contrib, P_L_ratio, def_benefits
        real(dp)            :: f_tau
        real(dp),intent(in) :: k, zeta
        if (def_contrib) then
            f_tau = tau
        else
            f_tau = def_benefits/f_grosswage(k,zeta)*P_L_ratio
        endif
    end function f_tau

!-------------------------------------------------------------------------------
! The next subroutine should be in aggregate_grid, but could not due to circular dependency: params_mod-> income -> aggregate_grid -> params_mod
! So I put it here although it has nothing to do here. It's written in a type-bound procedure way.
    subroutine calc_max_mu(this)
        use classes_mod ,only : tAggGrids
        use params_mod, only: pi_z
        use fun_zbrent     ! NR: Brent method
        use sub_zbrac      ! NR: outwards bracketing
        use sub_zbrak      ! NR: inwards bracketing
        use sub_broyden

        class(tAggGrids), intent(inout)  :: this
        real(dp) ,parameter :: tolerance = 1e-10_dp
        logical  ,parameter :: use_broyden =.false.
        real(dp) :: max_mu_temp(1), fvals(1), brack1, brack2
        real(dp) ,dimension(:) ,allocatable :: kappal, kappau !sub_zbrak: lower and upper bounds of segment containing a root
        logical  :: not_converged, bracket_found

        if (.not.(this%fixed .and. this%infer_max_mu)) return

        if (use_broyden) then
            max_mu_temp(1) = this%max_mu
            call s_broyden(return_diff2,max_mu_temp,fvals,not_converged, get_fd_jac_o=.true.,tolf_o=tolerance)
            if (.not. not_converged) this%max_mu = max_mu_temp(1)

        else ! use Brent's algorithm
            bracket_found = .false.
            brack1 = this%min_mu
            brack2 = this%max_mu
            if (return_diff(brack1) * return_diff(brack2) <0.0 ) then
                bracket_found = .true.
            else
                call s_zbrak(return_diff,brack1,brack2,1000,kappal,kappau)
                if (size(kappal) > 0) then
                    bracket_found=.true.
                    brack1    = kappal(1)
                    brack2    = kappau(1)
                    deallocate(kappal,kappau)
                endif
                if (.not. bracket_found) then
                    ! Look 'outwards', i.e. extend the bounds
                    call s_zbrac(return_diff,brack1,brack2,bracket_found)
                endif
            endif

            if (bracket_found) this%max_mu= f_zbrent(return_diff,brack1,brack2,tolerance)
        endif

        contains

            pure real(dp) function return_diff(mu)
                real(dp) ,intent(in) :: mu
                real(dp) :: k, rf
                real(dp), parameter :: rf_above_rp1 = 1e-2_dp ! This is important

                k = this%max_k
                rf = f_riskfree_rate(k,mu,pi_z(1,:))
                return_diff = rf - f_stock_return(k,zeta(1), delta(1),rf) - rf_above_rp1
            end function return_diff

            function return_diff2(mu)
                real(dp) , dimension(1) :: return_diff2
                real(dp) ,intent(in) :: mu(1)
                real(dp) :: k, rf
                real(dp), parameter :: rf_above_rp1 = 1e-2_dp ! This is important

                k = this%max_k
                rf = f_riskfree_rate(k,mu(1),pi_z(1,:))
                return_diff2 = [rf - f_stock_return(k,zeta(1), delta(1),rf) - rf_above_rp1]
            end function return_diff2

    end subroutine calc_max_mu


end module income
