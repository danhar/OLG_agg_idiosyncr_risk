! Copyright (C) 2016 Daniel Harenberg - All Rights Reserved
module income
! contains functions for net wages, net returns, and pensions
    use kinds
    use params_mod,     only: alpha, tau, zeta, delta
    implicit none

contains

    elemental real(dp) function f_net_mpk(k,zeta, delta)
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

end module income
