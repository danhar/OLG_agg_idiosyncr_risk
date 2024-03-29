module fun_aggregate_diff
    use kinds,  only: dp
    implicit none

contains

    pure real(dp) function f_aggregate_diff(output,  invest, C,bequests)
        use params_mod ,only: bequests_to_newborn
        real(dp), intent(in) :: output, invest, C, bequests

        if (bequests_to_newborn) then
            f_aggregate_diff = (output - C - invest)/output
        else
            f_aggregate_diff = (output - C - invest - bequests)/output
        endif

    end function f_aggregate_diff

    pure function f_bequests(rf, r, stocks, ap, Phi)
        use params_mod, only: L_N_ratio, surv, g, n
        real(dp)                               :: f_bequests
        real(dp)                   ,intent(in) :: rf, r
        real(dp) ,dimension(:,:,:) ,intent(in) :: ap, stocks, Phi
        real(dp) ,dimension(:,:,:) ,allocatable:: assets_cum !asset cum returns
        integer :: nj
        nj = size(ap,3)

        assets_cum = ap*(1.0+rf)+stocks*(r-rf)
        f_bequests = sum((1.0-surv(2:))* sum(sum(assets_cum(:,:,:nj-1)*Phi(:,:,:nj-1),1),1))/(L_N_ratio*(1.0+g)*(1.0+n))

    end function f_bequests

    pure real(dp) function f_income_diff(K, zeta, r, rf, delta)
        ! Euler's formula
        use params_mod, only: alpha, de_ratio
        use income, only: f_grosswage
        real(dp), intent(in) :: K, zeta, r, rf, delta
        real(dp) :: output, total_factor_income

        output   = zeta*K**alpha
        total_factor_income = (r/(1.0 + de_ratio) + rf*de_ratio/(1.0 + de_ratio) + delta)*K + f_grosswage(K, zeta) !  simvars%wage(1) is after tau
        f_income_diff = (output - total_factor_income)/output
    end function f_income_diff

    ! ************************************************************
    ! Maybe later also implement a pension check, but now senseless.
    ! ************************************************************


end module fun_aggregate_diff
