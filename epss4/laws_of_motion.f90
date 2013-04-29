module laws_of_motion
    use kinds      ,only: dp
    use classes_mod,only: tCoeffs
    use params_mod ,only: pooled_regression, loms_in_logs, lom_k_version, lom_mu_version

    implicit none
    private
    public Initialize, Regression, Forecast_k, Forecast_mu

    interface Initialize
        module procedure initialize_coeffs_ms, initialize_coeffs_guess, initialize_coeffs_simvars
    end interface Initialize

    type tMask
        logical, private, dimension(:), allocatable :: z, zp
    end type tMask

contains
!-------------------------------------------------------------------------------
! Module procedures in order: (would be nice in submodules)
! - pure real(dp) function Forecast(coeffs, k_in, mu_in_o)
! - pure subroutine Regression(simvars,coeffs)
! - function initialize_coeffs_ms() result(coeffs)
! - function initialize_coeffs_simvars(simvars) result(coeffs)
! - function initialize_coeffs_guess(agg_grid) result(coeffs)
!-------------------------------------------------------------------------------

pure real(dp) function Forecast_k(coeffs, k_in, mu)
    real(dp) ,intent(in)           :: coeffs(:), k_in, mu
    real(dp)                       :: k

    if (loms_in_logs) then
        k = log(k_in)
    else
        k = k_in
    endif
! Taking logs of mu didn't seemed counterproductive in tests

    select case (lom_k_version)
    case(1)
        Forecast_k = coeffs(1) + coeffs(2)*k
    case(2)
        Forecast_k = coeffs(1) + coeffs(2)*k + coeffs(3)*k**2
    case(3)
        Forecast_k = coeffs(1) + coeffs(2)*k + coeffs(3)*k**2 + coeffs(4)*mu
    case(4)
        Forecast_k = coeffs(1) + coeffs(2)*k + coeffs(3)*mu
    case(5)
        Forecast_k = coeffs(1) + coeffs(2)*k + coeffs(3)*mu + coeffs(4)*k*mu
    case(6)
        Forecast_k = coeffs(1) + coeffs(2)*k    + coeffs(3)*mu    &
                               + coeffs(4)*k**2 + coeffs(5)*mu**2
    case(7)
        Forecast_k = coeffs(1) + coeffs(2)*k    + coeffs(3)*mu    &
                               + coeffs(4)*k**2 + coeffs(5)*mu**2 + coeffs(6)*k*mu
    case default
        Forecast_k = coeffs(1) + coeffs(2)*k + coeffs(3)*k**2
    end select

    if (loms_in_logs) Forecast_k = exp(Forecast_k)

end function Forecast_k
!-------------------------------------------------------------------------------

pure real(dp) function Forecast_mu(coeffs, kp_in, mu)
! kp_in comes from Forecast_k(). One could also pass the derived type coeffs and call Forecast_k here when needed.
    real(dp) ,intent(in)           :: coeffs(:), kp_in, mu
    real(dp)                       :: kp

    if (loms_in_logs) then
        kp = log(kp_in)
    else
        kp = kp_in
    endif
! Taking logs of mu didn't seemed counterproductive in tests

    if (coeffs(1) == 0.0 .and. coeffs(2) == 1.0) then
        Forecast_mu = mu ! Handles mean shock case and STY calibration
        return
    endif

    select case (lom_mu_version)
    ! At the moment, cases structured like in Forecast_k, but could change it.
    case(1)
        Forecast_mu = coeffs(1) + coeffs(2)*kp
    case(2)
        Forecast_mu = coeffs(1) + coeffs(2)*kp + coeffs(3)*kp**2
    case(3)
        Forecast_mu = coeffs(1) + coeffs(2)*kp + coeffs(3)*kp**2 + coeffs(4)*mu
    case(4)
        Forecast_mu = coeffs(1) + coeffs(2)*kp + coeffs(3)*mu
    case(5)
        Forecast_mu = coeffs(1) + coeffs(2)*kp + coeffs(3)*mu + coeffs(4)*kp*mu
    case(6)
        Forecast_mu = coeffs(1) + coeffs(2)*kp    + coeffs(3)*mu    &
                                + coeffs(4)*kp**2 + coeffs(5)*mu**2
    case(7)
        Forecast_mu = coeffs(1) + coeffs(2)*kp    + coeffs(3)*mu    &
                                + coeffs(4)*kp**2 + coeffs(5)*mu**2 + coeffs(6)*kp*mu
    case default
        Forecast_mu = coeffs(1) + coeffs(2)*kp + coeffs(3)*kp**2
    end select

end function Forecast_mu
!-------------------------------------------------------------------------------

pure subroutine Regression(simvars,coeffs)
    use params_mod   ,only: t_scrap
    use classes_mod  ,only: tSimvars
    use MKL95_LAPACK ,only: gels, gelsy ! Intel MKL/ LAPACK dgels. Explanation at the bottom.

    type(tSimvars)           ,intent(in)   :: simvars(:)
    type(tCoeffs)            ,intent(inout):: coeffs   ! inout coz some coeffs might not be updated
    type(tMask)                            :: mask(size(simvars))
    real(dp) ,dimension(:,:) ,allocatable  :: cov_k, cov_mu  ! covariates, i.e. rhs
    real(dp) ,dimension(:) ,allocatable    :: k_z, kp_z, kp_zp, mu_z, mu_zp, mup_zp              ! lhs
    real(dp) ,dimension(2) :: residual_sum_of_squares, total_sum_of_squares
    real(dp)               :: mean
    integer                :: n_covars_k, n_covars_mu, nz, nt, zc, n_zc(size(simvars)), n_zpc(size(simvars)), info, i, shape(2)
    logical, parameter :: use_gelsy = .false. ! instead of using gels (which had an error in Intel Composer XE 2013 original release)

    n_covars_k  = size(coeffs%k,1)
    n_covars_mu = size(coeffs%mu,1)
    nz       = size(coeffs%k,2)
    nt       = size(simvars(1)%z)

    do zc=1,nz
        if (pooled_regression) then
            n_zc = nt-t_scrap
            n_zpc = n_zc
            do i=1,size(simvars)
                mask(i)%z  = .true. ! select all
                mask(i)%zp = .true. ! select all
            enddo
        else
            do i=1,size(simvars)
                n_zc(i)    = count(simvars(i)%z(t_scrap:nt-1)==zc)
                n_zpc(i)   = count(simvars(i)%z(t_scrap+1:nt)==zc)
                mask(i)%z  = [simvars(i)%z(t_scrap:nt-1)==zc]
                mask(i)%zp = [simvars(i)%z(t_scrap+1:nt)==zc]
            enddo
            if (sum(n_zc) < 2) cycle     ! coeffs stay the same. Takes care of STY calibration, or degenerate cases.
        endif

!        allocate(kp_z(sum(n_zc)), mup_zp(sum(n_zpc)))

        k_z    = [(pack(simvars(i)%K (t_scrap:nt-1), mask(i)%z ), i=1,size(simvars))]
        kp_z   = [(pack(simvars(i)%K (t_scrap+1:nt), mask(i)%z ), i=1,size(simvars))]
        kp_zp  = [(pack(simvars(i)%K (t_scrap+1:nt), mask(i)%zp), i=1,size(simvars))]
        mu_z   = [(pack(simvars(i)%mu(t_scrap:nt-1), mask(i)%z ), i=1,size(simvars))]
        mu_zp  = [(pack(simvars(i)%mu(t_scrap:nt-1), mask(i)%zp), i=1,size(simvars))] ! Because lom for mu conditions on shock tomorrow!
        mup_zp = [(pack(simvars(i)%mu(t_scrap+1:nt), mask(i)%zp), i=1,size(simvars))]

        if (loms_in_logs) then
            k_z   = log(k_z)
            kp_z  = log(kp_z)
            kp_zp = log(kp_zp)
        ! Taking log of mu seemed counterproductive
        endif

        mean = sum(kp_z) /real(sum(n_zc),dp)
        total_sum_of_squares(1) = dot_product((kp_z -mean), (kp_z -mean))
        mean = sum(mup_zp)/real(sum(n_zpc),dp)
        total_sum_of_squares(2) = dot_product((mup_zp-mean), (mup_zp-mean))

        allocate(cov_k(sum(n_zc),n_covars_k), cov_mu(sum(n_zpc),n_covars_mu))

        shape = [size(cov_k,1),n_covars_k-1]
        cov_k(:,1) = 1.0
        select case(lom_k_version)
        case(1)
            cov_k(:,2:) = reshape([  k_z                                   ],shape)
        case(2)
            cov_k(:,2:) = reshape([  k_z, k_z**2                           ],shape)
        case(3)
            cov_k(:,2:) = reshape([  k_z, k_z**2, mu_z                     ],shape)
        case(4)
            cov_k(:,2:) = reshape([  k_z, mu_z                             ],shape)
        case(5)
            cov_k(:,2:) = reshape([  k_z, mu_z, k_z*mu_z                   ],shape)
        case(6)
            cov_k(:,2:) = reshape([  k_z, mu_z, k_z**2, mu_z**2            ],shape)
        case(7)
            cov_k(:,2:) = reshape([  k_z, mu_z, k_z**2, mu_z**2, k_z*mu_z  ],shape)
        case default
            cov_k(:,2:) = reshape([  k_z, k_z**2                           ],shape)
        end select

        shape = [size(cov_mu,1),n_covars_mu-1]
        cov_mu(:,1) = 1.0
        select case(lom_mu_version)
        case(1)
            cov_mu(:,2:) = reshape([  kp_zp                                          ],shape)
        case(2)
            cov_mu(:,2:) = reshape([  kp_zp, kp_zp**2                                ],shape)
        case(3)
            cov_mu(:,2:) = reshape([  kp_zp, kp_zp**2, mu_zp                         ],shape)
        case(4)
            cov_mu(:,2:) = reshape([  kp_zp, mu_zp                                   ],shape)
        case(5)
            cov_mu(:,2:) = reshape([  kp_zp, mu_zp, kp_zp*mu_zp                      ],shape)
        case(6)
            cov_mu(:,2:) = reshape([  kp_zp, mu_zp, kp_zp**2, mu_zp**2               ],shape)
        case(7)
            cov_mu(:,2:) = reshape([  kp_zp, mu_zp, kp_zp**2, mu_zp**2, kp_zp*mu_zp  ],shape)
        case default
            cov_mu(:,2:) = reshape([  kp_zp, kp_zp**2                                ],shape)
        end select

        if (use_gelsy) then ! This helps to check results. Intel Composer XE 13.0 had a bug in gels
            call gelsy(cov_k,kp_z    ,info=info) ! should print or return info
            call gelsy(cov_mu,mup_zp ,info=info)
        else
            call gels(cov_k,kp_z     ,info=info) ! should print or return info
            call gels(cov_mu,mup_zp  ,info=info)
        endif
        coeffs%k(:,zc)  = kp_z  (1:n_covars_k )
        coeffs%mu(:,zc) = mup_zp(1:n_covars_mu)

        residual_sum_of_squares(1) = dot_product(kp_z  (n_covars_k +1:),kp_z  (n_covars_k +1:))
        residual_sum_of_squares(2) = dot_product(mup_zp(n_covars_mu+1:),mup_zp(n_covars_mu+1:))
        if (total_sum_of_squares(1) == 0.0) then
            coeffs%r_squared(1,zc) = 1.0
        else
            coeffs%r_squared(1,zc) = 1.0-residual_sum_of_squares(1)/total_sum_of_squares(1)
        endif
        if (total_sum_of_squares(2) == 0.0) then
            coeffs%r_squared(2,zc) = 1.0
        else
            coeffs%r_squared(2,zc) = 1.0-residual_sum_of_squares(2)/total_sum_of_squares(2)
        endif
        deallocate(cov_mu, cov_k)
    enddo
! ---------------------Usage of LAPACK / Intel MKL routine gels --------------------
! gels(a, b [,trans] [,info])
! Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
! 1. If trans = 'N' and m ≥ n: find the least squares solution of an overdetermined system, that is, solve the least squares problem
! minimize ||b - A*x||_2
! INPUT
! a: Holds the matrix A of size (m,n).
! b: Holds the matrix of size max(m,n)-by-nrhs.
!       If trans = 'N', then, on entry, the size of b is m-by-nrhs,
!       If trans = 'T', then, on entry, the size of b is n-by-nrhs,
! trans:  Must be 'N' or 'T'. The default value is 'N'.
!
! OUTPUT
! a: On exit, overwritten by the factorization data as follows:
!        if m ≥ n, array a contains the details of the QR factorization of the matrix A as returned by ?geqrf;
!        if m < n, array a contains the details of the LQ factorization of the matrix A as returned by ?gelqf.
! b: If info = 0, b overwritten by the solution vectors, stored columnwise:
!        if trans = 'N' and m ≥ n, rows 1 to n of b contain the least squares solution vectors; the residual sum of squares for the solution in each column is given by the sum of squares of modulus of elements n+1 to m in that column;
!        if trans = 'N' and m < n, rows 1 to n of b contain the minimum norm solution vectors;
!        if trans = 'T' or 'C' and m ≥ n, rows 1 to m of b contain the minimum norm solution vectors; if trans = 'T' or 'C' and m < n, rows 1 to m of b contain the least squares solution vectors; the residual sum of squares for the solution in each column is given by the sum of squares of modulus of elements m+1 to n in that column.
end subroutine Regression
!-------------------------------------------------------------------------------

pure function initialize_coeffs_ms() result(coeffs)
    ! Degenerate coefficients for the mean shock equilibrium
    ! Can't make this a type-bound procedure of coeffs because then circular dependency
    type(tCoeffs)   :: coeffs

    call coeffs%allocate()
    coeffs%r_squared = 0.0
    coeffs%k         = 0.0
    coeffs%mu        = 0.0
    coeffs%k(2,:)  = 1.0     ! same for loms_in_logs
    coeffs%mu(2,:) = 1.0

end function initialize_coeffs_ms
!-------------------------------------------------------------------------------

pure function initialize_coeffs_simvars(simvars_old) result(coeffs)
    ! Can't make this a type-bound procedure of coeffs because then circular dependency
    use classes_mod   ,only: tSimvars

    type(tCoeffs)              :: coeffs
    type(tSimvars) ,intent(in) :: simvars_old(:)

    call coeffs%allocate()
    call Regression(simvars_old,coeffs)

end function initialize_coeffs_simvars
!-------------------------------------------------------------------------------

function initialize_coeffs_guess(agg_grid) result(coeffs)
    ! This is for putting hard-coded guesses for the coefficients
    ! Usually only necessary for very first Krusell-Smith runs.
    ! Can't make this a type-bound procedure of coeffs because then circular dependency
    use params_mod    ,only: pi1_delta, scale_IR
    use classes_mod   ,only: tAggGrids

    type(tCoeffs)                :: coeffs
    type(tAggGrids)  ,intent(in) :: agg_grid
    integer :: n_coeffs_k, n_coeffs_mu, nz

    call coeffs%allocate()
    coeffs%r_squared = 0.0
    coeffs%k         = 0.0
    coeffs%mu        = 0.0
    n_coeffs_k  = size(coeffs%k ,1)
    n_coeffs_mu = size(coeffs%mu,1)
    nz          = size(coeffs%k ,2)
!    print*, 'laws_of_motion:initialize_coeffs: using hard-coded guess'
    if (n_coeffs_k .ne. n_coeffs_mu) stop 'CRITICAL STOP: n_coeffs_k .ne. n_coeffs_mu'

coef:   select case (n_coeffs_k)
        case(2)
pi:     if (pi1_delta==1.0) then

            coeffs%k = reshape([ 3.241278E-01,  8.915865E-01, &
                                 0.000000E+00,  1.000000E+00, &
                                 0.000000E+00,  1.000000E+00, &
                                 3.287817E-01,  9.876539E-01],&
                       shape(coeffs%k))
            coeffs%mu= reshape([-1.191016E-02,  1.946109E+00, &
                                 0.000000E+00,  1.000000E+00, &
                                 0.000000E+00,  1.000000E+00, &
                                 5.598800E-03,  6.062592E-01],&
                       shape(coeffs%mu))

        else pi ! here (pi1_delta/=1.0)
pooled:     if (pooled_regression) then
ir:             if ((1.0 + scale_IR) > 0.1_dp) then
                        coeffs%k = reshape([ 5.846835E-01,  8.596145E-01, &
                                             5.846835E-01,  8.596145E-01, &
                                             5.846835E-01,  8.596145E-01, &
                                             5.846835E-01,  8.596145E-01],&
                                   shape(coeffs%k))
                        coeffs%mu= reshape([ 1.495434E-02,  6.518498E-01, &
                                             1.495434E-02,  6.518498E-01, &
                                             1.495434E-02,  6.518498E-01, &
                                             1.495434E-02,  6.518498E-01],&
                                   shape(coeffs%mu))
                    else
                        coeffs%k = reshape([ 1.691731E-01,  9.287697E-01, &
                                             1.691731E-01,  9.287697E-01, &
                                             1.691731E-01,  9.287697E-01, &
                                             1.691731E-01,  9.287697E-01],&
                                   shape(coeffs%k))
                        coeffs%mu= reshape([ 2.271671E-02,  4.942666E-01, &
                                             2.271671E-02,  4.942666E-01, &
                                             2.271671E-02,  4.942666E-01, &
                                             2.271671E-02,  4.942666E-01],&
                                   shape(coeffs%mu))
                    endif ir
                else pooled
pi2:                if (pi1_delta==0.5_dp) then
ir2:                    if ((1.0 + scale_IR) > 0.1_dp) then
                            coeffs%k = reshape([ 2.813136E-01,  7.990639E-01, &
                                                 6.045659E-01,  9.972293E-01, &
                                                 3.854341E-01,  7.827162E-01, &
                                                 5.677132E-01,  1.014907E+00],&
                                       shape(coeffs%k))
                            coeffs%mu= reshape([ 3.153609E-03,  9.444029E-01, &
                                                 6.153475E-03,  8.864563E-01, &
                                                 3.287029E-03,  9.491282E-01, &
                                                 7.719321E-03,  8.651943E-01],&
                                       shape(coeffs%mu))
                        else
                            coeffs%k = reshape([ 1.248526E-01,  8.248078E-01, &
                                                 2.138598E-01,  1.013803E+00, &
                                                 2.198917E-01,  7.827063E-01, &
                                                 2.327648E-01,  1.027819E+00],&
                                       shape(coeffs%k))
                            coeffs%mu= reshape([ 2.456737E-03,  9.265814E-01, &
                                                 7.047098E-04,  9.847655E-01, &
                                                 5.484127E-03,  8.878176E-01, &
                                                 1.032811E-02,  7.893054E-01],&
                                       shape(coeffs%mu))
                        endif ir2
                    else pi2 ! here (pi1_delta/=0.5_dp)
! Zurich
                        coeffs%k = reshape([ 3.804995E-01,  7.930654E-01, &
                                             9.636235E-01,  9.127569E-01, &
                                             3.680785E-01,  8.124554E-01, &
                                             1.194160E+00,  8.838562E-01],&
                                   shape(coeffs%k))
                        coeffs%mu= reshape([ 7.968296E-03,  8.921678E-01, &
                                             3.147255E-03,  9.670157E-01, &
                                             1.195201E-02,  7.841391E-01, &
                                             1.019259E-02,  8.145751E-01],&
                                   shape(coeffs%mu))
! neta5, psi15, tau12
!                            coeffs%k = reshape([ 4.630577E-01,  7.567895E-01, &
!                                                 7.327892E-01,  9.549323E-01, &
!                                                 5.834030E-01,  7.626787E-01, &
!                                                 9.587330E-01,  9.409641E-01],&
!                                       shape(coeffs%k))
!                            coeffs%mu= reshape([ 8.800742E-03,  8.832917E-01, &
!                                                 6.087253E-03,  9.195544E-01, &
!                                                 1.296528E-02,  7.778659E-01, &
!                                                 1.166864E-02,  7.996547E-01],&
!                                       shape(coeffs%mu))
                    endif pi2
                endif pooled
            endif pi

            if (loms_in_logs) then
                coeffs%k(1,:)  = log(agg_grid%k(1))*(1-coeffs%k(2,:))
                coeffs%mu(1,:) = log(agg_grid%mu(1))*(1-coeffs%mu(2,:))
            endif

        case(3:) coef ! all higher 3
            if (pooled_regression) then

                coeffs%k(2,:)  = 9.8E-01
                coeffs%k(3,:)  = 2.1e+02
                coeffs%mu(2,:) = 2.2e-05
                coeffs%mu(3,:) = 3.0E-01

                if (loms_in_logs) then
                    coeffs%k(1,:)  = 1.770695
                    coeffs%mu(1,:) = -7.241165E+01
                else
                    coeffs%k(1,:)  = -3.0
                    coeffs%mu(1,:) = 1.1E-02
                endif

            elseif (pi1_delta == 1.0) then ! STY
                if (loms_in_logs) then ! not implemented
                    coeffs%k(1,:)  = 1.770695 !log(agg_grid%k(1))*(1-coeffs%k(2,:))
                    coeffs%mu(1,:) = -7.241165E+01 !log(agg_grid%mu(1))*(1-coeffs%mu(3,:))
                else
                    coeffs%k = reshape([ 1.598852E-01,  8.942185E-01,  1.191438E+01, &
                                         0.000000E+00,  1.000000E+00,  0.000000E+00, &
                                         0.000000E+00,  1.000000E+00,  0.000000E+00, &
                                         1.090199E+00,  9.811942E-01, -4.155187E+01],&
                               shape(coeffs%k))
                    coeffs%mu= reshape([ 2.801304E-03, -1.496703E-04,  9.485543E-01, &
                                         0.000000E+00,  0.000000E+00,  1.000000E+00, &
                                         0.000000E+00,  0.000000E+00,  1.000000E+00, &
                                         5.306034E-03,  1.330638E-05,  6.162367E-01],&
                               shape(coeffs%mu))
                endif
            else ! here (pi1_delta /= 1.0)
                if (loms_in_logs) then
                    ! This is for calibration calib_ep/STY2cc
                    coeffs%k = reshape([ 2.200245E-01,  8.084008E-01,  2.839362E-02, &
                                         1.924221E-01,  9.414862E-01,  8.135341E-04, &
                                         2.258880E-01,  8.114570E-01,  2.707929E-02, &
                                         1.893091E-01,  9.499170E-01, -1.361134E-03],&
                               shape(coeffs%k))
                    coeffs%mu= reshape([ 2.465450E-02, -7.913885E-03,  1.833304E-03, &
                                         3.133807E-02, -1.456047E-02,  3.452852E-03, &
                                         2.630511E-02, -8.020631E-03,  1.855730E-03, &
                                         3.444827E-02, -1.608870E-02,  3.809692E-03],&
                               shape(coeffs%mu))
                else    ! not implemented
                    coeffs%k(1,:)  = 1.770695 !log(agg_grid%k(1))*(1-coeffs%k(2,:))
                    coeffs%mu(1,:) = -7.241165E+01 !log(agg_grid%mu(1))*(1-coeffs%mu(3,:))
                endif
            endif
        end select coef
 end function initialize_coeffs_guess
 !-------------------------------------------------------------------------------

end module laws_of_motion
