module laws_of_motion
    use kinds      ,only: dp
    use params_mod ,only: n_coeffs, nz, pooled_regression, loms_in_logs, pi1_delta, scale_IR

    implicit none
    private
    public tCoeffs, MakeType, MakeVector, Initialize, Regression, Forecast
    real(dp), parameter :: norm_factor_small = 10.0, norm_factor_large = 100.0 ! should be type component. Could normalize everytime with a different value (e.g. I could normalize all coeffs = 1.0 every time in the krusell-smith alg)

    type tCoeffs
        real(dp), dimension(:,:), allocatable :: k, mu, r_squared  ! coefficients for laws of motion, and residual sum of squares
    end type tCoeffs

    interface MakeType
        module procedure maketype_coeffs
    end interface

    interface MakeVector
        module procedure makevector_coeffs
    end interface

    interface Initialize
        module procedure initialize_coeffs
    end interface

contains
!-------------------------------------------------------------------------------
! Module procedures in order: (would be nice in submodules)
! - pure function makevector_coeffs(coeffs_in,o_norm) result(coeff_vec)
! - pure function maketype_coeffs(coeff_vec, o_norm) result(coeffs)
! - function initialize_coeffs(dir,n_coeffs, agg_grid_o) result(coeffs)
! - pure subroutine Regression(simvars,coeffs)
! - pure real(dp) function Forecast(coeffs, k_in, mu_in_o)
!-------------------------------------------------------------------------------

    pure function makevector_coeffs(coeffs_in,normalize) result(coeff_vec)
        real(dp), dimension(:), allocatable :: coeff_vec
        type(tCoeffs) ,intent(in) :: coeffs_in
        logical       ,intent(in) :: normalize
        type(tCoeffs) :: coeffs     ! need this only because coeffs_in has intent(in), and if I want to normalize
        integer :: zc, n_coeffs

        n_coeffs = size(coeffs_in%k,1)
        coeffs = coeffs_in

        if (normalize) then
            coeffs%mu(1,:) = coeffs%mu(1,:)*norm_factor_small
            coeffs%mu(2,:) = coeffs%mu(2,:)*norm_factor_large
            if (n_coeffs >= 3) then
                coeffs%k(3,:) = coeffs%k(3,:)*norm_factor_small
                coeffs%mu(3,:) = coeffs%mu(3,:)*norm_factor_large
            endif
        endif

        if (pooled_regression) then
            allocate(coeff_vec(2*n_coeffs))
            coeff_vec(1:n_coeffs) = coeffs%k(:,1)
            coeff_vec(n_coeffs+1:)= coeffs%mu(:,1)
        elseif (pi1_delta==1.0) then ! Alternatively, could have nz_1=count(stat_dist_z==0.0)
            allocate(coeff_vec(2*n_coeffs*nz/2))
            coeff_vec(1:n_coeffs) = coeffs%k(:,1)
            coeff_vec(n_coeffs+1:2*n_coeffs)= coeffs%k(:,4)
            coeff_vec(2*n_coeffs+1:3*n_coeffs)= coeffs%mu(:,1)
            coeff_vec(3*n_coeffs+1:4*n_coeffs)= coeffs%mu(:,4)
        else
            allocate(coeff_vec(2*n_coeffs*nz))
            do zc=1, nz
                coeff_vec(n_coeffs*(zc-1)+1:n_coeffs*zc) = coeffs%k(:,zc)
                coeff_vec(n_coeffs*(zc-1)+nz*n_coeffs+1:n_coeffs*zc+nz*n_coeffs)= coeffs%mu(:,zc)
            enddo
        endif

    end function makevector_coeffs
!-------------------------------------------------------------------------------

    pure function maketype_coeffs(coeff_vec, normalize) result(coeffs)
        type(tCoeffs) :: coeffs
        real(dp), dimension(:) ,intent(in) :: coeff_vec
        logical                ,intent(in) :: normalize
        integer :: zc

        allocate(coeffs%k(n_coeffs,nz), coeffs%mu(n_coeffs,nz), coeffs%r_squared(2,nz))
        coeffs%r_squared = 0.0

        if (pooled_regression) then
            do zc=1, nz
                coeffs%k(:,zc) = coeff_vec(1:n_coeffs)
                coeffs%mu(:,zc) = coeff_vec(n_coeffs+1:)
            enddo
        elseif (pi1_delta==1.0) then ! Alternatively, could have nz_1=count(stat_dist_z==0.0)
            coeffs%k(:,1) = coeff_vec(1:n_coeffs)
            coeffs%k(:,2) = 0.0
            coeffs%k(2,2) = 1.0
            coeffs%k(:,3) = coeffs%k(:,2)
            coeffs%k(:,4) = coeff_vec(n_coeffs+1:2*n_coeffs)
            coeffs%mu(:,1)= coeff_vec(2*n_coeffs+1:3*n_coeffs)
            coeffs%mu(:,2) = 0.0
            coeffs%mu(:,3) = coeffs%mu(:,2)
            coeffs%mu(:,4)= coeff_vec(3*n_coeffs+1:4*n_coeffs)
        else
            do zc=1, nz
                coeffs%k(:,zc)  = coeff_vec(n_coeffs*(zc-1)+1:n_coeffs*zc)
                coeffs%mu(:,zc) = coeff_vec(n_coeffs*(zc-1)+nz*n_coeffs+1:n_coeffs*zc+nz*n_coeffs)
            enddo
        endif

        if (normalize) then
	        coeffs%mu(1,:) = coeffs%mu(1,:)/norm_factor_small
	        coeffs%mu(2,:) = coeffs%mu(2,:)/norm_factor_large
	        if (n_coeffs >= 3) then
	            coeffs%k(3,:) = coeffs%k(3,:)/norm_factor_small
	            coeffs%mu(3,:) = coeffs%mu(3,:)/norm_factor_large
	        endif
        endif

    end function maketype_coeffs
!-------------------------------------------------------------------------------

    function initialize_coeffs(dir,n_coeffs,nz,estimate_from_simvars_o,agg_grid_o) result(coeffs)
        use classes_mod   ,only: tSimvars, tAggGrids
        use simvars_class ,only: read_unformatted

        type(tCoeffs)                          :: coeffs
        character(len=*) ,intent(in)           :: dir
        integer          ,intent(in)           :: n_coeffs, nz
        logical          ,intent(in) ,optional :: estimate_from_simvars_o
        type(tAggGrids)  ,intent(in) ,optional :: agg_grid_o
        real(dp) ,parameter :: coef_k2_init  = 0.95_dp ! could delete this
        real(dp) ,parameter :: coef_mu3_init = 0.70_dp
        integer             :: io_err
        type(tSimvars) ,allocatable :: simvars_old(:)

        allocate(coeffs%k(n_coeffs,nz), coeffs%mu(n_coeffs,nz), coeffs%r_squared(2,nz))
        coeffs%r_squared = 0.0
        coeffs%k         = 0.0
        coeffs%mu        = 0.0
ifdir:  if (dir == 'msge' .or. dir == 'mspe') then
            coeffs%k(2,:)  = 1.0     ! same for loms_in_logs

        elseif (dir == 'ge') then
            if (present(estimate_from_simvars_o)) then
		        if (estimate_from_simvars_o) then
		            print*, 'laws_of_motion:initialize_coeffs: estimating from old simvars'
		            call read_unformatted(simvars_old, io_err)
		            if (io_err == 0) then
		                call Regression(simvars_old,coeffs)
		                return
		            else
		                print*, 'WARNING: laws_of_motion:initialize_coeffs: could not read simvars'
		            endif
		        endif
            endif

            print*, 'laws_of_motion:initialize_coeffs: using hard-coded guess'
coef:       select case (n_coeffs)
            case(2)
pi:             if (pi1_delta==1.0) then

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
pooled:             if (pooled_regression) then
ir:                     if ((1.0 + scale_IR) > 0.1_dp) then
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
pi2:                    if (pi1_delta==0.5_dp) then
ir2:                       if ((1.0 + scale_IR) > 0.1_dp) then
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
                    coeffs%k(1,:)  = log(agg_grid_o%k(1))*(1-coeffs%k(2,:))
                    coeffs%mu(1,:) = log(agg_grid_o%mu(1))*(1-coeffs%mu(2,:))
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
                    if (loms_in_logs) then ! NEED TO IMPLEMENT
                        coeffs%k(1,:)  = 1.770695 !log(agg_grid_o%k(1))*(1-coeffs%k(2,:))
                        coeffs%mu(1,:) = -7.241165E+01 !log(agg_grid_o%mu(1))*(1-coeffs%mu(3,:))
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
                    else    ! NEED TO IMPLEMENT
                        coeffs%k(1,:)  = 1.770695 !log(agg_grid_o%k(1))*(1-coeffs%k(2,:))
                        coeffs%mu(1,:) = -7.241165E+01 !log(agg_grid_o%mu(1))*(1-coeffs%mu(3,:))
                    endif
                endif
            end select coef
        endif ifdir
     end function initialize_coeffs
!-------------------------------------------------------------------------------

pure subroutine Regression(simvars,coeffs)
    use params_mod   ,only: nt, t_scrap
    use classes_mod  ,only: tSimvars
    use MKL95_LAPACK ,only: gels, gelsy ! Intel MKL/ LAPACK dgels. Explanation at the bottom.

    type(tSimvars)           ,intent(in)   :: simvars(:)
    type(tCoeffs)            ,intent(inout):: coeffs   ! inout coz some coeffs might not be updated
    real(dp) ,dimension(:,:) ,allocatable  :: cov_k, cov_mu  ! covariates, rhs
    real(dp) ,dimension(:) ,allocatable    :: kp, mup              ! lhs
    real(dp) ,dimension(2) :: residual_sum_of_squares, total_sum_of_squares
    real(dp)               :: mean
    integer                :: n_covars, zc, n_zc(size(simvars)), n_zpc(size(simvars)), info, i
    logical, parameter :: use_gelsy = .false. ! instead of using gels (which had an error in Intel Composer XE 2013 original release)

    n_covars = size(coeffs%k,1)

    do zc=1,nz
        if (pooled_regression) then
            n_zc = nt-t_scrap
            n_zpc = n_zc

            allocate(cov_k (sum(n_zc), n_covars), kp (sum(n_zc) ))
            allocate(cov_mu(sum(n_zpc),n_covars), mup(sum(n_zpc)))
            cov_k(:,1)  = 1.0
            cov_k(:,2)  = [(simvars(i)%K (t_scrap:nt-1), i=1,size(simvars))]
            if (loms_in_logs) cov_k(:,2) = log(cov_k(:,2))
            if (n_covars > 2) cov_k(:,3) = cov_k(:,2)**2
            kp          = [(simvars(i)%K (t_scrap+1:nt), i=1,size(simvars))]

            cov_mu(:,1) = 1.0
            cov_mu(:,2) = [(simvars(i)%K (t_scrap+1:nt), i=1,size(simvars))]
            if (loms_in_logs) cov_mu(:,2) = log(cov_mu(:,2))
            if (n_covars > 2) cov_mu(:,3) = cov_mu(:,2)**2 !simvars%mu(t_scrap:nt-1)
            mup         = [(simvars(i)%mu(t_scrap+1:nt), i=1,size(simvars))]
        else
            do i=1,size(simvars)
                n_zc(i)   = count(simvars(i)%z(t_scrap:nt-1)==zc)
                n_zpc(i)  = count(simvars(i)%z(t_scrap+1:nt)==zc)
            enddo
            if (sum(n_zc) < 2) cycle     ! coeffs stay the same

            allocate(cov_k (sum(n_zc), n_covars), kp (sum(n_zc) ))
            allocate(cov_mu(sum(n_zpc),n_covars), mup(sum(n_zpc)))

            cov_k(:,1)  = 1.0
            cov_k(:,2)  = [(pack(simvars(i)%K (t_scrap:nt-1), simvars(i)%z(t_scrap:nt-1)==zc), i=1,size(simvars))]
            if (loms_in_logs) cov_k (:,2) = log(cov_k (:,2))
            if (n_covars > 2) cov_k(:,3)  = cov_k(:,2)**2
            kp          = [(pack(simvars(i)%K (t_scrap+1:nt), simvars(i)%z(t_scrap:nt-1)==zc), i=1,size(simvars))]

            cov_mu(:,1) = 1.0
            cov_mu(:,2) = [(pack(simvars(i)%K (t_scrap+1:nt), simvars(i)%z(t_scrap+1:nt)==zc), i=1,size(simvars))]
            if (loms_in_logs) cov_mu(:,2) = log(cov_mu(:,2))
            if (n_covars > 2) cov_mu(:,3) = cov_mu(:,2)**2 !pack(simvars%mu(t_scrap:nt-1), simvars%z(t_scrap+1:nt)==zc)
            mup         = [(pack(simvars(i)%mu(t_scrap+1:nt), simvars(i)%z(t_scrap+1:nt)==zc), i=1,size(simvars))]
        endif

        if (loms_in_logs) then
            kp            = log(kp)
!            mup           = log(mup)
        endif

        mean = sum(kp) /real(sum(n_zc),dp)
        total_sum_of_squares(1) = dot_product((kp -mean), (kp -mean))
        mean = sum(mup)/real(sum(n_zpc),dp)
        total_sum_of_squares(2) = dot_product((mup-mean), (mup-mean))

        if (n_covars == 4) then
            cov_k (:,4) = cov_k (:,2) *cov_k (:,3)
            cov_mu(:,4) = cov_mu(:,2) *cov_mu(:,3)
        elseif (n_covars == 5) then
            cov_k(:,4)  = cov_k (:,2)**2
            cov_k(:,5)  = cov_k (:,3)**2
            cov_mu(:,4) = cov_mu(:,2)**2
            cov_mu(:,5) = cov_mu(:,3)**2
        elseif (n_covars == 6) then
            cov_k(:,4)  = cov_k (:,2)**2
            cov_k(:,5)  = cov_k (:,3)**2
            cov_k(:,6)  = cov_k (:,2)*cov_k (:,3)
            cov_mu(:,4) = cov_mu(:,2)**2
            cov_mu(:,5) = cov_mu(:,3)**2
            cov_mu(:,6) = cov_mu(:,2)*cov_mu(:,3)
        endif

        if (use_gelsy) then ! This helps to check results. Intel Composer XE 13.0 had a bug in gels
            call gelsy(cov_k,kp, info=info) ! should print or return info
        else
            call gels(cov_k,kp, info=info) ! should print or return info
        endif
        coeffs%k(:,zc)  = kp(1:n_covars)
        if (use_gelsy) then
            call gelsy(cov_mu,mup, info=info)
        else
            call gels(cov_mu,mup, info=info)
        endif
        coeffs%mu(:,zc) = mup(1:n_covars)

        residual_sum_of_squares(1) = dot_product(kp (n_covars+1:),kp (n_covars+1:))
        residual_sum_of_squares(2) = dot_product(mup(n_covars+1:),mup(n_covars+1:))
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
        deallocate(cov_mu, cov_k,kp, mup)
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

pure real(dp) function Forecast(coeffs, k_in, mu_o)
    real(dp) ,intent(in)           :: coeffs(:), k_in
    real(dp) ,intent(in), optional :: mu_o
    real(dp)                       :: mu, k

    if (present(mu_o)) then
        mu= mu_o
    else
        mu = 0.0    ! ATTENTION, HACKISH
    endif

	if (loms_in_logs) then
	    k = log(k_in)
    else
        k = k_in
    endif

    select case (size(coeffs))
    case(2)
        Forecast = coeffs(1) + coeffs(2)*k
    case(3)
!        if (present(mu_in_o)) then
!            Forecast = coeffs(1) + coeffs(2)*k + coeffs(3)*mu
!        else
            Forecast = coeffs(1) + coeffs(2)*k + coeffs(3)*k**2
!        endif
    case(4)
        Forecast = coeffs(1) + coeffs(2)*k + coeffs(3)*mu + coeffs(4)*k*mu
    case(5)
        Forecast = coeffs(1) + coeffs(2)*k    + coeffs(3)*mu    &
                             + coeffs(4)*k**2 + coeffs(5)*mu**2
    case(6)
        Forecast = coeffs(1) + coeffs(2)*k    + coeffs(3)*mu    &
                             + coeffs(4)*k**2 + coeffs(5)*mu**2 + coeffs(6)*k*mu
    end select

    if (loms_in_logs .and. .not. present(mu_o)) Forecast = exp(Forecast)    ! VERY HACKISH to avoid exp of mu
    if (present(mu_o) .and. all(coeffs == 0.0) ) Forecast = mu ! This is hackish, should do better. It handles the cases ms and STY

end function Forecast


end module laws_of_motion
