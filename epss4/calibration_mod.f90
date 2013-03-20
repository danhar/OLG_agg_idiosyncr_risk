module calibration_mod
! calibration of 'endogenous parameters'
    use kinds

    implicit none
    private
    public calibrate, read_data_targets, model_targets, distance_norm

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - subroutine calibrate(projectname, calib_name)
! -- (internal) function krusellsmith(coeffvec) result(distance)
! - pure function get_params()
! - subroutine set_params(param_vec)
! - subroutine write2file(not_converged,fvals, calib_name)
!-------------------------------------------------------------------------------
    subroutine calibrate(projectname, calib_name)
    ! Sets up environment to use rootfinder,
    ! then pass the function krusell_smith as a function argument to a root finder.
    ! In this version, the function argument is an internal procedure, which is a thread-safe Fortran 2008 feature implemented
    ! in the Intel Fortran Compiler >= 11.0 and in gfortran >= 4.5
    ! This subroutine doesn't have a return value, since it sets the parameters globally.

        use params_mod   ,only: n_end_params, calib_targets, tol_calib
        use numrec_utils ,only: put_diag
        use sub_alg_qn   ,only: s_alg_qn
        use sub_zbrac    ,only: s_zbrac_array
        use fun_zbrent   ,only: zbrent
        !use sub_broyden

        character(len=*) ,intent(in) :: projectname, calib_name
        real(dp) ,dimension(:)   ,allocatable :: xvals, fvals, data_targets, norm_vector
        real(dp) ,dimension(:,:) ,allocatable :: Rmat, QTmat    ! QR decomposition in s_alg_qn
        real(dp)            :: maxstp, brack1, brack2
        logical             :: intialize_jacobi, not_converged, bracket_found
        integer             :: it
        integer  ,parameter :: max_iterations = 500
        real(dp) ,parameter :: brac_cover=0.005_dp
        logical  ,parameter :: use_brent_1D   = .false. ! Brent doesn't work very well b/c of KS guesses and aggregate grids.
        logical  ,parameter :: norm_params_to_1 = .true.  ! for Broyden

        xvals = get_params(n_end_params, trim(adjustl(calib_targets)))
        allocate(fvals(n_end_params))

        print *
        print '(t2,a)','- calibration: starting root finder'

        call read_data_targets(n_end_params, data_targets, trim(adjustl(calib_targets)))
        it=0

alg:    if (n_end_params == 1 .and. use_brent_1D) then ! Use a bracketing algorithm, i.e. Brent's Method (s_zbrent)

            print '(t2,a)','- calibration: WARNING: using Brent unreliable with K/S'
            ! First try to bracket a root by extending the bounds
            brack1 = xvals(1)*(1.0-brac_cover)
            brack2 = xvals(1)*(1.0+brac_cover)
            call s_zbrac_array(calibration_step, brack1, brack2, bracket_found)

            if (bracket_found) then
                print '(t2,a)','- calibration: bracketed a root, starting Brent'
                ! Now Brent's bracketing algorithm
                call zbrent(calibration_step, brack1, brack2,tol_calib, bracket_found, fvals(1))
                xvals(1)= brack2
            else
                print '(t2,a)','- calibration: ERROR, could not bracket root'
            endif
            not_converged = .not. bracket_found

        else alg ! Use a Broyden algorithm (s_alg_qn)

            ! Initialize root finder
            if (norm_params_to_1) then
                norm_vector = xvals
                xvals = xvals / norm_vector
            endif

            maxstp=.02_dp     ! this is crucial. 0.2_dp for IES=0.5, .02_dp for IES=1.5
            intialize_jacobi=.true.
            allocate(Rmat(n_end_params,n_end_params), QTmat(n_end_params,n_end_params))
            Rmat  = 0.0
            QTmat = 0.0
            call put_diag(1.0/0.5_dp,Rmat)

            ! Start root finder
            !call s_broyden(solve_krusellsmith, xvals, fvals,not_converged, tolf_o=tol_coeffs, maxstp_o = 0.5_dp, maxlnsrch_o=5) !df_o=Rmat,get_fd_jac_o=.true.
            call s_alg_qn(calibration_step,fvals,xvals,n_end_params,QTmat,Rmat,intialize_jacobi, &
                 reevalj=.true.,check=not_converged,rstit0=10,MaxLns=5,max_it=max_iterations,maxstp=maxstp,tol_f=tol_calib)
        endif alg

        print *, ''
        if (not_converged) then
            print *,' CRITICAL ERROR: Calibration didnt converge ' !, fvals =', fvals
        else
            print *,' *** Calibration converged *** ' !  fvals =', fvals
        endif
        print *, ''

        call write2file(not_converged,fvals, calib_name)

    contains

        function calibration_step(param_vec) result(distance)
            ! The whole model is solved and the distance is calculated from the moments of equilibrium simulations.
            ! This function has side effects, because it writes into the global parameters.

            use run_model_mod
            use classes_mod ,only: tSimvars

            real(dp), dimension(:), intent(in) :: param_vec
            real(dp), dimension(size(param_vec)):: distance
            type(tSimvars), allocatable :: simvars(:)
            real(dp) :: welfare_temp

            it = it+1
            if (norm_params_to_1) then
                call set_params(param_vec*norm_vector, trim(adjustl(calib_targets)))
            else
                call set_params(param_vec,  trim(adjustl(calib_targets)))
            endif
            print *,''
            print '(a,i3.3)','Calibration iteration ', it

            call run_model(projectname, calib_name, welfare_temp, simvars)

            distance = data_targets - model_targets(n_end_params, simvars, trim(adjustl(calib_targets)))

            ! if (save_all_iterations) call save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, calib_name, projectname)

        end function calibration_step

    end subroutine calibrate
!-------------------------------------------------------------------------------

    pure function get_params(n, targetname)
        use params_mod ,only: beta, theta, del_mean, del_std, pi1_delta, zeta_std
        real(dp) ,allocatable ,dimension(:) :: get_params
        integer          ,intent(in) :: n
        character(len=*) ,intent(in) :: targetname

        allocate(get_params(n))

        get_params(1) = beta
        select case (targetname)
        case('nosharpe','no_ep')
            if (n > 1) get_params(2) = del_std
            if (n > 2) get_params(3) = del_mean
            if (n > 3) get_params(4) = pi1_delta
            if (n > 4) get_params(5) = zeta_std
        case default
            if (n > 1) get_params(2) = del_std
            if (n > 2) get_params(3) = theta
            if (n > 3) get_params(4) = del_mean
            if (n > 4) get_params(5) = pi1_delta
            if (n > 5) get_params(6) = zeta_std
        end select

    end function get_params
!-------------------------------------------------------------------------------

    subroutine set_params(param_vec, targetname)
        use params_mod ,only: params_set, calibration_set_derived_params
        real(dp) ,dimension(:) ,intent(in) :: param_vec
        character(len=*)       ,intent(in) :: targetname
        integer :: n

        n=size(param_vec)

        call params_set('beta',param_vec(1))

        select case (targetname)
        case('nosharpe','no_ep')
            if (n > 1) call params_set('del_std',param_vec(2))
            if (n > 2) call params_set('del_mean',param_vec(3))
            if (n > 3) call params_set('pi1_delta',param_vec(4))
            if (n > 4) call params_set('zeta_std',param_vec(5))
        case default
            if (n > 1) call params_set('del_std',param_vec(2))
            if (n > 2) call params_set('theta',param_vec(3))
            if (n > 3) call params_set('del_mean',param_vec(4))
            if (n > 4) call params_set('pi1_delta',param_vec(5))
            if (n > 5) call params_set('zeta_std',param_vec(6))
        end select

        ! The following two calls set 'derived' parameters, e.g. gamma, which is a function of theta
        call calibration_set_derived_params()

    end subroutine set_params
!-------------------------------------------------------------------------------

    pure function model_targets(n, simvars, targetname)
        use classes_mod ,only: tSimvars, tStats
        use statistics  ,only: corr

        real(dp) ,allocatable:: model_targets(:)
        integer        ,intent(in) :: n
        type(tSimvars) ,intent(in) :: simvars(:)
        character(len=*)      ,intent(in) :: targetname
        type(tStats) :: K_Y, ex_ret, r, rf, zeta, netwage, cons_grow

        allocate(model_targets(n))

        K_Y%name='K_Y'; call K_Y%calc_stats(simvars)
        ex_ret%name='ex_ret'; call ex_ret%calc_stats(simvars)
        r%name='r'; call r%calc_stats(simvars)
        rf%name='rf'; call rf%calc_stats(simvars)
        zeta%name='zeta'; call zeta%calc_stats(simvars)
        netwage%name='netwage'; call netwage%calc_stats(simvars)
        cons_grow%name='cons_grow'; call cons_grow%calc_stats(simvars)

        select case (targetname)
        case('sharpe')
            model_targets(1) = rf%avg_exerr_()
            if (n > 1) model_targets(2) = cons_grow%std_()
            if (n > 2) model_targets(3) = ex_ret%avg_exerr_()/ex_ret%std_()
            if (n > 3) model_targets(4) = K_Y%avg_exerr_()
            if (n > 4) model_targets(5) = corr(zeta,r)
        case('nosharpe')
            model_targets(1) = rf%avg_exerr_()
            if (n > 1) model_targets(2) = cons_grow%std_()
            if (n > 2) model_targets(3) = K_Y%avg_exerr_()
            if (n > 3) model_targets(4) = corr(zeta,r)
        case('no_ep')
            model_targets(1) = K_Y%avg_exerr_()
            if (n > 1) model_targets(2) = r%std_()
            if (n > 2) model_targets(3) = r%avg_exerr_()
            if (n > 3) model_targets(4) = corr(zeta,r)
        case default
            model_targets(1) = K_Y%avg_exerr_()
            if (n > 1) model_targets(2) = r%std_()
            if (n > 2) model_targets(3) = ex_ret%avg_exerr_()
            if (n > 3) model_targets(4) = r%avg_exerr_()
            if (n > 4) model_targets(5) = corr(zeta,r)
            if (n > 5) model_targets(6) = netwage%cv_()
        end select

    end function model_targets
!-------------------------------------------------------------------------------

    subroutine read_data_targets(n, data_targets, filename)

        real(dp) ,allocatable ,intent(out):: data_targets(:)
        integer               ,intent(in) :: n
        character(len=*)      ,intent(in) :: filename
        integer  :: io_stat, line
        character(len=80) :: val, param, description

        print '(t2,a)','- calibration: reading calibration targets from model_input/data/calibration_targets/'//filename//'.txt'
        allocate(data_targets(n))
        line = 1

        open(unit=301, file='model_input/data/calibration_targets/'//filename//'.txt', status='OLD', form='formatted',iostat=io_stat, action='read')
        if (io_stat==0) then
            do
                if (line > n) exit
                read (301,*,iostat=io_stat) val, param, description
                if (scan(val,'!')>0) cycle
                if (io_stat/=0) then
                    print '(a,i6)', 'calibration_mod:read_data_targets: An error occured reading line', line
                    exit
                endif
                read (val,*) data_targets(line)
                line = line + 1
            enddo
        end if
        close (unit=301)

        if (io_stat .ne. 0) then
            print '(a,i6)', 'I/O ERROR: IOSTAT=',io_stat
            stop 'STOP in in calibration_mod:read_data_targets'
        endif

    end subroutine read_data_targets
!-------------------------------------------------------------------------------

    pure function distance_norm(data_targets, model_targets)
        ! This function is not needed here!
        real(dp) :: distance_norm
        real(dp), dimension(:), intent(in) :: data_targets, model_targets

        distance_norm = sum(abs(data_targets - model_targets))
    end function distance_norm
!-------------------------------------------------------------------------------

    subroutine write2file(not_converged, fvals, calib_name)
        use params_mod ,only: cal_id

        logical                ,intent(in) :: not_converged
        real(dp), dimension(:) ,intent(in) :: fvals
        character(len=*)       ,intent(in) :: calib_name
        character(:) ,allocatable :: prefix

        if (not_converged) then
            prefix = 'ERROR_'
        else
            prefix = 'success_'
        endif

        open(unit=301, file='model_output/'//prefix//'calibration_'//cal_id(calib_name)//'.txt', status = 'replace')
        write(301,*) 'not_converged = ', not_converged
        write(301,*) 'fvals = ', fvals
        close(301)

    end subroutine write2file
!-------------------------------------------------------------------------------

end module calibration_mod
