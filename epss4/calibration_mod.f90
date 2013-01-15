module calibration_mod
! calibration of 'endogenous parameters'
    use kinds

    implicit none
    private
    public calibrate

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

        use params_mod            ,only: n_end_params, tol_calib
        use numrec_utils          ,only: put_diag
        use sub_alg_qn
        !use sub_broyden

        character(len=*) ,intent(in) :: projectname, calib_name
        real(dp) ,allocatable :: xvals(:), fvals(:), Rmat(:,:), QTmat(:,:)    ! QR decomposition in s_alg_qn
        real(dp)              :: maxstp
        logical               :: intialize_jacobi, not_converged
        integer               :: it
        integer, parameter    :: max_iterations = 3

        xvals = get_params()
        allocate(fvals(n_end_params))

        print *
        print '(t2,a)','- calibration: starting root finder'

        ! Initialize root finder
        it=0
        maxstp=1.0     ! this is large and arbitrary
        intialize_jacobi=.true.
        allocate(Rmat(n_end_params,n_end_params), QTmat(n_end_params,n_end_params))
        Rmat  = 0.0
        QTmat = 0.0
        call put_diag(1.0/0.5_dp,Rmat)

        ! Start root finder
        !call s_broyden(solve_krusellsmith, xvals, fvals,not_converged, tolf_o=tol_coeffs, maxstp_o = 0.5_dp, maxlnsrch_o=5) !df_o=Rmat,get_fd_jac_o=.true.
        call s_alg_qn(calibration_step,fvals,xvals,n_end_params,QTmat,Rmat,intialize_jacobi, &
             reevalj=.true.,check=not_converged,rstit0=10,MaxLns=5,max_it=max_iterations,maxstp=maxstp,tol_f=tol_calib)

        if (not_converged) then
            print *,' CRITICAL ERROR: Calibration didnt converge, fvals =', fvals
        else
            print *,' Calibration converged, fvals =', fvals
        endif

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
            call set_params(param_vec)
            print *,''
            print '(a,i3.3)','Calibration iteration ', it
            call run_model(projectname, calib_name, welfare_temp, simvars)

            distance = 1.0

            ! if (save_all_iterations) call save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, calib_name, projectname)

        end function calibration_step

    end subroutine calibrate
!-------------------------------------------------------------------------------

    pure function get_params()
        use params_mod ,only: n_end_params, beta, theta
        real(dp) ,allocatable ,dimension(:) :: get_params

        allocate(get_params(n_end_params))
        get_params(1) = beta
        if (n_end_params > 1) get_params(2) = theta
        ! if (n_end_params > 2) get_params(3) =

    end function get_params
!-------------------------------------------------------------------------------

    subroutine set_params(param_vec)
        use params_mod ,only: params_set
        real(dp) ,dimension(:) ,intent(in):: param_vec
        integer :: n

        n=size(param_vec)

        call params_set('beta',param_vec(1))
        if (n > 1) call params_set('theta',param_vec(2))

    end subroutine set_params
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
            prefix = ''
        endif

        open(unit=301, file='model_output/'prefix//cal_id(calib_name)//'_calibration.txt', status = 'replace')
        write(301,*) 'not_converged = ', not_converged
        write(301,*) 'fvals = ', fvals
        close(301)

    end subroutine write2file
!-------------------------------------------------------------------------------

end module calibration_mod
