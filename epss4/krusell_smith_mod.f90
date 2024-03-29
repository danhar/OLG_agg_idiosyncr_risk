!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

module krusell_smith_mod
    use kinds
    use classes_mod ,only: tSimvars, tLifecycle,tErrors, tAggGrids, tPolicies, tCoeffs, t_timer

    implicit none
    private
    public solve_krusellsmith

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - subroutine solve_krusellsmith(grids, projectname, calib_name, output_path, it, coeffs, simvars, Phi, policies, value, lifecycles, err, calibrating)
! -- (internal) function krusellsmith(coeffvec) result(distance)
! - subroutine save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, secs, dir, calib_name)
!-------------------------------------------------------------------------------
    subroutine solve_krusellsmith(grids, projectname, calib_name, output_path, input_path, it, coeffs, simvars, Phi, xgrid_ms, policies, value, lifecycles, err, calibrating, calc_euler_err)
    ! Set up environment to use rootfinder for coefficients of loms (KS),
    ! then pass the function krusell_smith as a function argument to a root finder.
    ! In this version, the function argument is an internal procedure, which is a thread-safe Fortran 2008 feature implemented
    ! in the Intel Fortran Compiler >= 11.0 and in gfortran >= 4.5

        use params_mod            ,only: normalize_coeffs, tol_coeffs, partial_equilibrium, save_all_iterations, construct_path, maxstp_ks
        use numrec_utils          ,only: put_diag
        use sub_alg_qn
        use alg_gs_mod
        !use sub_broyden

        type(tAggGrids) ,intent(in) :: grids
        character(len=*),intent(in)    :: projectname, calib_name, output_path, input_path
        type(tCoeffs)   ,intent(inout) :: coeffs
        type(tSimvars)  ,intent(inout) :: simvars(:)
        real(dp),allocatable ,intent(inout) :: Phi(:,:,:)
        real(dp)        ,intent(in)    :: xgrid_ms(:,:,:) ! Could remove if Phi was derived type carrying its own grid.
        logical         ,intent(in)    :: calibrating, calc_euler_err
        type(tPolicies) ,intent(out)   :: policies
        type(tLifecycle),intent(out)   :: lifecycles
        type(tErrors)   ,intent(out)   :: err
        integer         ,intent(out)   :: it
        real(dp) ,allocatable ,intent(out) :: value(:,:,:,:,:,:)

        real(dp) ,allocatable :: xvals(:), fvals(:), Rmat(:,:), QTmat(:,:)    ! QR decomposition in s_alg_qn
        real(dp) ,allocatable :: xgrid_mean_new(:,:,:)
        real(dp)              :: maxstp
        logical               :: intialize_jacobi
        integer               :: n, i, max_iter, recomp_jacobian
        logical, parameter :: ks_newton_alg = .true.

        coeffs%normalize = normalize_coeffs ! Can change it here or in params_mod (hidden from calibration file)
        call coeffs%save_initial_values()

        n= size(coeffs%makevector())
        allocate(xvals(n), fvals(n))

        it = 0
        xvals = coeffs%makevector()
        xgrid_mean_new = xgrid_ms ! First grid is mean shock grid. Could remove if Phi was derived type carrying its own grid.

        ! The following block is useful for debugging Krusell Smith
        if (save_all_iterations) then   ! write the headers into the file for saving intermediate coeffs
            open(unit=132, file=construct_path(calib_name)//'/loms_it.txt', status = 'replace')
            write(132,*)
            write(132,'(t8,a)') 'coeffs_k(1) , coeffs_k(2) , coeffs_k(3) , ...     ---    coeffs_mu(1), coeffs_mu(2), coeffs_mu(3) , ...'
            write(132,*)
            close(132)
        endif

        if (partial_equilibrium) then
            fvals = krusellsmith(xvals)

        else
            print *
            print '(t2,a)','- krusell_smith: starting root finder'

            ! Initialize root finder
            if (coeffs%normalize) then ! instead of if, could put maxstp in calibration file
                call coeffs%save_initial_values()
                maxstp=maxstp_ks      ! This makes sense, because coeffs are normalized to lie between 0.1 and 1.0. However, maxstp > 1 possible.
            else
                maxstp=10.0     ! this is large and arbitrary. Not =maxstp_ks here because coeffs%normalize not controlled from calibration file.
            endif

            if (calibrating) then
                recomp_jacobian = 15
                max_iter = 89        ! recomp_jacobian*6 -1
            else
                recomp_jacobian = 15
                max_iter = 499
            endif
            if (ks_newton_alg) then
                intialize_jacobi=.true.
                allocate(Rmat(n,n), QTmat(n,n))
                Rmat  = 0.0
                QTmat = 0.0
                call put_diag(1.0/0.5_dp,Rmat)

                ! Start root finder over coefficients of laws of motion
                !call s_broyden(solve_krusellsmith, xvals, fvals,not_converged, tolf_o=tol_coeffs, maxstp_o = 0.5_dp, maxlnsrch_o=5) !df_o=Rmat,get_fd_jac_o=.true.
                call s_alg_qn(krusellsmith,fvals,xvals,n,QTmat,Rmat,intialize_jacobi, &
                     reevalj=.true.,check=err%not_converged,rstit0=recomp_jacobian,MaxLns=5,max_it=max_iter,maxstp=maxstp,tol_f=tol_coeffs)
            else
                maxstp=.1_dp ! Not =maxstp_ks here because coeffs%normalize not controlled from calibration file.
                recomp_jacobian = 0
                max_iter = 400
                intialize_jacobi=.false.
                allocate(Rmat(n,n), QTmat(n,n))
                Rmat  = 0.0
                QTmat = 0.0
                call put_diag(1.0_dp,Rmat)

                ! Start root finder over coefficients of laws of motion
                call s_alg_qn(krusellsmith,fvals,xvals,n,QTmat,Rmat,intialize_jacobi, &
                     reevalj=.false.,check=err%not_converged,rstit0=recomp_jacobian,MaxLns=5,max_it=max_iter,maxstp=maxstp,tol_f=tol_coeffs)

!                call alg_gs(krusellsmith,fvals,xvals,n,err%not_converged,maxstp,max_iter,tol_coeffs)
            endif

            call coeffs%maketype(xvals) ! Need to get the xvals that led to the result, not the updated after the last regression.

            if (err%not_converged) call err%write2file(fvals, output_path)

        endif

    contains

        function krusellsmith(coeffvec) result(distance)
            ! Here we first solve for the policyfunctions, then simulate, then update the regression coefficients.
            ! This function has many side-effects. In particular, it writes directly into host's coeffs, so that in the end we have the latest update and also the R2 in that type.
            use lifecycles_class      ,only: average
            use simvars_class         ,only: print_error
            use params_mod            ,only: exogenous_xgrid, nx_factor
            use household_solution_mod,only: olg_backwards_recursion
            use laws_of_motion        ,only: Regression
            use simulation_mod        ,only: simulate
            use interpolate_xgrid

            real(dp), dimension(:), intent(in) :: coeffvec
            real(dp), dimension(size(coeffvec)):: distance
            type(tPolicies) :: pol_newx    ! if (exogenous_xgrid) then this will hold interpolated policies
            type(tCoeffs)   :: coeffs_old, coeff_dif
            type(tLifecycle) :: lifecycles_array(size(simvars))
            real(dp), allocatable :: val_newx(:,:,:,:,:,:), xgrid_mean_old(:,:,:), Phi_spread(:,:,:,:)
            type(t_timer):: timer
            logical ,parameter:: debugging_this = .false., timing = .false.

            call coeffs%maketype(coeffvec)
            it = it+1

            print '(t2,a43,i3.3)','- krusell_smith: solving for policies,  it = ', it
            if (timing) call timer%start()
            call olg_backwards_recursion(policies,coeffs, grids, value, err, input_path, 'ge')
            call err%print2stderr
            if (timing) then
                call timer%stop()
                print *, "Solution time = ", timer%elapsedTime()," sec"
                call timer%start()
            endif

            print *,'- krusell_smith: simulating'
            if (exogenous_xgrid) then
                ! This is the standard case which should always be used, because we make the xgrid much finer
                call InterpolateXgrid(nx_factor, policies, value, pol_newx, val_newx)
                ! We also want the initial Phi that we take from previous simulations to be defined over the new xgrid
                xgrid_mean_old = xgrid_mean_new ! Would be nicer to have a derived type Phi which carries its own xgrid.
                xgrid_mean_new = sum(pol_newx%xgrid(:,:,1,:,:,1),4)/size(pol_newx%xgrid,5) ! only an approximation of the grid over which Phi is defined
                call InterpolateXgrid(Phi, xgrid_mean_old, xgrid_mean_new)
                Phi_spread = spread(Phi,4,size(simvars))
                !$OMP  PARALLEL DO IF(size(simvars)>1)
                do i=1,size(simvars)
                    call simulate(pol_newx, val_newx, grids, coeffs, calc_euler_err, simvars(i), Phi_spread(:,:,:,i), lifecycles_array(i))
                enddo
                !$OMP END PARALLEL DO
            else
                ! In this case, we simultaneously solve for the rf(t+1) (actually the mu(t+1)) and the corresponding xgrid(since it depends on mu).
                ! That is very costly, so we do not refine xgrid. While more correct theoretically, the coarse xgrid makes the solution less precise.
                Phi_spread = spread(Phi,4,size(simvars))
                !$OMP  PARALLEL DO IF(size(simvars)>1)
                do i=1,size(simvars)
                    call simulate(policies, value, grids, coeffs, calc_euler_err, simvars(i), Phi_spread(:,:,:,i), lifecycles_array(i))
                enddo
                !$OMP END PARALLEL DO
            endif
            lifecycles=average(lifecycles_array)
            Phi = sum(Phi_spread,4)/size(Phi_spread,4)
            call print_error(simvars)

            if (timing) then
                call timer%stop()
                print *, "Simulation time = ", timer%elapsedTime()," sec"
                call timer%start()
            endif

            ! Here, one could make a (dampened) update of grids using the mean from simvars, but this turns out to be too complex for rootfinder

            coeffs_old    = coeffs
            call Regression(simvars,coeffs)

            coeff_dif%k  = (coeffs%k  - coeffs_old%k)  !/(coeffs_old%k+1.0) - this normalization is not necessary with root-finding. Also, it doesn't correspond to Judd (1998), p. 43: there, one takes the norm of the denominator.
            coeff_dif%mu = (coeffs%mu - coeffs_old%mu) !/(coeffs_old%mu+1.0)

            distance = coeff_dif%makevector()

            if (timing) then
                call timer%stop()
                print *, "Regression time = ", timer%elapsedTime()," sec"
            endif

            if (save_all_iterations) then
                print *,'- krusell_smith: WARNING: saving intermediate results (disk space)'
                call save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, calib_name, projectname)
            endif

        end function krusellsmith

    end subroutine solve_krusellsmith
!-------------------------------------------------------------------------------

    subroutine save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, calib_name, projectname)
        use save_results_mod
        use ifport     ,only: system     ! Intel Fortran portability library
        use params_mod ,only: pooled_regression, construct_path

        integer         ,intent(in) :: it
        real(dp)        ,intent(in) :: distance(:)
        type(tCoeffs)   ,intent(in) :: coeffs, coeffs_old
        type(tAggGrids) ,intent(in) :: grids
        character(len=*),intent(in) :: calib_name, projectname
        type(tSimvars)  ,intent(in) :: simvars(:)
        real(dp)        ,intent(in) :: Phi(:,:,:)
        type(tPolicies) ,intent(in) :: policies
        type(tLifecycle),intent(in) :: lifecycles
        type(tErrors)   ,intent(in) :: err
        integer            :: syserr, zc
        real(dp)           :: secs
        character(:), allocatable   :: outpath
        character(len=5) :: dir

        outpath= construct_path(calib_name)

        open(unit=132, file=outpath//'/loms_it.txt', status = 'old', position='append')
        write(132,'(a10,i3.3,a24,es13.6,a29,es13.6)')   &
        'iteration ', it, ' :  max(abs(distance) = ', maxval(abs(distance)), ',  0.5*(distance*distance) = ', 0.5_dp*dot_product(distance,distance)
        write(132,333) 'in:   ', coeffs_old%k(:,1),'   ---   ', coeffs_old%mu(:,1)
333     format(a6,<size(coeffs_old%k,1)>(es13.6,x),a9,<size(coeffs_old%mu,1)>(es13.6,x))
        if (.not. pooled_regression) then
            do zc=2,size(coeffs_old%k,2)
                write(132,333) '      ', coeffs_old%k(:,zc),'   ---   ', coeffs_old%mu(:,zc)
            enddo
            write(132,*) ''
        endif
        write(132,333) 'out:  ', coeffs%k(:,1),'   ---   ', coeffs%mu(:,1)
        if (.not. pooled_regression) then
            do zc=2,size(coeffs_old%k,2)
                write(132,333) '      ', coeffs%k(:,zc),'   ---   ', coeffs%mu(:,zc)
            enddo
        endif
        write(132,'(a99)') '---------------------------------------------------------------------------------------------------'
        close(132)

        write(dir, '(a2,i3.3)') 'it',it
        outpath= construct_path(calib_name,dir)
        syserr = system('mkdir '//outpath)
        secs= 0.0        ! No need to calc here
        call save_results(Phi, simvars, coeffs, grids,lifecycles,&
                         policies, secs, it, projectname, calib_name, dir, err)

    end subroutine save_intermediate_results

end module krusell_smith_mod
