module krusell_smith_mod
    use kinds
    use types                  ,only: tSimvars, tLifecycle
    use error_class            ,only: tErrors
    use laws_of_motion         ,only: tCoeffs
    use aggregate_grids        ,only: tAggGrids
    use household_solution_mod ,only: tPolicies

    implicit none
    private
    public solve_krusellsmith

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - subroutine solve_krusellsmith(grids, projectname, calib_name, output_path, it, coeffs, simvars, Phi, policies, value, lifecycles, err)
! -- (internal) function solve_krusellsmith(coeffvec) result(distance)
! - subroutine save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, secs, dir, calib_name)
!-------------------------------------------------------------------------------
    subroutine solve_krusellsmith(grids, projectname, calib_name, output_path, it, coeffs, simvars, Phi, policies, value, lifecycles, err)
    ! Set up environment to use rootfinder for coefficients of loms (KS)

        use household_solution_mod,only: olg_backwards_recursion
        use laws_of_motion        ,only: MakeType, MakeVector, Regression
        use params_mod            ,only: exogenous_xgrid, save_all_iterations, normalize_coeffs, tol_coeffs, partial_equilibrium
        use numrec_utils          ,only: put_diag
        use sub_alg_qn
        !use sub_broyden

        type(tAggGrids) ,intent(in)    :: grids
        character(len=*),intent(in)    :: projectname, calib_name, output_path
        type(tCoeffs)   ,intent(inout) :: coeffs
        type(tSimvars)  ,intent(inout) :: simvars
        real(dp)        ,intent(inout) :: Phi(:,:,:)
        type(tPolicies) ,intent(out)   :: policies
        type(tLifecycle),intent(out)   :: lifecycles
        type(tErrors)   ,intent(out)   :: err
        integer         ,intent(out)   :: it
        real(dp) ,allocatable ,intent(out) :: value(:,:,:,:,:,:)

        real(dp) ,allocatable :: xvals(:), fvals(:), Rmat(:,:), QTmat(:,:)    ! QR decomposition in s_alg_qn
        real(dp)              :: maxstp
        logical               :: intialize_jacobi
        integer               :: n

        if (normalize_coeffs) then ! instead of if, could put maxstp in calibration file
            maxstp=1.0_dp      ! This makes sense, because coeffs are normalized to lie between 0.1 and 1.0
        else
            maxstp=10.0     ! this is large and arbitrary
        endif

        intialize_jacobi=.true.
        n= size(MakeVector(coeffs, normalize_coeffs))
        allocate(xvals(n), fvals(n), Rmat(n,n), QTmat(n,n))
        Rmat  = 0.0
        QTmat = 0.0
        call put_diag(1.0/0.5_dp,Rmat)

        it = 0
        xvals = MakeVector(coeffs, normalize_coeffs)

        if (partial_equilibrium) then
            fvals = krusellsmith(xvals)
        else
            !call s_broyden(solve_krusellsmith, xvals, fvals,not_converged, tolf_o=tol_coeffs, maxstp_o = 0.5_dp, maxlnsrch_o=5) !df_o=Rmat,get_fd_jac_o=.true.
            call s_alg_qn(krusellsmith,fvals,xvals,n,QTmat,Rmat,intialize_jacobi, &
                 reevalj=.true.,check=err%not_converged,rstit0=10,MaxLns=5,max_it=100,maxstp=maxstp,tol_f=tol_coeffs) ! maxstp=1.0_dp

            if (err%not_converged) call err%write2file(fvals, output_path)
            coeffs = MakeType(xvals, normalize_coeffs)

        endif

    contains

        function krusellsmith(coeffvec) result(distance)
            use simulate_economy
            use interpolate_xgrid

            real(dp), dimension(:), intent(in) :: coeffvec
            real(dp), dimension(size(coeffvec)):: distance
            type(tPolicies) :: pol_newx    ! if (exogenous_xgrid) then this will hold interpolated policies
            type(tCoeffs)   :: coeffs_old, coeff_dif
            real(dp), allocatable :: val_newx(:,:,:,:,:,:)

            coeffs = MakeType(coeffvec, normalize_coeffs)
            it = it+1

            print '(t2,a43,i3.3)','- krusell_smith: solving for policies,  it = ', it
            call olg_backwards_recursion(policies,coeffs, grids, value, err)
            call err%print2stderr

            print *,'- krusell_smith: simulating'
            if (exogenous_xgrid) then
                call InterpolateXgrid(policies, value, pol_newx, val_newx) ! No need to interpolate Phi
                call simulate(pol_newx, val_newx, grids, simvars, Phi, lifecycles)
            else
                call simulate(policies, value, grids, simvars, Phi, lifecycles)
            endif
            call print_error_msg(simvars)
            ! One could make a (dampened) update of grids using the mean from simvars, but too complex for rootfinder

            coeffs_old    = coeffs
            call Regression(simvars,coeffs)

            coeff_dif%k  = (coeffs%k  - coeffs_old%k)  !/(coeffs_old%k+1.0)
            coeff_dif%mu = (coeffs%mu - coeffs_old%mu) !/(coeffs_old%mu+1.0)

            distance = MakeVector(coeff_dif, normalize_coeffs)

            if (save_all_iterations) call save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, calib_name, projectname)

        end function krusellsmith

    end subroutine solve_krusellsmith
!-------------------------------------------------------------------------------

    subroutine save_intermediate_results(it, distance, coeffs, coeffs_old, Phi, simvars, grids, lifecycles, policies, err, calib_name, projectname)
        use save_results_mod
        use ifport     ,only: system     ! Intel Fortran portability library
        use params_mod ,only: pooled_regression, n_coeffs, construct_path

        integer         ,intent(in) :: it
        real(dp)        ,intent(in) :: distance(:)
        type(tCoeffs)   ,intent(in) :: coeffs, coeffs_old
        type(tAggGrids) ,intent(in) :: grids
        character(len=*),intent(in) :: calib_name, projectname
        type(tSimvars)  ,intent(in) :: simvars
        real(dp)        ,intent(in) :: Phi(:,:,:)
        type(tPolicies) ,intent(in) :: policies
        type(tLifecycle),intent(in) :: lifecycles
        type(tErrors)   ,intent(in) :: err
        integer            :: syserr, zc
        real(dp)           :: secs
        character(:), allocatable   :: outpath
        character(len=5) :: dir

        outpath= construct_path('ge',calib_name)

        open(unit=132, file=outpath//'/loms_it.txt', status = 'old', position='append')
        write(132,'(a10,i3.3,a24,es13.6,a29,es13.6)')   &
        'iteration ', it, ' :  max(abs(distance) = ', maxval(abs(distance)), ',  0.5*(distance*distance) = ', 0.5_dp*dot_product(distance,distance)
        write(132,333) 'in:   ', coeffs_old%k(:,1),'   ---   ', coeffs_old%mu(:,1)
333     format(a6,<n_coeffs>(es13.6,x),a9,<n_coeffs>(es13.6,x))
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
        outpath= construct_path(dir,calib_name)
        syserr = system('mkdir '//outpath)
        secs= 0.0        ! No need to calc here
        call save_results(Phi, simvars, coeffs, grids,lifecycles,&
                         policies, secs, it, projectname, calib_name, dir, err)

    end subroutine save_intermediate_results

end module krusell_smith_mod
