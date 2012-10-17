module krusell_smith

    use kinds
    use types
    use aggregate_grids,only: tAggGrids
    use laws_of_motion ,only: tCoeffs, MakeType, MakeVector, Regression
    use params_mod     ,only: exogenous_xgrid, save_all_iterations, normalize_coeffs
    use error_class
    implicit none

    private
    public solve_krusellsmith, set_krusellsmith, get_krusellsmith

    type(tCoeffs)    :: coeffs
    type(tPolicies)  :: policies
    type(tAggGrids)  :: grids
    type(tLifecycle) :: lifecycles
    type(tSimvars)   :: simvars
    type(tErrors)    :: err
    integer          :: it
    real(dp)     ,allocatable :: Phi_ms(:,:,:), Phi(:,:,:), value(:,:,:,:,:,:)
    character(:) ,allocatable :: projectname, calib_name

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - function solve_krusellsmith(coeffvec) result(distance)
! - subroutine set_krusellsmith(p,gr,lc,sv,e)
! - pure subroutine get_krusellsmith(p, v, ph, lc,sv,e)
!-------------------------------------------------------------------------------

	function solve_krusellsmith(coeffvec) result(distance)
        use simulate_economy
        use policyfunctions
        use interpolate_xgrid

	    real(dp), dimension(:), intent(in) :: coeffvec
	    real(dp), dimension(size(coeffvec)):: distance
	    type(tPolicies) :: pol_newx    ! if (exogenous_xgrid) then this will hold interpolated policies
	    type(tCoeffs)   :: coeffs_old, coeff_dif
	    real(dp), allocatable :: val_newx(:,:,:,:,:,:)

	    coeffs = MakeType(coeffvec, normalize_coeffs)
	    it = it+1

	    print '(t2,a43,i3.3)','- krusell_smith: solving for policies,  it = ', it
	    call calc_policyfunctions(coeffs, grids, policies, value, err)
	    call err%print2stderr

	    print *,'- krusell_smith: simulating'
	    Phi = Phi_ms
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

        if (save_all_iterations) call save_intermediate_results()

    contains
	!-------------------------------------------------------------------------------
	! Internal procedures in order:
	! - subroutine save_intermediate_results()
	!-------------------------------------------------------------------------------
	    subroutine save_intermediate_results()
	        use save_results_mod
            use ifport    ,only:  system     ! Intel Fortran portability library
	        use params_mod ,only:  pooled_regression, n_coeffs, construct_path
	        integer            :: syserr, zc
	        real(dp)           :: secs
	        character(:), allocatable   :: outpath
	        character(len=5) :: dir

	        outpath= construct_path('ge',calib_name)

	        open(unit=132, file=outpath//'/loms_it.txt', status = 'old', position='append')
	        write(132,'(a10,i3.3,a24,es13.6,a29,es13.6)')   &
	        'iteration ', it, ' :  max(abs(distance) = ', maxval(abs(distance)), ',  0.5*(distance*distance) = ', 0.5_dp*dot_product(distance,distance)
	        write(132,333) 'in:   ', coeffs_old%k(:,1),'   ---   ', coeffs_old%mu(:,1)
333         format(a6,<n_coeffs>(es13.6,x),a9,<n_coeffs>(es13.6,x))
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

	end function solve_krusellsmith

!-------------------------------------------------------------------------------
    subroutine set_krusellsmith(p,gr,ph,lc,sv,projname,calname)
        intent(in)        :: p,gr,ph,lc,sv,projname,calname
	    type(tPolicies)   :: p
	    type(tAggGrids)   :: gr
	    type(tLifecycle)  :: lc
	    type(tSimvars)    :: sv
	    real(dp)          :: ph(:,:,:)
	    character(len=*)  :: projname, calname

        policies   = p
        grids      = gr
        Phi_ms     = ph
        lifecycles = lc
        simvars    = sv
        projectname= projname
        calib_name = calname
        it         = 0
    end subroutine set_krusellsmith

!-------------------------------------------------------------------------------
    pure subroutine get_krusellsmith(p, v, ph, lc, sv, r2, e, its)
        intent(out)       :: p, v, ph, lc, sv, r2, its
        intent(inout)     :: e ! this is only because I also have e%not_converged
        type(tPolicies)   :: p
        real(dp), allocatable :: v(:,:,:,:,:,:), ph(:,:,:), r2(:,:)
        type(tLifecycle)  :: lc
        type(tSimvars)    :: sv
        type(tErrors)     :: e
        integer           :: its

        p  = policies
        v  = value
        ph = Phi
        lc = lifecycles
        sv = simvars
        r2 = coeffs%r_squared
        ! The trouble with err is that I do not want to overwrite e%not_converged
        e%asset = err%asset
	    e%cons  = err%cons
	    e%kp    = err%kp
	    e%mup   = err%mup
	    e%rfp   = err%rfp
        its = it
    end subroutine get_krusellsmith

end module krusell_smith
