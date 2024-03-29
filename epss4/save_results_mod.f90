!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

module save_results_mod
    implicit none
contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - subroutine save_results(Phi, simvars, coeffs, grids, lc, pol, secs, it, projectname, calib_name, dir, err)
! - subroutine plot_results(dir,mfile)
!-------------------------------------------------------------------------------
subroutine save_results(Phi, simvars, coeffs, grids, lc, &
                          pol, secs, it, projectname, calib_name, dir, err, cal_iter_o)

    use kinds
    use classes_mod     ,only: tPolicies, tAggGrids, tErrors, tSimvars, tLifecycle, tCoeffs, tStats, tStats_logical, tStats_integer
    use statistics      ,only: cov, corr
    use params_mod      ,only: n_eta, nj,nz, pop_frac, construct_path, cal_id, calc_euler_errors, check_dynamic_efficiency

    intent(in):: Phi, simvars, coeffs, grids, lc, pol, err, secs, it, projectname, calib_name, dir, cal_iter_o
    optional:: cal_iter_o

    type(tSimvars)   :: simvars(:) ! (kt, mut, bt,...), first element contains starting values
    type(tCoeffs)    :: coeffs
	type(tAggGrids)  :: grids   ! grids for aggregate states k and mu
    type(tLifecycle) :: lc      ! lifecycle profiles
    type(tPolicies)  :: pol
    type(tErrors)    :: err
	real(dp)         :: Phi(:,:,:), secs !distribution, seconds
	integer          :: it, i ! it = total iterations in Krusell Smith
    character(len=*) :: dir, projectname, calib_name, cal_iter_o ! current iteration in calibration routine

!    real(dp),dimension(nx,n_eta,nz,nj,size(pol%apgrid,5),size(pol%apgrid,6)):: cons
    real(dp), dimension(size(pol%apgrid,1),size(pol%apgrid,4)) :: apgrid_mean, stocks_mean, kappa_mean, xgrid_mean !, cons_mean
    type(tStats) :: K, mu, output, stock, bonds, invest, cons, cons_grow, netwage, transfers, pension, tau, net_mpk, r, rf, r_pf, r_pf_median, r_pf_kappa_med, zeta, delta, I_Y, K_Y, welfare, &
                    Phi_1, Phi_nx, err_aggr,B, err_inc, eul_err_max, eul_err_avg, bequest_rate, ex_ret, &
                    gini_income, gini_assets, gini_stocks, gini_consumption, cv_income, cv_assets, cv_stocks, cv_consumption ! inequality measures
    type(tStats_logical) :: err_K, err_mu, dyn_eff_a
    type(tStats_integer) :: dyn_eff_b
    character(:), allocatable :: path
    logical :: calibrating

    if (present(cal_iter_o)) then
        calibrating = .true.
        path = construct_path(calib_name)
    else
        calibrating = .false.
        path = construct_path(calib_name,dir)
    endif

!    cons = pol%xgrid -pol%apgrid

    call calc_stats()

    if (.not. calibrating) call write_stats('short')
    call write_stats('long')

    if (.not. calibrating) call write_output_files()

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - subroutine calc_stats()
! - subroutine write_stats(prec)
! - subroutine write_output_files()
!-------------------------------------------------------------------------------
    subroutine calc_stats()
	    integer  :: nk,nmu, xc, jc

        ! Assign character name and calculate the statistics
        K%name='K'; call K%calc_stats(simvars)
        mu%name='mu'; call mu%calc_stats(simvars)
        ex_ret%name='ex_ret'; call ex_ret%calc_stats(simvars)
        net_mpk%name='net_mpk'; call net_mpk%calc_stats(simvars)
        r%name='r'; call r%calc_stats(simvars)
        rf%name='rf'; call rf%calc_stats(simvars)
        r_pf%name='r_pf'; call r_pf%calc_stats(simvars)
        r_pf_median%name='rpf_med'; call r_pf_median%calc_stats(simvars)
        r_pf_kappa_med%name='rpf_kapm'; call r_pf_kappa_med%calc_stats(simvars)
        output%name='output'; call output%calc_stats(simvars)
        stock%name='stock'; call stock%calc_stats(simvars)
        bonds%name='bonds'; call bonds%calc_stats(simvars)
        invest%name='invest'; call invest%calc_stats(simvars)
        cons%name='cons'; call cons%calc_stats(simvars)
        netwage%name='netwage'; call netwage%calc_stats(simvars)
        transfers%name='transfers'; call transfers%calc_stats(simvars)
        pension%name='pension'; call pension%calc_stats(simvars)
        tau%name='tau'; call tau%calc_stats(simvars)
        K_Y%name='K_Y'; call K_Y%calc_stats(simvars)
        I_Y%name='I_Y'; call I_Y%calc_stats(simvars)
        welfare%name='welfare'; call welfare%calc_stats(simvars)
        Phi_1%name='Phi_1'; call Phi_1%calc_stats(simvars)
        Phi_nx%name='Phi_nx'; call Phi_nx%calc_stats(simvars)
        gini_income%name='gini_income'; call gini_income%calc_stats(simvars)
        gini_assets%name='gini_assets'; call gini_assets%calc_stats(simvars)
        gini_stocks%name='gini_stocks'; call gini_stocks%calc_stats(simvars)
        gini_consumption%name='gini_consumption'; call gini_consumption%calc_stats(simvars)
        cv_income%name='cv_income'; call cv_income%calc_stats(simvars)
        cv_assets%name='cv_assets'; call cv_assets%calc_stats(simvars)
        cv_stocks%name='cv_stocks'; call cv_stocks%calc_stats(simvars)
        cv_consumption%name='cv_consumption'; call cv_consumption%calc_stats(simvars)
        err_aggr%name='err_aggr'; call err_aggr%calc_stats(simvars)
        B%name='B'; call B%calc_stats(simvars)
        err_inc%name='err_inc'; call err_inc%calc_stats(simvars)
        eul_err_max%name = 'eul_err_max'; call eul_err_max%calc_stats(simvars)
        eul_err_avg%name = 'eul_err_avg'; call eul_err_avg%calc_stats(simvars)
        bequest_rate%name='bequests,%'; call bequest_rate%calc_stats(simvars)
        cons_grow%name='cons_grow'; call cons_grow%calc_stats(simvars)
        zeta%name='zeta'; call zeta%calc_stats(simvars)
        delta%name='delta'; call delta%calc_stats(simvars)
        err_K%name ='err_K' ; call err_K %calc_stats(simvars)
        err_mu%name='err_mu'; call err_mu%calc_stats(simvars)
        if (check_dynamic_efficiency) then
            dyn_eff_a%name ='dyn_eff_a' ; call dyn_eff_a%calc_stats(simvars)
            dyn_eff_b%name ='dyn_eff_b' ; call dyn_eff_b%calc_stats(simvars)
        endif

       ! The next holds for mean shock only if wm=0.25, which is ususally true
        nmu = size(pol%apgrid,6)
        nk  = size(pol%apgrid,5)
        do jc=1,size(pol%apgrid,4)
            do xc =1,size(pol%apgrid,1)
                apgrid_mean(xc,jc) = sum(pol%apgrid(xc,:,:,jc,:,:))/(nk*nmu*nz*n_eta)
                stocks_mean(xc,jc) = sum(pol%stocks(xc,:,:,jc,:,:))/(nk*nmu*nz*n_eta)
                xgrid_mean(xc,jc)  = sum(pol%xgrid(xc,:,:,jc,:,:)) /(nk*nmu*nz*n_eta)
!               cons_mean(xc,jc)   = sum(cons(xc,:,:,jc,:,:))  /(nk*nmu*nz*n_eta)
            enddo
        enddo
        where (apgrid_mean .ne. 0.0)
            kappa_mean = stocks_mean/apgrid_mean
        elsewhere
            kappa_mean = 0.0
        end where

    end subroutine calc_stats
!-------------------------------------------------------------------------------

    subroutine write_stats(prec)
    use params_mod, only: redistribute2workers
    character(len=*) ,intent(in) :: prec
    character(:), allocatable :: fmt1
    integer :: nl, show_digits
    integer ,parameter :: n_stats = 6 ! actually number of displayed statistics minus one

    nl = K%get_namelength()
    show_digits = K%get_digits2display(prec)
    if (calibrating) then
        open(unit=21, file=path//'/equilibrium_'//prec//cal_iter_o//'.txt', status = 'replace', action='write')
    else
        open(unit=21, file=path//'/equilibrium_'//prec//'.txt', status = 'replace', action='write')
    endif
    write(21,'(a12,a,",",a13,a)') ' Project    ', projectname, ' calibration ', calib_name
    write(21,*)
    write(21,*) 'Aggregate statistics'
    write(21,*) repeat('-',63)

    if (prec=='long') then
    write(21,123) repeat(' ',nl+2),' Average','  Median', 'Std.Dev.', 'Coef.Var', 'Autocorr', 'Abs.Max.', 'Abs.Min.'
    else
    write(21,123) repeat(' ',nl+2),' Average', 'Std.Dev.', 'Coef.Var', 'Autocorr'
    endif
123 format(a<nl+2>,a<show_digits+6>,<n_stats>(<show_digits+2>x,a8))
    fmt1 = K%writing_format(show_digits,nl) ! To guarantee compatibility when call K%write(21,'K',prec)

    call K%write(21,'K',prec)   ! in units of efficient labor, so that stationary
    call output%write(21,'output',prec)
    call invest%write(21,'investm',prec)
    call mu%write(21,'mu',prec)
    call ex_ret%write(21,'Ex_ret',prec)
    write(21,fmt1)' Sharpe    ', ex_ret%avg_exerr_()/ex_ret%std_()
    call net_mpk%write(21,'net_mpk',prec)
    call r%write(21,'r',prec)
    call rf%write(21,'rf',prec)
    if (prec == 'long') then
        call r_pf%write(21,'r_pf',prec)
        call r_pf_median%write(21,'rpf_med',prec)
        call r_pf_kappa_med%write(21,'rpf_kapm',prec)
    endif
    call cons%write(21,'cons',prec)   ! in units of efficient labor, so that stationary
    call cons_grow%write(21,'cons_grow',prec)
    call netwage%write(21,'netwage',prec)
    if (redistribute2workers) call transfers%write(21,'transfers',prec)
    call pension%write(21,'pension',prec)
    call tau%write(21,'tau',prec)
    call zeta%write(21,'zeta',prec)
    call I_Y%write(21,'I_Y',prec)
    call K_Y%write(21,'K_Y',prec)
    call welfare%write(21,'welfare',prec)

    if (prec == 'long') then
        write(21,*)
        call stock%write(21,'stock',prec)
        call bonds%write(21,'bonds',prec)
    endif

    write(21,*)
    write(21,*) 'Correlations'
    write(21,*) repeat('-',63)
    write(21,123)'            ',       '       r',             ' netwage',          '    zeta',            '  output',            '  invest',            '  ex_ret'
    write(21,fmt1)' r         ', 1.0
    write(21,fmt1)' netwage   ', corr(netwage  ,r), 1.0
    write(21,fmt1)' zeta      ', corr(zeta     ,r), corr(zeta     ,netwage), 1.0
    write(21,fmt1)' output    ', corr(output   ,r), corr(output   ,netwage), corr(output   ,zeta), 1.0
    write(21,fmt1)' invest    ', corr(invest   ,r), corr(invest   ,netwage), corr(invest   ,zeta), corr(invest   ,output), 1.0
    write(21,fmt1)' cons      ', corr(cons     ,r), corr(cons     ,netwage), corr(cons     ,zeta), corr(cons     ,output), corr(cons     ,invest), corr(cons     ,ex_ret)
    write(21,fmt1)' cons_grow ', corr(cons_grow,r), corr(cons_grow,netwage), corr(cons_grow,zeta), corr(cons_grow,output), corr(cons_grow,invest), corr(cons_grow,ex_ret)
    write(21,fmt1)' delta     ', corr(delta    ,r), corr(delta    ,netwage), corr(delta    ,zeta), corr(delta    ,output), corr(delta    ,invest), corr(delta    ,ex_ret)

    if (prec == 'long') then
        write(21,*)
        write(21,*) 'Covariances'
        write(21,*) repeat('-',63)
        write(21,123)'            ',       '       r',             ' netwage',          '    zeta',            '  output',            '  invest',            '  ex_ret'
        write(21,fmt1)' r         ', r%std_()**2
        write(21,fmt1)' netwage   ', cov(netwage  ,r), netwage%std_()**2
        write(21,fmt1)' zeta      ', cov(zeta     ,r), cov(zeta     ,netwage), zeta%std_()**2
        write(21,fmt1)' output    ', cov(output   ,r), cov(output   ,netwage), cov(output   ,zeta), output%std_()**2
        write(21,fmt1)' invest    ', cov(invest   ,r), cov(invest   ,netwage), cov(invest   ,zeta), cov(invest   ,output), invest%std_()**2
        write(21,fmt1)' cons      ', cov(cons     ,r), cov(cons     ,netwage), cov(cons     ,zeta), cov(cons     ,output), cov(cons     ,invest), cov(cons     ,ex_ret)
        write(21,fmt1)' cons_grow ', cov(cons_grow,r), cov(cons_grow,netwage), cov(cons_grow,zeta), cov(cons_grow,output), cov(cons_grow,invest), cov(cons_grow,ex_ret)
        write(21,fmt1)' delta     ', cov(delta    ,r), cov(delta    ,netwage), cov(delta    ,zeta), cov(delta    ,output), cov(delta    ,invest), cov(delta    ,ex_ret)

        write(21,*)
        write(21,*) 'Gini and Coefficient of Variation'
        write(21,*) repeat('-',63)
        write(21,123)'            ',         '  assets',         '  stocks',         '  income',          '    cons'
        write(21,fmt1)' Gini      ', gini_assets%avg_(), gini_stocks%avg_(), gini_income%avg_(), gini_consumption%avg_()
        write(21,fmt1)' Coef Var. ',   cv_assets%avg_(),   cv_stocks%avg_(),   cv_income%avg_(),   cv_consumption%avg_()
    endif

    write(21,*)
    write(21,*) 'Laws of motion'
    write(21,*) repeat('-',63)
    call coeffs%write(21, show_digits)

    write(21,*)
    write(21,*) 'Lifecycle statistics   (conditional on survival)' ! in per capita units
    write(21,*) repeat('-',63)
    write(21,123) repeat(' ',nl+2),' Average', 'Std.Dev.'
    write(21,fmt1) ' r_portf   ', sum(lc%return), sum(lc%return_var) !!! Attn: the second is Rpvar_avg!!!
    write(21,fmt1) ' cons      ', sum(lc%cons), sum(lc%cons_var)     !!! Attn: the second is cons_avg!!!
    if (prec=='long') then
        write(21,*) 'Report more lifecycle statistics in equil_long'
    endif

    write(21,*)
    write(21,*) 'Checks   (ExEr = exclude hits of aggregate grid bounds)'
    write(21,*) repeat('-',63)
    if (prec=='long') then
    write(21,123) repeat(' ',nl+2), 'Abs.Max.', 'Max_ExEr', ' Average', 'Avg_ExEr', 'Abs.Min.', 'Min_ExEr'
    else
    write(21,123) repeat(' ',nl+2), 'Abs.Max.', 'Max_ExEr', ' Average', 'Avg_ExEr'
    end if
    call B%write(21,'B',prec//'_max')
    call Phi_1%write(21,'Phi_1',prec//'_max')
    call Phi_nx%write(21,'Phi_nx',prec//'_max')
    if (prec=='long' .or. (err_inc%max_exerr_() > 1.0e-7)) call err_inc%write(21,'err_inc',prec//'_max')
    if (prec=='long' .or. (eul_err_max%max_exerr_() > 1.0e-1) .or. (eul_err_avg%max_exerr_() > 1.0e-2)) then
        call eul_err_max%write(21,'euler_max',prec//'_max')
        call eul_err_avg%write(21,'euler_avg',prec//'_max')
    endif
    call err_aggr%write(21,'err_aggr',prec//'_max')
    call bequest_rate%write(21,'bequests,%',prec//'_max')

    write(21,*)
125 format(3(a7,i6,a3,f4.1,a5))

    write(21,*)   'Warnings in solution'
    write(21,125) ' kp    ', count(err%kp) , '  (',real(count(err%kp ),dp)/real(size(err%kp ),dp)*100.0,'%)   ', &
                  '   mup ', count(err%mup), '  (',real(count(err%mup),dp)/real(size(err%mup),dp)*100.0,'%)   ', &
                  '   rfp ', count(err%rfp), '  (',real(count(err%rfp),dp)/real(size(err%rfp),dp)*100.0,'%)   '
    write(21,125) ' asset ', count(err%asset), '  (',real(count(err%asset),dp)/real(size(err%asset),dp)*100.0,'%)   ', &
                  '  cons ', count(err%cons) , '  (',real(count(err%cons ),dp)/real(size(err%cons ),dp)*100.0,'%)   '
    if (err%not_converged) write(21,*)'*** WARNING: root finder did not converge! ***'


    write(21,*)   'Warnings in simulation'
    write(21,125) ' K     ', err_K%count  , '  (',err_K%percent ,  '%)   ', &
                  '    mu ', err_mu%count , '  (',err_mu%percent,  '%)   '
    if (check_dynamic_efficiency) then
        write(21,*)   'Dynamic efficiency (Cond. A as violations, Cond. B as max)'
        write(21,125) ' CondA ', dyn_eff_a%count , '  (',dyn_eff_a%percent ,'%)   ' ,&
                      ' CondB ', dyn_eff_b%absmax, '  (',dyn_eff_b%avg     ,' avg)'
    endif

    write(21,*)   'Duration measures'
    write(21,'(a7,i6, 14x, a5, es9.2)') ' KS_it ', it, 'secs ', secs

    close(21)
    end subroutine write_stats
!-------------------------------------------------------------------------------

    subroutine write_output_files()
    ! Since all arrays follow Fortran's natural storage order (i.e. column-major),
    ! no (implied) do-loops necessary when writing to file.
    ! The following will write all nx in one line, then change zc, then jc, ...
    use params_mod, only: save_all_to_txt, redistribute2workers
    integer :: nx, i2

! The following seems the best format, because it will not produce ****, which can't be read by Matlab
301 format(<nx> (es13.6,1x))
302 format(<nj> (es13.6,1x))
303 format(<nx> (es14.6e3,1x))
!360 format(<nj> (f0.6,x))

    nx = size(pol%apgrid,1)
!    open(20, file=path//'/cons.txt', status = 'replace')
!    write(20,302) cons
!    close(20)

    if (save_all_to_txt) then
        open(40, file=path//'/apgrid.txt',  status = 'replace')
        write(40,301) pol%apgrid
        close(40)

        open(50, file=path//'/xgrid.txt',   status = 'replace')
        write(50,301) pol%xgrid
        close(50)

        open(60, file=path//'/kappa.txt',   status = 'replace')
        write(60,301) pol%kappa
        close(60)

        ! The next is important, but also large.
        call err%write2file(path)
    else
        print*, 'Note: not saving complete policies and all errors to txt-file.'
    endif

    open(40, file=path//'/apgrid_mean.txt',  status = 'replace')
    write(40,301) apgrid_mean
    close(40)

    open(50, file=path//'/xgrid_mean.txt',   status = 'replace')
    write(50,301) xgrid_mean
    close(50)

!    open(50, file=path//'/cons_mean.txt',   status = 'replace')
!    write(50,301) cons_mean
!    close(50)

    open(60, file=path//'/kappa_mean.txt',   status = 'replace')
    write(60,301) kappa_mean
    close(60)

    open(60, file=path//'/stocks_mean.txt',   status = 'replace')
    write(60,301) stocks_mean
    close(60)

    nx = size(Phi,1)
    open(unit=21, file=path//'/Phi_tilde.txt', status = 'replace')
    write(21,303) sum(Phi,2) ! This is a linear approximation, which corresponds to the xgrid below
    close(21)

    open(unit=21, file=path//'/Phi.txt', status = 'replace')
    write(21,303) Phi
    close(21)

    open(unit=21, file=path//'/exp_val_tilde.txt', status = 'replace')
    write(21,303) sum(lc%exp_value,2)
    close(21)

    open(unit=21, file=path//'/xgrid_tilde.txt', status = 'replace')
    write(21,301) sum(lc%xgrid,2)/size(lc%xgrid,2)
    close(21)

    open(unit=21, file=path//'/ap_lc.txt', status = 'replace')
    write(21,302) lc%ap
    close(21)

    open(unit=21, file=path//'/kappa_lc.txt', status = 'replace')
    write(21,302) lc%kappa
    close(21)

    open(unit=21, file=path//'/cons_lc.txt', status = 'replace')
    write(21,302) lc%cons
    close(21)

    open(unit=21, file=path//'/stock_lc.txt', status = 'replace')
    write(21,302) lc%stock
    close(21)

    open(unit=21, file=path//'/consvar_lc.txt', status = 'replace')
    write(21,302) lc%cons_var
    close(21)

    open(unit=21, file=path//'/return_lc.txt', status = 'replace')
    write(21,302) lc%return
    close(21)

    open(unit=21, file=path//'/return_var_lc.txt', status = 'replace')
    write(21,302) lc%return_var
    close(21)

    open(unit=21, file=path//'/log_cons_lc.txt', status = 'replace')
    write(21,302) lc%log_cons
    close(21)

    open(unit=21, file=path//'/var_log_cons_lc.txt', status = 'replace')
    write(21,302) lc%var_log_cons
    close(21)

    open(unit=21, file=path//'/pop_frac.txt', status = 'replace')
    write(21,302) pop_frac
    close(21)

    open(unit=21, file=path//'/agg_grids.txt', status = 'replace')
    write(21,*) ' grids%k   =  '
    write(21,'(<size(grids%k)>(es13.6,1x))') grids%k
    write(21,*) ' grids%mu  =  '
    write(21,'(<size(grids%mu)>(es13.6,1x))') grids%mu
    close(21)

    open(unit=21, file=path//'/simvars.txt', status = 'replace')
    do i=1,size(simvars)
        write(21,'(a,i3,a,i3)') 'Parallel simulation ',i,' of ',size(simvars)
        write(21,368) ' t         ', [(i2,i2=1,size(simvars(i)%z))]
    368 format(a11,<size(simvars(i)%z)> (i6.6,4x))
        write(21,369) ' z        ', simvars(i)%z
    369 format(a10,<size(simvars(i)%z)> (i2,8x))
        write(21,371) ' K        ', simvars(i)%K
    371 format(a10,<size(simvars(i)%K)> (f9.5,1x))
        write(21,370) ' mu       ', simvars(i)%mu
    370 format(a10,<size(simvars(i)%z)> (f9.5,1x))
        write(21,370) ' output   ', simvars(i)%output
        write(21,370) ' invest   ', simvars(i)%invest
        write(21,370) ' stock    ', simvars(i)%stock
        write(21,370) ' bonds    ', simvars(i)%bonds
        write(21,370) ' cons     ', simvars(i)%C
        write(21,370) ' net_mpk  ', simvars(i)%net_mpk
        write(21,370) ' r        ', simvars(i)%r
        write(21,371) ' rf       ', simvars(i)%rf
        write(21,370) ' netwage  ', simvars(i)%wage
        if (redistribute2workers) write(21,370) ' transfer ', simvars(i)%trans
        write(21,370) ' pension  ', simvars(i)%pens
        write(21,370) ' tau      ', simvars(i)%tau
        write(21,373) ' welf     ', simvars(i)%welf
        write(21,372) ' err_K    ', simvars(i)%err_K
        write(21,372) ' err_mu   ', simvars(i)%err_mu
    372 format(a10,<size(simvars(i)%z)> (1x,l1,8x))
        write(21,373) ' Phi_1    ', simvars(i)%Phi_1
        write(21,373) ' Phi_nx   ', simvars(i)%Phi_nx
        write(21,370) ' B        ', simvars(i)%B
        write(21,373) ' err_inc  ', simvars(i)%err_income
        write(21,373) ' err_aggr ', simvars(i)%err_aggr
        if (calc_euler_errors) then
            write(21,373) ' eul_max  ', simvars(i)%eul_err_max
            write(21,373) ' eul_avg  ', simvars(i)%eul_err_avg
        endif
        if (check_dynamic_efficiency) then
            write(21,372)  ' dyn_eff_a', simvars(i)%dyn_eff_a
            write(21,368) ' dyn_eff_b ', simvars(i)%dyn_eff_b
        endif
        write(21,373) ' bequests ', simvars(i)%bequests
        write(21,373) ' gini_inc ', simvars(i)%gini_income
        write(21,373) ' gi_asset ', simvars(i)%gini_assets
        write(21,373) ' gi_stock ', simvars(i)%gini_stocks
        write(21,373) ' gi_cons  ', simvars(i)%gini_consumption
        write(21,373) ' cv_inc   ', simvars(i)%cv_income
        write(21,373) ' cv_asset ', simvars(i)%cv_assets
        write(21,373) ' cv_stock ', simvars(i)%cv_stocks
        write(21,373) ' cv_cons  ', simvars(i)%cv_consumption
    373 format(a10,<size(simvars(i)%z)> (es9.2,1x))
        write(21,*)
    enddo
    close(21)

    end subroutine write_output_files

end subroutine save_results

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine plot_results(path,mfile)
    ! Call Matlab and plot
    use ifport            ,only: system  ! Intel Fortran portability library
    use params_mod        ,only: no_plotting
    character(len=*), intent(in)   :: path
    character(len=*), intent(in)   :: mfile
    integer                        :: err_matl
    character(:), allocatable      :: epssdir

    if (no_plotting) return
    !print*, ' '
    !print*, 'Starting MATLAB to plot ',mfile,'(',dir,')'
    err_matl = system('mkdir '//path//'/graphs  > /dev/null 2>&1')
    epssdir = path(1+scan(path,'/'):)
    err_matl = system('EPSSDIR='//epssdir//' && export EPSSDIR && cd src_matlab && matlab -nodesktop -nosplash -nodisplay -r '//mfile//' > cl_output.txt')
    if (err_matl ==-1) print*, 'Warning in save_results:plot: An error occured while plotting graphs'
end subroutine plot_results

end module save_results_mod
