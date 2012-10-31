module save_results_mod
    implicit none
contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - subroutine save_results(Phi, simvars, coeffs, grids, lc, pol, secs, it, projectname, calib_name, dir, err)
! - subroutine plot_results(dir,mfile)
!-------------------------------------------------------------------------------
subroutine save_results(Phi, simvars, coeffs, grids, lc, &
                          pol, secs, it, projectname, calib_name, dir, err)

    use kinds
    use types
    use policyfunctions
    use statistics      ,only: tStats, covariance
    use aggregate_grids ,only: tAggGrids
    use laws_of_motion  ,only: tCoeffs
    use error_class      ,only: tErrors
    use params_mod      ,only: n_eta, nj,nx,nz, n_coeffs, pop_frac, construct_path, loms_in_logs

    intent(in):: Phi, simvars, coeffs, grids, lc, pol, err, secs, it, projectname, calib_name, dir

    type(tSimvars)   :: simvars ! (kt, mut, bt,...), first element contains starting values
    type(tCoeffs)    :: coeffs
	type(tAggGrids)  :: grids   ! grids for aggregate states k and mu
    type(tLifecycle) :: lc      ! lifecycle profiles
    type(tPolicies)  :: pol
    type(tErrors)    :: err
	real(dp)         :: Phi(nx,n_eta,nj), secs !distribution, seconds
	integer          :: it
    character(len=*) :: dir, projectname, calib_name

!    real(dp),dimension(nx,n_eta,nz,nj,size(pol%apgrid,5),size(pol%apgrid,6)):: cons
    real(dp), dimension(nx,nj) :: apgrid_mean, stocks_mean, kappa_mean, xgrid_mean !, cons_mean
    type(tStats) :: K, mu, output, stock, bonds, invest, cons, cg, netwage, pens, tau, r, rf, r_pf_median, r_pf_kappa_med, zeta, delta, K_Y, welfare, &
                    Phi_1, Phi_nx, err_aggr,B, err_inc, bequest_rate, ex_ret
    real(dp) :: percent_err_K, percent_err_mu, cov_w_r, cov_zeta_r, cov_zeta_w, cov_output_r, cov_output_w, cov_output_zeta, cov_inv_r, cov_inv_w, cov_inv_zeta, cov_inv_output, &
                cov_cons_r, cov_cons_w, cov_cons_zeta, cov_cons_output, cov_cons_invest, cov_cons_ex_ret, cov_cg_r, cov_cg_w, cov_cg_zeta, cov_cg_output, cov_cg_invest, cov_cg_ex_ret, &
                cov_delta_r, cov_delta_w, cov_delta_zeta, cov_delta_output, cov_delta_invest, cov_delta_ex_ret

    integer  :: count_err_K, count_err_mu, i
    character(:), allocatable :: path

    path = construct_path(dir, calib_name)
!    cons = pol%xgrid -pol%apgrid

    call calc_stats()

    call write_stats('short')
    call write_stats('long')

    call write_output_files()

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - subroutine calc_stats()
! - subroutine write_stats()
! - subroutine write_output_files()
!-------------------------------------------------------------------------------
    subroutine calc_stats()
	    use params_mod, only: alpha, nt, t_scrap
	    use income, only: zetaval => zeta, deltaval => delta
	    real(dp), dimension(:), allocatable :: simvars_zeta, simvars_delta
	    real(dp), dimension(size(simvars%C)-1) :: C_growth
	    real(dp), dimension(size(simvars%r))   :: excess_return
	    integer  :: nk,nmu, xc, jc, lb

        associate(err_mu => simvars%err_mu, err_k => simvars%err_K)

        excess_return = simvars%r - simvars%rf

        call K%calc_stats(simvars%K,   err_mu, err_K)
        call mu%calc_stats(simvars%mu, err_mu, err_K)
        call ex_ret%calc_stats(excess_return, err_mu, err_K)
        call r%calc_stats(simvars%r, err_mu, err_K)
        call rf%calc_stats(simvars%rf, err_mu, err_K)
        call r_pf_median%calc_stats(simvars%r_pf_median, err_mu, err_K)
        call r_pf_kappa_med%calc_stats(simvars%r_pf_kappa_med, err_mu, err_K)
        call output%calc_stats(simvars%output,   err_mu, err_K)
        call stock%calc_stats(simvars%stock,   err_mu, err_K)
        call bonds%calc_stats(simvars%bonds,   err_mu, err_K)
        call invest%calc_stats(simvars%invest,   err_mu, err_K)
        call cons%calc_stats(simvars%C,   err_mu, err_K)
        call netwage%calc_stats(simvars%wage,   err_mu, err_K)
        call pens%calc_stats(simvars%pens,   err_mu, err_K)
        call tau%calc_stats(simvars%tau,   err_mu, err_K)
        call K_Y%calc_stats(simvars%K**(1.0-alpha),   err_mu, err_K)
        call welfare%calc_stats(simvars%welf,   err_mu, err_K)
        call Phi_1%calc_stats(simvars%Phi_1,   err_mu, err_K)
        call Phi_nx%calc_stats(simvars%Phi_nx,   err_mu, err_K)
        call err_aggr%calc_stats(abs(simvars%err_aggr),   err_mu, err_K)
        call B%calc_stats(simvars%B,   err_mu, err_K)
        call err_inc%calc_stats(abs(simvars%err_income),   err_mu, err_K)
        call bequest_rate%calc_stats(simvars%bequests/simvars%K**alpha,   err_mu, err_K)

        allocate(simvars_zeta(size(simvars%z)))
        if (dir=='msge' .or. dir=='mspe') then
            lb = 2
            C_growth    = 0.0
            call cg%set_number(0.0_dp)
            simvars_zeta  = [1.0_dp ,  zetaval(simvars%z(lb:))]
            simvars_delta = [sum(deltaval)/size(deltaval) , deltaval(simvars%z(lb:))]
        else
            lb = t_scrap +1
            C_growth    = (simvars%C(2:) - simvars%C(1:size(simvars%C)-1))/simvars%C(1:size(simvars%C)-1)
            call cg%calc_stats(C_growth, err_mu, err_K)
            simvars_zeta =  zetaval(simvars%z)
            simvars_delta= deltaval(simvars%z)
        endif
        call  zeta%calc_stats(simvars_zeta, err_mu, err_K)
        call delta%calc_stats(simvars_delta,err_mu, err_K)

        cov_w_r         = covariance(simvars%wage,simvars%r,err_mu, err_K)
        cov_zeta_r      = covariance(simvars_zeta,simvars%r,err_mu, err_K)
        cov_zeta_w      = covariance(simvars_zeta,simvars%wage,err_mu, err_K)
        cov_output_r    = covariance(simvars%output,simvars%r,err_mu, err_K)
        cov_output_w    = covariance(simvars%output,simvars%wage,err_mu, err_K)
        cov_output_zeta = covariance(simvars%output,simvars_zeta,err_mu, err_K)
        cov_inv_r       = covariance(simvars%invest,simvars%r,err_mu, err_K)
        cov_inv_w       = covariance(simvars%invest,simvars%wage,err_mu, err_K)
        cov_inv_zeta    = covariance(simvars%invest,simvars_zeta,err_mu, err_K)
        cov_inv_output  = covariance(simvars%invest,simvars%output,err_mu, err_K)
        cov_cons_r      = covariance(simvars%C,simvars%r,err_mu, err_K)
        cov_cons_w      = covariance(simvars%C,simvars%wage,err_mu, err_K)
        cov_cons_zeta   = covariance(simvars%C,simvars_zeta,err_mu, err_K)
        cov_cons_output = covariance(simvars%C,simvars%output,err_mu, err_K)
        cov_cons_invest = covariance(simvars%C,simvars%invest,err_mu, err_K)
        cov_cons_ex_ret = covariance(simvars%C,excess_return,err_mu, err_K)
        cov_cg_r        = covariance(C_growth,simvars%r(2:size(simvars%r)),err_mu, err_K)
        cov_cg_w        = covariance(C_growth,simvars%wage(2:size(simvars%wage)),err_mu, err_K)
        cov_cg_zeta     = covariance(C_growth,simvars_zeta(2:size(simvars_zeta)),err_mu, err_K)
        cov_cg_output   = covariance(C_growth,simvars%output(2:size(simvars%output)),err_mu, err_K)
        cov_cg_invest   = covariance(C_growth,simvars%invest(2:size(simvars%invest)),err_mu, err_K)
        cov_cg_ex_ret   = covariance(C_growth,excess_return,err_mu, err_K)
        cov_delta_r     = covariance(simvars_delta,simvars%r,err_mu, err_K)
        cov_delta_w     = covariance(simvars_delta,simvars%wage,err_mu, err_K)
        cov_delta_zeta  = covariance(simvars_delta,simvars_zeta,err_mu, err_K)
        cov_delta_output= covariance(simvars_delta,simvars%output,err_mu, err_K)
        cov_delta_invest= covariance(simvars_delta,simvars%invest,err_mu, err_K)
        cov_delta_ex_ret= covariance(simvars_delta,excess_return,err_mu, err_K)

        deallocate(simvars_zeta, simvars_delta)

        count_err_K    = count(err_K(lb:))
        count_err_mu   = count(err_mu(lb:))
        percent_err_K  = real(count_err_K,dp) / real(size(simvars%err_K(lb:)),dp) *100_dp
        percent_err_mu = real(count_err_mu,dp)/ real(size(simvars%err_mu(lb:)),dp)*100_dp
        end associate
       ! The next holds for mean shock only if wm=0.25, which is ususally true
        nk  = size(pol%apgrid,5)
        nmu = size(pol%apgrid,6)
        do jc=1,nj
            do xc =1,nx
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
    character(len=*) ,intent(in) :: prec
    character(:), allocatable :: fmt1
    integer :: nl, show_digits

    nl = K%get_namelength()
    show_digits = K%get_digits2display(prec)
    open(unit=21, file=path//'/equilibrium_'//prec//'.txt', status = 'replace', action='write')
    write(21,'(a12,a,",",a13,a)') ' Project    ', projectname, ' calibration ', calib_name
    write(21,*)
    write(21,*) 'Aggregate statistics'
    write(21,*) repeat('-',63)

    if (prec=='long') then
    write(21,123) repeat(' ',nl+2),' Average', 'Std.Dev.', 'Coef.Var', 'Autocorr', 'Abs.Max.', 'Abs.Min.'
    else
    write(21,123) repeat(' ',nl+2),' Average', 'Std.Dev.', 'Coef.Var', 'Autocorr'
    endif
123 format(a<nl+2>,a<show_digits+6>,5(<show_digits+2>x,a8))
    fmt1 = K%writing_format(show_digits,nl) ! To guarantee compatibility when call K%write(21,'K',prec)

    call K%write(21,'K',prec)   ! in units of efficient labor, so that stationary
    call output%write(21,'output',prec)
    call invest%write(21,'investm',prec)
    call mu%write(21,'mu',prec)
    call ex_ret%write(21,'Ex_ret',prec)
    call r%write(21,'r',prec)
    call rf%write(21,'rf',prec)
    if (prec == 'long') call r_pf_median%write(21,'rpf_med',prec)
    if (prec == 'long') call r_pf_kappa_med%write(21,'rpf_kapm',prec)
    call cons%write(21,'cons',prec)   ! in units of efficient labor, so that stationary
    call cg%write(21,'cons_grow',prec)
    call netwage%write(21,'netwage',prec)
    call pens%write(21,'pension',prec)
    call tau%write(21,'tau',prec)
    call zeta%write(21,'zeta',prec)
    call K_Y%write(21,'K_Y',prec)
    call welfare%write(21,'welfare',prec)

    if (prec == 'long') then
        write(21,*)
        call stock%write(21,'stock',prec)
        call bonds%write(21,'bonds',prec)
    endif

    write(21,*)
    write(21,*) 'Covariances'
    write(21,*) repeat('-',63)
    write(21,123) '            ',   '       r',   ' netwage',       '    zeta',      '  output',         '  invest',     '  Ex_ret'
    write(21,fmt1)' r         ' ,   r%std_()**2
    write(21,fmt1)' netwage   ' ,     cov_w_r , netwage%std_()**2
    write(21,fmt1)' zeta      ' ,  cov_zeta_r ,  cov_zeta_w , zeta%std_()**2
    write(21,fmt1)' output    ' ,cov_output_r ,cov_output_w , cov_output_zeta , output%std_()**2
    write(21,fmt1)' invest    ' ,   cov_inv_r ,   cov_inv_w ,    cov_inv_zeta ,  cov_inv_output , invest%std_()**2
    write(21,fmt1)' cons      ' ,  cov_cons_r ,  cov_cons_w ,   cov_cons_zeta , cov_cons_output ,  cov_cons_invest , cov_cons_ex_ret
    write(21,fmt1)' cons_grow ' ,    cov_cg_r ,    cov_cg_w ,     cov_cg_zeta ,   cov_cg_output ,    cov_cg_invest ,   cov_cg_ex_ret
    write(21,fmt1)' delta     ' , cov_delta_r , cov_delta_w ,  cov_delta_zeta ,cov_delta_output , cov_delta_invest ,cov_delta_ex_ret

    write(21,*)
    write(21,*) 'Laws of motion'
    write(21,*) repeat('-',63)
    if (n_coeffs == 2) then
        if (loms_in_logs) then
            write(21,123) '            ',   'constant',   '  log(k)',   '   R^2  '
        else
            write(21,123) '            ',   'constant',   '      k ',   '   R^2  '
        endif
    elseif (n_coeffs == 3) then
        if (loms_in_logs) then
            write(21,123) '            ',   'constant',   '  log(k)',   '      mu',   '   R^2  '
        else
            write(21,123) '            ',   'constant',   '      k ',   '      mu',   '   R^2  '
        endif
    endif

167 format(a<nl+1>,<n_coeffs>(es<show_digits+7>.<show_digits>,3x),'|',f<show_digits+4>.<show_digits+2>)
    if (loms_in_logs) then
        write(21,167)    " log(k')   "    , coeffs%k(:,1), coeffs%r_squared(1,1)
    else
        write(21,167)    " k'        "    , coeffs%k(:,1), coeffs%r_squared(1,1)
    endif
    do i=2,nz
        write(21,167) repeat(' ',nl+1), coeffs%k(:,i), coeffs%r_squared(1,i)
    enddo
    write(21,*)
    write(21,167)    ' mu        '    , coeffs%mu(:,1), coeffs%r_squared(2,1)
    do i=2,nz
        write(21,167) repeat(' ',nl+1), coeffs%mu(:,i), coeffs%r_squared(2,i)
    enddo

    write(21,*)
    write(21,*) 'Lifecycle statistics   (Cond. = conditional on survival)' ! in per capita units
    write(21,*) repeat('-',63)
    write(21,123) repeat(' ',nl+2),' Average', 'Std.Dev.', 'Cond.Avg', 'Cond.Std'
    write(21,fmt1) ' R_portf   ', sum(lc%return)/nj, sum(lc%return_var)/nj,  & !!! Attn: the second is Rpvar_avg!!!
                    sum(lc%return/pop_frac*pop_frac(1))/nj, sum(lc%return_var/pop_frac*pop_frac(1))/nj
    write(21,fmt1) ' Cons      ', sum(lc%cons)/nj, sum(lc%cons_var)/nj,  & !!! Attn: the second is Rpvar_avg!!!
                    sum(lc%cons/pop_frac*pop_frac(1))/nj, sum(lc%cons_var/pop_frac*pop_frac(1))/nj
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
    if (err%not_converged) write(21,*)' WARNING: root finder did not converge'


    write(21,*)   'Warnings in simulation'
    write(21,125) ' K     ', count_err_K, '  (',percent_err_K ,'%)   ', &
                  '    mu ', count_err_mu,'  (',percent_err_mu,'%)   '

    write(21,*)   'Duration measures'
    write(21,'(a7,i6, 14x, a5, es9.2)') ' KS_it ', it, 'secs ', secs

    close(21)
    end subroutine write_stats
!-------------------------------------------------------------------------------

    subroutine write_output_files()
    ! Since all arrays follow Fortran's natural storage order (i.e. column-major),
    ! no (implied) do-loops necessary when writing to file.
    ! The following will write all nx in one line, then change zc, then jc, ...
201 format(<nx> (f10.6,1x))
!    open(20, file=path//'/cons.txt', status = 'replace')
!    write(20,201) cons
!    close(20)

    open(40, file=path//'/apgrid.txt',  status = 'replace')
    write(40,201) pol%apgrid
    close(40)

    open(50, file=path//'/xgrid.txt',   status = 'replace')
    write(50,201) pol%xgrid
    close(50)

    open(60, file=path//'/kappa.txt',   status = 'replace')
    write(60,160) pol%kappa
160 format(<nx> (f0.6,x))
    close(60)

    open(40, file=path//'/apgrid_mean.txt',  status = 'replace')
    write(40,201) apgrid_mean
    close(40)

    open(50, file=path//'/xgrid_mean.txt',   status = 'replace')
    write(50,201) xgrid_mean
    close(50)

!    open(50, file=path//'/cons_mean.txt',   status = 'replace')
!    write(50,201) cons_mean
!    close(50)

    open(60, file=path//'/kappa_mean.txt',   status = 'replace')
    write(60,160) kappa_mean
    close(60)

    open(60, file=path//'/stocks_mean.txt',   status = 'replace')
    write(60,160) stocks_mean
    close(60)

    open(unit=21, file=path//'/Phi_tilde.txt', status = 'replace')
    write(21,220) sum(Phi,2) ! This is rough approximation, since different xgrids for each eta.
220 format(<nx> (es13.6,1x))
    close(21)

    open(unit=21, file=path//'/Phi.txt', status = 'replace')
    write(21,220) Phi
    close(21)

    open(unit=21, file=path//'/ap_lc.txt', status = 'replace')
    write(21,301) lc%ap
301 format(<nj> (f10.6,1x))
    close(21)

    open(unit=21, file=path//'/kappa_lc.txt', status = 'replace')
    write(21,360) lc%kappa
360 format(<nj> (f0.6,x))
    close(21)

    open(unit=21, file=path//'/cons_lc.txt', status = 'replace')
    write(21,301) lc%cons
    close(21)

    open(unit=21, file=path//'/stock_lc.txt', status = 'replace')
    write(21,301) lc%stock
    close(21)

    open(unit=21, file=path//'/consvar_lc.txt', status = 'replace')
    write(21,301) lc%cons_var
    close(21)

    open(unit=21, file=path//'/return_lc.txt', status = 'replace')
    write(21,301) lc%return
    close(21)

    open(unit=21, file=path//'/return_var_lc.txt', status = 'replace')
    write(21,302) lc%return_var
302 format(<nj> (es9.2,1x))
    close(21)

    open(unit=21, file=path//'/pop_frac.txt', status = 'replace')
    write(21,302) pop_frac
    close(21)

    open(unit=21, file=path//'/agg_grids.txt', status = 'replace')
    write(21,*) ' grid%k    =  '
    write(21,201) grids%k
    write(21,*) ' grid%mu   =  '
    write(21,201) grids%mu
    close(21)

    call err%write2file(path)

    open(unit=21, file=path//'/simvars.txt', status = 'replace')
    write(21,369) ' z        ', simvars%z
369 format(a10,<size(simvars%z)> (i2,8x))
    write(21,371) ' K        ', simvars%K
371 format(a10,<size(simvars%K)> (f9.6,1x))
    write(21,370) ' mu       ', simvars%mu
370 format(a10,<size(simvars%z)> (f9.6,1x))
    write(21,370) ' output   ', simvars%output
    write(21,370) ' invest   ', simvars%invest
    write(21,370) ' stock    ', simvars%stock
    write(21,370) ' bonds    ', simvars%bonds
    write(21,370) ' cons     ', simvars%C
    write(21,370) ' r        ', simvars%r
    write(21,371) ' rf       ', simvars%rf
    write(21,370) ' wage     ', simvars%wage
    write(21,370) ' pens     ', simvars%pens
    write(21,370) ' tau      ', simvars%tau
    write(21,373) ' welf     ', simvars%welf
    write(21,372) ' err_K    ', simvars%err_K
    write(21,372) ' err_mu   ', simvars%err_mu
372 format(a10,<size(simvars%z)> (1x,l1,8x))
    write(21,373) ' Phi_1    ', simvars%Phi_1
    write(21,373) ' Phi_nx   ', simvars%Phi_nx
    write(21,370) ' B        ', simvars%B
    write(21,373) ' err_inc  ', simvars%err_income
    write(21,373) ' err_aggr ', simvars%err_aggr
    write(21,373) ' bequests ', simvars%bequests
373 format(a10,<size(simvars%z)> (es9.2,1x))
    close(21)

    end subroutine write_output_files

end subroutine save_results

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine plot_results(path,mfile)
    ! Call Matlab and plot
    use ifport            ,only: system  ! Intel Fortran portability library
    character(len=*), intent(in)   :: path
    character(len=*), intent(in)   :: mfile
    integer                        :: err_matl
    character(:), allocatable      :: epssdir
    !print*, ' '
    !print*, 'Starting MATLAB to plot ',mfile,'(',dir,')'
    err_matl = system('mkdir '//path//'/graphs')
    epssdir = path(1+scan(path,'/'):)
    err_matl = system('EPSSDIR='//epssdir//' && export EPSSDIR && cd Matlab && matlab -nodesktop -nosplash -r '//mfile//' > cl_output.txt')
    if (err_matl ==-1) print*, 'Warning in save_results:plot: An error occured while plotting graphs'
end subroutine plot_results

end module save_results_mod
