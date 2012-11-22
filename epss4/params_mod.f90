module params_mod
! Contains all economic and numerical parameters.
! All 'variables' in this module have either parameter or protected attribute.

	use kinds
	use aggregate_grids_class, only: tAggGrids
	implicit none
!-------------------------------------------------------------------------------------------------
! The following are set and explained in the calibration file (see select_calibration_here.txt)
!-------------------------------------------------------------------------------------------------
	real(dp),protected :: theta, psi, beta, alpha, g, de_ratio, zeta_mean, zeta_std, del_mean, del_std,&
	                      pi1_zeta, pi1_delta, nu_sigma_h, nu_sigma_l, rho, n, tau, scale_AR, scale_IR, &
	                      factor_k, factor_mu, cover_k, cover_mu, apmax_factor, cmin, kappamax, &
	                      apmax_curv, tol_coeffs, tol_asset_eul, r_ms_guess, mu_ms_guess
    integer ,protected :: nj, jr, econ_life_start, nap, n_eta, n_zeta, n_delta, nk, nmu,&
                          n_coeffs, nt, t_scrap, nx_factor, opt_initial_ms_guess, run_n_times, run_counter_start
    logical ,protected :: ccv, surv_rates, def_contrib, partial_equilibrium, twosided_experiment, collateral_constraint, kappa_in_01,&
                          loms_in_logs, pooled_regression, estimate_from_simvars, exogenous_xgrid, &
                          save_all_iterations, detailed_euler_errs, normalize_coeffs, opt_zbren, opt_zbrak, tau_experiment

!-------------------------------------------------------------------------------------------------
! The following are 'derived parameters' calculated from the previous (read) parameters
!-------------------------------------------------------------------------------------------------
	! (Some more) demographics and pension system
    real(dp) ,dimension(:) ,allocatable ,protected :: &
            pop_frac,       &   ! Mass of generation as fraction of pop size N
            surv,           &   ! age-specific conditional survival rates
            ej                  ! age-specific labor productivity
    real(dp) ,protected ::  &
            P_L_ratio,      &   ! economic dependency ratio (Pensioners / Labor)
            L_N_ratio,      &   ! Labor / Population
            def_benefits,   &   ! defined benefits (amount)
            gamm                ! Recursive utility parameter: gamm=(1.0-theta)/(1.0-1.0/psi)

	! Grids (number of discrete grid points)
    integer ,protected ::   &
            nx,             &   ! grid for Cash-at-Hand, nx = nap (= grid for savings aprime)
            nz                  ! Number of aggregate states, nz = n_zeta*n_delta
    type(tAggGrids) ,protected :: ms_guess      ! Takes the values depending on opt_initial_ms_guess, see subroutine SetRemainingParams() below
    real(dp) ,dimension(:,:,:) ,allocatable ,protected :: apmax ! maximum a prime (savings) for each generation

	! Idiosyncratic and aggregate uncertainty
    real(dp) ,dimension(:,:) ,allocatable ,protected :: &
            pi_eta,         &   ! Markov transition matrix for idiosyncratic state
            etagrid,        &   ! idiosyncratic stochastic state
            pi_z                ! Markov transition matrix for aggregate shocks: pi(z1,z2)=Prob(z2|z1)
    real(dp) ,dimension(:)   ,allocatable ,protected :: &
            zeta,           &   ! aggregate technology shocks
            delta,          &   ! Shock to depreciation - REVERSE ORDER so that zc=1 is lowest return, zc=4 highest
            stat_dist_eta,  &   ! stationary distribution of idiosyncratic state
            stat_dist_z         ! stationary distribution of aggregate shock z
    integer, dimension(:)    ,allocatable ,protected :: &
            seed                ! seed to keep the same random variables over multiple calibration runs
    real(dp), protected ::  &
            scale_IR_orig, scale_AR_orig ! keep the original calibration values

interface params_set
    module procedure params_set_real, params_set_integer, params_set_logical
end interface params_set

! maybe implement the F2008 submodule feature: put the following, very large procedures in separate submodules
contains
!-------------------------------------------------------------------------------------------------
! - subroutine ReadCalibration()
! - subroutine SetRemainingParams() contains the following internal subroutines
!   -- subroutine set_demographics()
!   -- pure function set_pi_z(pi1_zeta, p1_delta)
!   -- subroutine set_idiosync_shocks(etagrid, pi_eta, stat_dist_eta, n_eta, ccv)
!   -- subroutine set_benefits
!   -- pure subroutine set_ms_guess(ms_guess, r_ms_guess, ccv, tau)
! - subroutine set_apmax(k, factor, curv_o)
! - subroutine CheckParams(zeta,delta)
! - pure function cal_id(calib_name)
! - pure function construct_path(dir, calib_name)
! - subroutine SaveParams(projectname)
! - subroutine params_set_real(param_name, new_value)
! - subroutine params_set_logical(param_name, new_value)
!-------------------------------------------------------------------------------------------------

subroutine SetDefaultValues()
! Some variables might be set that are not explicitly mentioned in the calibration file
    ! Reals
    theta=8.0; psi=0.5_dp; beta=0.98_dp; alpha=0.33_dp; g=0.01_dp; de_ratio=0.66_dp; zeta_mean=1.0; zeta_std=0.02_dp; del_mean=0.06_dp; del_std=0.06_dp
    pi1_zeta=0.7_dp; pi1_delta=.5_dp; nu_sigma_h=0.211_dp; nu_sigma_l=0.125_dp; rho=0.952_dp; n=0.01_dp; tau=0.0; scale_AR=0.0; scale_IR = 0.0
    factor_k=1.1_dp; factor_mu=1.1_dp; cover_k=0.8_dp; cover_mu=0.7_dp; apmax_factor=18.0_dp; cmin=1.0e-6_dp; kappamax=1000.0_dp
    apmax_curv=1.0; tol_coeffs=1e-4_dp; tol_asset_eul=1e-8_dp; r_ms_guess=3.0e-3_dp; mu_ms_guess=1.9e-2_dp
    ! Integers
    nj=64; jr=44; econ_life_start=22; nap=20; n_eta=2; n_zeta=2; n_delta=2; nk=10; nmu=8; n_coeffs=3; nt=5000; nx_factor=1; t_scrap=nt/10; opt_initial_ms_guess=0
    run_n_times=1; run_counter_start=1
    ! Logicals
    ccv=.true.; surv_rates=.false.; def_contrib=.true.; partial_equilibrium=.false.; twosided_experiment=.false.; collateral_constraint=.false.; kappa_in_01=.false.
    loms_in_logs=.true.; pooled_regression=.false.; estimate_from_simvars=.true.; exogenous_xgrid=.true.
    save_all_iterations=.false.; detailed_euler_errs=.false.; normalize_coeffs=.false.; opt_zbren=.true.; opt_zbrak=.false.; tau_experiment=.false.
end subroutine SetDefaultValues

subroutine ReadCalibration(calib_name)
    character(len=*) ,intent(in) :: calib_name
    character(len=80) :: parname, parval, parval2
    integer           :: istat, line, pos

    open (unit=21, file='model_input/'//calib_name, status='OLD', action='READ', iostat=istat)
    openif: if (istat==0) then
        line=0
        do
            line=line+1 ! blank lines in txt-file are not counted
            read (21,*,iostat=istat) parname, parval, parval2
            if (istat/=0) exit
            if (scan(parname,'!')>0) cycle

            select case (parname)
            case ('theta'); read (parval,*) theta
            case ('psi')
                pos = scan(parval,'/')
                if (pos>0) then
                    read (parval(:pos-1),*) psi
                    psi = psi / theta
                else
                    read (parval,*) psi
                endif
            case ('beta')
                read (parval,*) beta
            case ('alpha')
                read (parval,*) alpha
            case ('g')
                read (parval,*) g
            case ('de_ratio')
                read (parval,*) de_ratio
            case ('zeta_mean')
                read (parval,*) zeta_mean
            case ('zeta_std')
                read (parval,*) zeta_std
            case ('del_mean')
                read (parval,*) del_mean
            case ('del_std')
                read (parval,*) del_std
            case ('pi1_zeta')
                read (parval,*) pi1_zeta
            case ('pi1_delta')
                read (parval,*) pi1_delta
            case ('rho')
                read (parval,*) rho
            case ('nu_sigma_h')
                read (parval,*) nu_sigma_h
            case ('nu_sigma_l')
                read (parval,*) nu_sigma_l
            case ('nj')
                read (parval,*) nj
            case ('jr')
                read (parval,*) jr
            case ('econ_life_start')
                read (parval,*) econ_life_start
            case ('n')
                read (parval,*) n
            case ('tau')
                read (parval,*) tau
            case ('tau_experiment')
                read (parval,*) tau_experiment
            case ('ccv')
                read (parval,*) ccv
            case ('scale_AR')
                read (parval,*) scale_AR
            case ('scale_IR')
                read (parval,*) scale_IR
            case ('run_n_times')
                read (parval,*) run_n_times
            case ('run_counter_start')
                read (parval,*) run_counter_start
            case ('surv_rates')
                read (parval,*) surv_rates
            case ('def_contrib')
                read (parval,*) def_contrib
            case ('partial_equilibrium')
                read (parval,*) partial_equilibrium
            case ('twosided_experiment')
                read (parval,*) twosided_experiment
            case ('collateral_constraint')
                read (parval,*) collateral_constraint
            case ('kappa_in_01')
                read (parval,*) kappa_in_01
            case ('nap')
                read (parval,*) nap
            case ('n_eta')
                read (parval,*) n_eta
            case ('n_zeta')
                read (parval,*) n_zeta
            case ('n_delta')
                read (parval,*) n_delta
            case ('nk')
                read (parval,*) nk
            case ('nmu')
                read (parval,*) nmu
            case ('n_coeffs')
                read (parval,*) n_coeffs
            case ('nt')
                read (parval,*) nt
            case ('t_scrap')
                read (parval,*) t_scrap
            case ('nx_factor')
                read (parval,*) nx_factor
            case ('loms_in_logs')
                read (parval,*) loms_in_logs
            case ('pooled_regression')
                read (parval,*) pooled_regression
            case ('estimate_from_simvars')
                read (parval,*) estimate_from_simvars
            case ('opt_initial_ms_guess')
                read (parval,*) opt_initial_ms_guess
            case ('exogenous_xgrid')
                read (parval,*) exogenous_xgrid
            case ('save_all_iterations')
                read (parval,*) save_all_iterations
            case ('detailed_euler_errs')
                read (parval,*) detailed_euler_errs
            case ('opt_zbren')
                read (parval,*) opt_zbren
            case ('opt_zbrak')
                read (parval,*) opt_zbrak
            case ('factor_k')
                read (parval,*) factor_k
            case ('factor_mu')
                read (parval,*) factor_mu
            case ('cover_k')
                read (parval,*) cover_k
            case ('cover_mu')
                read (parval,*) cover_mu
            case ('apmax_factor')
                read (parval,*) apmax_factor
            case ('cmin')
                read (parval,*) cmin
            case ('kappamax')
                read (parval,*) kappamax
            case ('apmax_curv')
                read (parval,*) apmax_curv
            case ('tol_coeffs')
                read (parval,*) tol_coeffs
            case ('tol_asset_eul')
                read (parval,*) tol_asset_eul
            case ('r_ms_guess')
                read (parval,*) r_ms_guess
            case ('mu_ms_guess')
                read (parval,*) mu_ms_guess
            case default
                print '(a,a)', 'Unknown parameter: ',parname
            end select
        end do
        if (istat>0) print '(a,i6)', 'An error occured reading line', line
    else openif
        print '(a,i6)', 'ERROR opening calibration file: IOSTAT=',istat
        stop '*********STOP********* in Params:ReadCalib'
    end if openif
    close (unit=21)

end subroutine ReadCalibration
!-------------------------------------------------------------------------------------------------

subroutine SetRemainingParams()
! This is called from main.f90 right in the beginning
    use markov_station_distr
    integer, parameter :: iounit=124
    integer :: seedsize

    gamm = (1.0-theta)/(1.0-1.0/psi)
    nz   = n_zeta*n_delta
    nx   = nap
    scale_IR_orig= scale_IR
    if (scale_IR .ne. -1.0) scale_IR = 0.0 ! for the first run of a calibration
    scale_AR_orig= scale_AR
    if (scale_AR .ne. -1.0) scale_AR = 0.0 ! for the first run of a calibration

    if (.not.allocated(seed)) then ! this is done only once even in multiple calibration runs
        call random_seed(size=seedsize) ! the reason is that we want the same realization of variables over all runs
        allocate(seed(seedsize))
        call random_seed(get=seed)
    endif

    call set_demographics()

    pi_z = set_pi_z(pi1_zeta,pi1_delta, n_zeta, nz)

    stat_dist_z= f_markov_statdist(pi_z,400)       ! need to check eigenvalue calc!

    if (.not. def_contrib) call set_defined_benefits()

    select case (opt_initial_ms_guess)
    case (0) ! use previous ms equilibrium values saved in ./input/last_results/
        call ms_guess%read_unformatted('ms')
    case (1) ! use parameter-sensitive hard-coded guesses
        call set_ms_guess(ms_guess, r_ms_guess, ccv, scale_IR, tau)
    case (2) ! use user-supplied guess in this calibration file
        ms_guess%mu = [mu_ms_guess]
        ms_guess%k  = [(alpha/(r_ms_guess+del_mean))**(1.0/(1.0-alpha))] ! should I set ms_guess%k directly instead of r_ms_guess?
    case default
        print '(a,i1)', 'Warning in params_mod:SetRemainingParams(): opt_initial_ms_guess has value ', opt_initial_ms_guess
        print *, 'Admissible values are 0, 1, 2. Setting to 1 (i.e. use parameter-sensitive hard-coded guesses)'
        call set_ms_guess(ms_guess, r_ms_guess, ccv, scale_IR, tau)
    end select

contains

    subroutine set_demographics()
    ! Sets productivity profile, survival rates, and calculates pop ratios

        real(dp),allocatable,dimension(:,:) :: age_prod_profile, cond_mort_rates
        real(dp),allocatable,dimension(:)   :: mass_j ! Mass of generation
        real(dp)          :: P, L , Pop ! Pensioners, Labor (efficiency units), Total population
        integer           :: jc

        if (allocated(ej)) deallocate(ej,surv)
        allocate(ej(nj), surv(nj), mass_j(nj))
        allocate(age_prod_profile(2,nj+econ_life_start), cond_mort_rates(2,nj+econ_life_start))

        ! input of conditional mortality rates
        if (surv_rates) then
            open(iounit,file='model_input/data/mortality_rates_US.txt', action='read') ! starts at age 0
            read(iounit,*) cond_mort_rates
            close(iounit)
            surv(1)    = 1.0    ! agents born into model survive birth
            surv(2:nj) = 1.0 - cond_mort_rates(2,econ_life_start + 2 : econ_life_start + nj)  ! coz cond_mort_rates(econ_life_start) is for being born
        else
            surv       = 1.0
        endif

        ! input of age-earnings profiles of Hugget/Ventura/Yaron
        open(iounit,file='model_input/data/hvyageearn.txt', action='read')
        read(iounit,*) age_prod_profile
        close(iounit)

        ! Translate productivity profiles into model
        ej=0.0
        do jc=1,jr-1
            ej(jc)  = age_prod_profile(2, jc + econ_life_start) ! because age_prod_profile(econ_life_start) would correspond to j=0
        end do
        ! Normalization, so that centered around 1 and sum to jr-1
        ej=ej/sum(ej)*real(jr-1,dp)

        deallocate(age_prod_profile, cond_mort_rates)

        ! Calculate population ratios and fractions
        mass_j(1)=(1.0 + n)**(nj-1)
        do jc=2,nj
            mass_j(jc) = mass_j(jc-1)/(1.0 + n) * surv(jc)
        enddo
        Pop = sum(mass_j)                       ! total Population size
        L   = sum(mass_j(1:jr-1)*ej(1:jr-1))    ! Labor (in efficiency units)
        P   = sum(mass_j(jr:nj))                ! Pensioners

        P_L_ratio = P/L
        L_N_ratio = L/Pop
        pop_frac  = mass_j/Pop

    end subroutine set_demographics
!-------------------------------------------------------------------------------------------------

    pure function set_pi_z(pi1_zeta, p1_delta, n_zeta, nz)
        use fun_kronprod
        real(dp) ,dimension(nz,nz):: set_pi_z
        real(dp) ,intent(in)      :: pi1_zeta, p1_delta
        integer  ,intent(in)      :: n_zeta, nz
        real(dp)                  :: pi_zeta(n_zeta,n_zeta), pi_delta(n_delta,n_delta)
        real(dp)                  :: pi_1, w_KK, pi_delta_vec(nz) ! only for Kubler Kruger spec
        real(dp) ,dimension(nz,nz):: pi_KK, id, pi_zeta_kron, pi_delta_kron
        integer                   :: i
        real(dp)                  :: unit_mat(n_zeta,n_zeta), unit_vec(nz)
        logical  ,parameter       :: KK_setup = .false. ! Set pi_z like in Kubler/Krueger 2006 AER

        unit_mat =1.0
        unit_vec = 1.0

        if (KK_setup) then
            w_KK = .78_dp !.15_dp !68_dp !
            pi_1=.12_dp !.22_dp !12_dp

            pi_KK(:,1)=pi_1
            pi_KK(:,2)=0.5_dp -pi_1
            pi_KK(:,3)=pi_KK(:,2)
            pi_KK(:,4)=pi_KK(:,1)

            id=0.0; do i=1,nz; id(i,i)=1.0; enddo
            set_pi_z=(1.0-w_KK)*pi_KK + w_kk*id

        else    ! This is our standard, Harenberg/Ludwig way: includes STY, GM, iid
            pi_zeta = reshape([ pi1_zeta,     1.0-pi1_zeta, &
                                1.0-pi1_zeta, pi1_zeta],&
                      shape(pi_zeta), order=[2,1])
            pi_zeta_kron = f_kronprod(pi_zeta, unit_mat)

            pi_delta_vec = [pi1_delta, 1.0 - pi1_delta, 1.0 - pi1_delta, pi1_delta]
            do i=1, nz
                pi_delta_kron(:,i) = pi_delta_vec
            enddo

            set_pi_z = pi_zeta_kron * transpose(pi_delta_kron)
        endif

    end function set_pi_z
!-------------------------------------------------------------------------------

    subroutine set_defined_benefits
        open(iounit, file='model_input/last_results/def_benefits.unformatted', form='unformatted', action='read')
        read(iounit) def_benefits
        close(iounit)
        print*, 'sub_setvars: setting def_benefits = ', def_benefits
    end subroutine set_defined_benefits
!-------------------------------------------------------------------------------------------------

    pure subroutine set_ms_guess(ms_guess, r_ms_guess, ccv, scale_IR, tau)
    ! set the guesses specifically for each calibration
    ! might need to distinguish surv_rates
    type(tAggGrids), intent(out) :: ms_guess
    real(dp), intent(out) :: r_ms_guess
    logical, intent(in) :: ccv
    real(dp), intent(in) :: scale_IR, tau

    call ms_guess%allocate(1,1) ! needed because assigning scalars, not vectors
ccvif:  if (ccv) then
        if (tau >= .10_dp) then
            if (del_std >= .15_dp) then
                r_ms_guess  =  7.83e-2 !7.7e-2!8.2e-2 !8.4e-2 !6.6e-2 !7.5e-2 !8.8e-2 !6.0e-2 !8.2e-2 !8.4e-2 !5.7e-2
                ms_guess%mu =4.8e-2 !4.5e-2 !3.3e-2 !5.6e-2 !6.6e-2 !5.6e-2 !4.8e-2 !3.2e-2 ! 5.2e-2 !5.1e-2 !2.5e-2
            elseif (del_std >= .10_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 4.1e-2
            elseif (del_std >= .05_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 4.0e-2
            elseif (del_std >= .01_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 0.1e-2
            else
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            endif

        elseif (tau >= .05_dp) then
            if (del_std >= .15_dp) then
                r_ms_guess  =  6.3e-2!8.2e-2 !8.4e-2 !6.6e-2 !7.5e-2 !8.8e-2 !6.0e-2 !8.2e-2 !8.4e-2 !5.7e-2
                ms_guess%mu =4.3e-2 !3.3e-2 !5.6e-2 !6.6e-2 !5.6e-2 !4.8e-2 !3.2e-2 ! 5.2e-2 !5.1e-2 !2.5e-2
            elseif (del_std >= .10_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 3.0e-2
            elseif (del_std >= .05_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 1.0e-2
            elseif (del_std >= .01_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 0.1e-2
            else
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            endif
        elseif (tau >= .02_dp) then
            if (del_std >= .15_dp) then
                r_ms_guess  = 6.6e-2 !8.2e-2 !8.4e-2 !6.6e-2 !7.5e-2 !8.8e-2 !6.0e-2 !8.2e-2 !8.4e-2 !5.7e-2
                ms_guess%mu =5.0e-2 !3.3e-2 !5.6e-2 !6.6e-2 !5.6e-2 !4.8e-2 !3.2e-2 ! 5.2e-2 !5.1e-2 !2.5e-2
            elseif (del_std >= .10_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 3.0e-2
            elseif (del_std >= .05_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 1.0e-2
            elseif (del_std >= .01_dp) then
                r_ms_guess  = 6.4e-2
                ms_guess%mu = 0.1e-2
            else
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            endif
        elseif (tau >= 0.0) then
            if (del_std >= .15_dp) then
                r_ms_guess  = 6.0e-2 !8.4e-2 !5.7e-2
                ms_guess%mu = 3.2e-2 !5.1e-2 !2.5e-2
            elseif (del_std == .06_dp) then
                r_ms_guess  = 7.4e-2
                ms_guess%mu = 1.1e-2
            elseif (del_std == .05_dp) then
                r_ms_guess  =  5.6e-2
                ms_guess%mu = 1.6e-2
            elseif (del_std == .04_dp) then
                r_ms_guess  =  8.0e-2
                ms_guess%mu = 1.2e-2
            elseif (del_std == .03_dp) then
                r_ms_guess  =  7.9e-2
                ms_guess%mu = 0.8e-2
            else
                r_ms_guess  =  3.83e-2
                ms_guess%mu = 2.23e-2
            endif
        else
            if (del_std == .06_dp) then
                r_ms_guess  = 3.83e-2
                ms_guess%mu = 2.23e-2
            elseif (del_std == .05_dp) then
                r_ms_guess  =  6.2e-2
                ms_guess%mu = 1.7e-2
            else
                r_ms_guess  =  3.83e-2
                ms_guess%mu = 2.23e-2
            endif
        endif

    else ccvif
        if (tau >= .02_dp) then
            if (del_std >= .15_dp) then
                r_ms_guess  = 8.7e-2 !5.7e-2
                ms_guess%mu = 4.7e-2 !2.5e-2
            elseif (del_std == .06_dp) then
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            elseif (del_std == .05_dp) then
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            else
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            endif
        elseif (tau >= 0.0) then
           if (del_std >= .15_dp) then
                r_ms_guess  = 8.5e-2 !5.7e-2
                ms_guess%mu = 5.0e-2 !2.5e-2
            elseif (del_std == .06_dp) then
                r_ms_guess  = 3.83e-2
                ms_guess%mu = 2.23e-2
            elseif (del_std == .05_dp) then
                r_ms_guess  =  6.2e-2
                ms_guess%mu = 1.7e-2
            else
                r_ms_guess  =  3.83e-2
                ms_guess%mu = 2.23e-2
            endif
        else
            if (del_std == .06_dp) then
                r_ms_guess  = 3.83e-2
                ms_guess%mu = 2.23e-2
            elseif (del_std == .05_dp) then
                r_ms_guess  =  6.2e-2
                ms_guess%mu = 1.7e-2
            else
                r_ms_guess  =  3.83e-2
                ms_guess%mu = 2.23e-2
            endif
        endif
    endif ccvif

scaleIR:if (scale_IR < -0.8) then
        if (tau >= .02_dp) then
            if (del_std >= .15_dp) then
                r_ms_guess  = 6.3e-2 !13.0e-2
                ms_guess%mu = 4.3e-2 !4.1e-2 !4.5e-2 !3.8e-2
            elseif (del_std == .06_dp) then
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            elseif (del_std == .05_dp) then
                r_ms_guess  = 1.78e-1      ! is this correct?
                ms_guess%mu = 2.0e-2       ! is this correct?
            else
                r_ms_guess  = 8.1e-2
                ms_guess%mu = 5.0e-2
            endif
        elseif (tau >= 0.0) then
            if (del_std >= .15_dp) then
                r_ms_guess  = 11.0e-2
                ms_guess%mu = 3.9e-2
            elseif (del_std == .06_dp) then
                r_ms_guess  = 7.08e-2
                ms_guess%mu = 2.28e-2
            elseif (del_std == .05_dp) then
                r_ms_guess  = 7.08e-2
                ms_guess%mu = 2.28e-2
            else
                r_ms_guess  = 7.08e-2
                ms_guess%mu = 2.28e-2
            endif
        else
            if (del_std == .06_dp) then
                r_ms_guess  = 3.83e-2
                ms_guess%mu = 2.23e-2
            elseif (del_std == .05_dp) then
                r_ms_guess  =  6.2e-2
                ms_guess%mu = 1.7e-2
            else
                r_ms_guess  =  3.83e-2
                ms_guess%mu = 2.23e-2
            endif
        endif
    endif scaleIR

    ms_guess%k= (alpha/(r_ms_guess+del_mean))**(1.0/(1.0-alpha)) ! should I set ms_guess%k directly instead of r_ms_guess?
!    ms_guess%k=3.2
    end subroutine set_ms_guess

end subroutine SetRemainingParams

!-------------------------------------------------------------------------------------------------
subroutine params_set_thisrun()
    real(dp) :: zeta_std_scaled, del_std_scaled

    zeta_std_scaled = zeta_std*(1.0 + scale_AR)
    del_std_scaled  = del_std *(1.0 + scale_AR)
    zeta =[zeta_mean-zeta_std_scaled, zeta_mean-zeta_std_scaled, zeta_mean+zeta_std_scaled, zeta_mean+zeta_std_scaled]
    delta=[del_mean + del_std_scaled, del_mean - del_std_scaled, del_mean + del_std_scaled, del_mean - del_std_scaled]

    call set_idiosync_shocks(etagrid, pi_eta, stat_dist_eta, n_eta, nz, ccv)

    call set_apmax(ms_guess%k(1), apmax_factor, scale_IR ,apmax_curv)

contains
!-------------------------------------------------------------------------------

    subroutine set_idiosync_shocks(etagrid, pi_eta, stat_dist_eta, n_eta, nz, ccv)
        use markov_chain_approx
        use twostate_exact

        integer  ,intent(in)  :: n_eta, nz
        logical  ,intent(in)  :: ccv
        real(dp) ,allocatable ,intent(out) :: etagrid(:,:), pi_eta(:,:), stat_dist_eta(:)
        real(dp) :: sigma(nz)
        integer  :: zc, ec
        logical  :: converged

        allocate(etagrid(n_eta,nz), pi_eta(n_eta,n_eta), stat_dist_eta(n_eta))

        if (ccv) then
            sigma = [nu_sigma_h,nu_sigma_h,nu_sigma_l,nu_sigma_l]
        else
            sigma = (nu_sigma_h + nu_sigma_l)/2.0
        endif

        sigma = sigma * (1.0 + scale_IR)

        if (n_eta == 2) then
            call calibridiorisk2(rho,sigma,pi_eta,etagrid,stat_dist_eta, converged)
            if (.not. converged) then
                print*, 'CRITICAL ERROR in params_mod:setremainingparams:set_idiosync_shocks'
                print*, 'could not set etagrid or pi_eta'
                stop
            endif
        else
            ! Using Tauchen not a good idea because different pi_eta in case ccv
            do zc = 1,nz
                call rouwenhorst(rho, sigma(zc), pi_eta, etagrid(:,zc),stat_dist_eta)
                etagrid(:,zc) = exp(etagrid(:,zc))
                etagrid(:,zc) = etagrid(:,zc)/dot_product((etagrid(:,zc)),stat_dist_eta(:))
            enddo
        endif

    end subroutine set_idiosync_shocks
end subroutine params_set_thisrun
!-------------------------------------------------------------------------------------------------

subroutine set_apmax(k, factor_o, scale_IR_o, curv_o)
! This is a very adhoc function, based on trials
    use makegrid_mod
    real(dp), intent(in)    :: k
    real(dp), intent(in), optional :: factor_o, scale_IR_o, curv_o
    real(dp) :: factor, scale_IR_loc, curv, &
                guess ,&             ! guess a constant as basis for min and max of grid
                shrink_factor_eta,&  ! shrink apmax in n_eta by shrink_factor percent to get apmax in etac = 1
                shrink_factor_z  ,&  ! shrink apmax in nz by shrink_factor percent to get apmax in zc = 1
                zc_fact, etac_fact
    integer  :: jmax, etac, zc

    if (present(factor_o)) then
        factor = factor_o
    else
        factor = apmax_factor
    endif
    if (present(scale_IR_o)) then
        scale_IR_loc = scale_IR_o
    else
        scale_IR_loc = scale_IR
    endif
    if (present(curv_o)) then
        curv = curv_o
    else
        curv = apmax_curv
    endif


    shrink_factor_eta = 0.02_dp*(1.0 + scale_IR_loc)  ! make dependent on realizations of eta (etagrid)?!

    if (zeta_std == 0.0) then  ! should I put scale_AR = -1.0 here?
        shrink_factor_z = 0.0
    else
        shrink_factor_z = 0.1_dp !0.25_dp  ! make dependent on realizations of z grid?!
    endif

    if (allocated(apmax)) deallocate(apmax)
    allocate(apmax(n_eta,nz,nj))
    jmax =jr-1

    guess         = k
    apmax(n_eta,nz,1)      = guess/2.0    ! calc wage?
    apmax(n_eta,nz,nj)     = guess*2.0    ! determins max cons of nj, and apmax(:,nj-1)
    apmax(n_eta,nz,jmax)   = guess*factor
    apmax(n_eta,nz,jmax:1:-1)     = -MakeGrid(-apmax(n_eta,nz,jmax),-apmax(n_eta,nz,1),jmax, 1.0/curv)
    apmax(n_eta,nz,jmax+1:nj-1)   = MakeGrid(apmax(n_eta,nz,jmax),apmax(n_eta,nz,nj),nj-jmax-1, curv) !curv, 2.0_dp

    do zc=1,nz
        zc_fact = (1.0 - shrink_factor_z/real(nz-1,dp)*real(nz -zc,dp))
	    do etac=1,n_eta
	       etac_fact = (1.0 - shrink_factor_eta/real(n_eta-1,dp)*real(n_eta -etac,dp))
	       apmax(etac,zc,:) = apmax(n_eta,nz,:) * etac_fact * zc_fact
	    enddo
    enddo
end subroutine set_apmax

!-------------------------------------------------------------------------------------------------
subroutine CheckParams()
use omp_lib           ,only: OMP_get_max_threads
! Very simple error checks on the parameters.
! Important because of different scenarios/calibrations, saved in modules calib_**.f90
! Put checks into the modules / objects they belong to? Would probably destroy pure in some cases?
    real(dp), parameter :: crit = 0.0000001_dp

    if (ms_guess%k(1) < 1e-6) then
        print*, 'ERROR: ms_guess%k(1)<1e-6'
        call critical_stop
    elseif (ms_guess%k(1)>1000.0) then
        print*, 'ERROR: ms_guess%k(1)>1000'
        call critical_stop
    elseif (ms_guess%k(1)>50.0) then
        print*, 'WARNING: ms_guess%k(1)>50'
    endif

    if (ms_guess%mu(1)<-1.0) then
        print*, 'ERROR: ms_guess%mu(1)<-1.0'
        call critical_stop
    elseif (ms_guess%mu(1) < 1e-4) then
        print*, 'WARNING: ms_guess%mu(1) <= 1e-4'
    elseif (ms_guess%mu(1)>1.0) then
        print*, 'ERROR: ms_guess%mu(1)>1.0'
        call critical_stop
    elseif (ms_guess%mu(1)>0.2_dp) then
        print*, 'WARNING: ms_guess%mu(1)>0.2'
    endif

    if (n_coeffs>3) then
        print*, 'ERROR: n_coeffs>3, check laws_of_motion:initialize_coeffs AND :Regression!'
        call critical_stop
    elseif (n_coeffs<2) then
        print*, 'ERROR: n_coeffs<2, check laws_of_motion:initialize_coeffs AND :Regression!'
        call critical_stop
    endif

    if (theta == 1.0) then
        print*, 'WARNING: theta == 1.0, Taylor expansion in household_solution_mod:asseteuler doesnt work!'
    endif

    if (nz /= 4) then
        print*, 'ERROR: nz /= 4'
        call critical_stop
    endif

    if (n_eta /= 2) then
        print*, 'Warning: n_eta /= 2'
    endif

    if (nt<=t_scrap) then
        print*, 'ERROR: nt<=t_scrap'
        call critical_stop
    endif

    if ((nt-t_scrap)*OMP_get_max_threads() < 5000) then
        print*, 'Warning: too few simulation periods? (nt-t_scrap)*OMP_get_max_threads() < 5000'
    endif

    if (nx_factor < 1) then
        print*, 'ERROR: nx_factor < 1'
        call critical_stop
    elseif (nx_factor > 100) then
        print*, 'WARNING: nx_factor > 100 : might run out of mem/ take long time'
    endif

    if (jr>nj+1) then
        print*, 'ERROR: jr > nj+1'
        call critical_stop
    endif

    if (n<0.0) then
        print*, 'WARNING: n < 0'
    endif

    if (any(surv>1.0) .or. any(surv<=0.0)) then
        print*, 'ERROR: surv > 1 or surv <= 0'
        call critical_stop
    endif

    if (any(ej<0.0)) then
        print*, 'ERROR: ej < 0'
        call critical_stop
    endif
    if (abs(sum(ej) - real(jr-1,dp)) > crit) print*, 'WARNING: sum(ej) .ne. jr-1'

    if (abs(sum(pop_frac)-1.0) > crit) then
        print*, 'ERROR: sum(pop_frac) \= 1'
        call critical_stop
    endif

    if (n_coeffs > 6) then
        print*, 'ERROR: n_coeffs > 6'
        print*, '       covariates not defined in sub_regression'
        call critical_stop
    elseif(n_coeffs < 2) then
        print*, 'ERROR: n_coeffs < 2'
        call critical_stop
    endif

    if (detailed_euler_errs) print*, 'WARNING: detailed_euler_errs=.true., but not implemented, check error_class'

    if (zeta_mean .ne. 1.0) print*, 'WARNING: zeta_mean .ne. 1.0'
    if (zeta_std > 0.1_dp) print*, 'WARNING: zeta_std > 0.1'
    if (zeta_mean - zeta_std <= 0.0) then
        print*, 'ERROR: zeta_mean - zeta_std <= 0.0'
        call critical_stop
    endif
    if (del_std  > 0.2_dp) print*, 'WARNING: del_std > 0.2'

    if (pi1_zeta > 1.0 .or. pi1_zeta < 0.0) then
        print*, 'ERROR: pi1_zeta not in [0.0, 1.0]'
        call critical_stop
    endif

    if (pi1_delta > 1.0 .or. pi1_delta <= 0.0) then
        print*, 'ERROR: pi1_delta not in (0.0, 1.0]'
        call critical_stop
    endif

    if (any(abs(stat_dist_z-1.0/real(nz,dp))>crit)) then
        print '(a,<nz>f7.4)', ' WARNING: stat_dist_z = ', stat_dist_z
    endif

    if (abs(dot_product(stat_dist_z,zeta)-1.0)>crit) then
        print*, 'WARNING: dot_product(stat_dist_z,zeta) .ne. 1'
        call critical_stop
    endif

    if (abs(dot_product(stat_dist_z,delta)-del_mean)>crit) then
        print*, 'WARNING: dot_product(stat_dist_z,delta) .ne. del_mean'
        call critical_stop
    endif

    if (abs(sum(stat_dist_eta) -1.0)>crit) then
        print*, 'Error: sum(stat_dist_eta) .ne. 1'
        call critical_stop
    endif

    if (n_eta == 2 .and. any(abs(stat_dist_eta-1.0/real(n_eta,dp))>crit)) then
        print '(a,<n_eta>f7.4)', ' WARNING: stat_dist_eta = ', stat_dist_eta
    endif

    if (any(abs(matmul(stat_dist_eta,etagrid)-1.0)>crit)) then
        print*, 'WARNING: matmul(stat_dist_eta,etagrid) .ne. 1 : ', matmul(stat_dist_eta,etagrid)
        call critical_stop
    endif

    if (r_ms_guess > .2_dp) then
        print*, 'ERROR: r_ms_guess > .2'
        call critical_stop
    endif

    if (ms_guess%mu(1) > .2_dp) then
        print*, 'ERROR: ms_guess%mu > .2'
        call critical_stop
    endif

    if (ms_guess%mu(1) < 0.0) then
        print*, 'ERROR: ms_guess%mu < 0.0'
        call critical_stop
    endif

    if (ms_guess%k(1) < 0.0) then
        print*, 'ERROR: ms_guess%k < 0.0'
        call critical_stop
    endif

    if (any(apmax > 1e6) .or. any(apmax < 0.0)) then
        print*, 'ERROR: some apmax unreasonably large or below zero'
        call critical_stop
    endif

    if (tau > 1.0 .or. tau < 0.0) then
        print*, 'ERROR: tau not in [0.0 , 1.0]'
        call critical_stop
    endif

    if (cover_k > 1.0 .or. cover_k < 0.0) then
        print*, 'ERROR: cover_k not in [0.0 , 1.0]'
        call critical_stop
    endif

    if (cover_mu > 1.0 .or. cover_mu < 0.0) then
        print*, 'ERROR: cover_mu not in [0.0 , 1.0]'
        call critical_stop
    endif

    if (psi > 5.0 .or. psi < 0.0) then
        print*, 'ERROR: psi out of range, psi = ', psi
        call critical_stop
    endif

    if (rho > 1.0 .or. rho < 0.0) then
        print*, 'ERROR: rho out of range, rho = ', rho
        call critical_stop
    endif

    if (scale_IR_orig .ne. 0.0 .and. scale_IR_orig .ne. -1.0 .and. scale_AR_orig .ne. 0.0 .and. scale_AR_orig .ne. -1.0) then
        print*, 'ERROR: scale_IR_orig .ne. 0.0 .and. scale_IR_orig .ne. -1.0 .and. scale_AR_orig .ne. 0.0  .and. scale_AR_orig .ne. -1.0'
        print*, '       -> Set either scale_IR to (-1,0) or scale_AR to (-1,0).'
        call critical_stop
    endif

    if (scale_IR > 1.0e2 .or. scale_IR < -1.0) then
        print*, 'ERROR: scale_IR out of range, scale_IR = ', scale_IR
        call critical_stop
    endif

    if (scale_AR > 1.0e2 .or. scale_AR < -1.0) then
        print*, 'ERROR: scale_AR out of range, scale_AR = ', scale_AR
        call critical_stop
    elseif (scale_AR > 5.0 .or. (scale_AR < -0.5_dp .and. scale_AR > -1.0)) then
        print*, 'WARNING: extreme value for scale_AR = ', scale_AR
    elseif (scale_AR == -1.0) then
        print*, 'ATTENTION: scale_AR = -1.0, i.e. no aggregate risk'
    endif

    if (run_counter_start < 1 .and. .not. twosided_experiment) then
        print*, 'Warning: run_counter_start < 1 .and. .not. twosided_experiment, setting to 1'
        run_counter_start = 1
    endif

    if (run_n_times <  run_counter_start) then
        print*, 'Warning: run_n_times < run_counter_start, setting to run_n_times = run_counter_start'
        run_n_times = run_counter_start
    endif

    if ((scale_IR_orig == 0.0 .or. scale_IR_orig == -1.0) .and. (scale_AR_orig == 0.0 .or. scale_AR_orig == -1.0) .and. run_n_times > 1 ) then
        print*, 'ERROR: scale_IR =0 or -1, scale_AR = 0 or -1, and run_n_times > 1.'
        print*, '       -> set run_n_times = 1 or specify experiment.'
        call critical_stop
    endif

    if (scale_AR_orig /= -1) then
	    if (scale_AR_orig*real(run_n_times-1,dp) > 1.0e2 .or. scale_AR_orig*real(run_n_times-1,dp) < -1.0) then
	        print*, 'ERROR: scale_AR * run_n_times out of range = ', scale_AR_orig*real(run_n_times-1,dp)
	        call critical_stop
	    endif
    endif

    if (scale_IR_orig /= -1) then
	    if (scale_IR_orig*real(run_n_times-1,dp) > 1.0e2 .or. scale_IR_orig*real(run_n_times-1,dp) < -1.0) then
	        print*, 'ERROR: scale_IR * (run_n_times-1) out of range = ', scale_IR_orig*real(run_n_times-1,dp)
	        call critical_stop
	    endif
    endif

    if ((abs(scale_IR_orig) < .01 .and. scale_IR_orig /= 0.0) .or. (abs(scale_AR_orig) < .01 .and. scale_AR_orig /= 0.0)) then
	    print*, 'ERROR: scale_IR or scale_AR < 0.01: folder structure not fine enough in main'
	    call critical_stop
    endif

    if (scale_AR == -1.0 .and. .not. partial_equilibrium) then
        print*, 'ERROR: scale_AR = -1.0, ie. no aggregate risk, not yet implemented with general equilibrium'
        call critical_stop
    endif

contains
    subroutine critical_stop()
        stop '*********STOP********* in params_mod:CheckParams'
    end subroutine critical_stop

end subroutine CheckParams

!-------------------------------------------------------------------------------------------------
pure function cal_id(calib_name)
	character(:), allocatable :: cal_id, temp
	character(len=*), intent(in) :: calib_name
	integer :: pos, pos2
    pos  = scan(calib_name,'/')
    temp = calib_name(pos+1:)
    pos  = scan(temp,'.')
    cal_id = temp(:pos-1)
    pos2 = index(temp,',')
    if (pos2 > pos) cal_id = cal_id//temp(pos2:)
end function cal_id
!-------------------------------------------------------------------------------------------------
pure function construct_path(dir, calib_name)
!	use params_mod, only: cal_id
	character(:), allocatable :: construct_path
	character(len=*), intent(in) :: dir, calib_name
    construct_path= 'model_output/'//cal_id(calib_name)//'/'//dir
end function construct_path

!-------------------------------------------------------------------------------------------------
subroutine SaveParams(projectname, calib_name)
! Save parameters and exogenous variables
    use omp_lib           ,only: OMP_get_max_threads
    character(len=*), intent(in)    :: projectname, calib_name
    character(:), allocatable :: path

    path = 'model_output/'//cal_id(calib_name)

    open(unit=21, file=path//'/params.txt', status = 'replace', action='write')
    write(21,'(a9,a,",",a13,a)') ' Project ', projectname, ' calibration ', calib_name
    write(21,*) '-------------------------- Prefs and Tech --------------------------'
    write(21,211) ' theta        = ', theta
211 format(a16, f0.4)
    write(21,211) ' psi          = ', psi
    write(21,211) ' gamma        = ', gamm
    write(21,211) ' beta         = ', beta
    write(21,211) ' alpha        = ', alpha
    write(21,211) ' de_ratio     = ', de_ratio
    write(21,211) ' g            = ', g
    write(21,211) ' zeta_mean    = ', zeta_mean
    write(21,211) ' zeta_std     = ', zeta_std
    write(21,211) ' del_mean     = ', del_mean
    write(21,211) ' del_std      = ', del_std
    write(21,*)
    write(21,*) '------------------- Soc. Sec. & Demographics -----------------------'
    write(21,211) ' tau          = ', tau
    write(21,211) ' def_benefits = ', def_benefits
    write(21,218) ' nj           = ', nj
    write(21,218) ' jr           = ', jr
    write(21,211) ' n            = ', n
    write(21,211) ' P_L_ratio    = ', P_L_ratio
    write(21,211) ' L_N_ratio    = ', L_N_ratio
    write(21,*)
    write(21,*) '----------------------- Stochastic processes -----------------------'
    write(21,212) ' pi_zeta(1,1) = ', pi1_zeta
    write(21,212) ' pi_delta(1,1)= ', pi1_delta
212 format(a16, f0.4)
    write(21,*)
    write(21,215) ' pi_z         = ', transpose(pi_z)
215 format(a16,<nz>(f0.4,x),/,<nz-1>(t17,<nz>(f0.4,x),/))
    write(21,'(a16,<nz>(f0.4,x))') ' stat_dist_z  = ', stat_dist_z
    write(21,*)
    write(21,209) ' pi_eta       = ', transpose(pi_eta)
209 format(a16,<n_eta>(f0.4,x),/,<n_eta-1>(t17,<n_eta>(f0.4,x),/))
    write(21,214) ' etagri       = ', transpose(etagrid)
214 format(a16,<nz>(f0.4,x),/,<n_eta-1>(t17,<nz>(f0.4,x),/))
    write(21,'(a16,<n_eta>(f0.4,x))') ' stat_dist_eta= ', stat_dist_eta
    write(21,*)
    write(21,*) '------------------------------ Grids -------------------------------'
    write(21,218) ' nap          = ', nap
218 format(a16, i3)
    write(21,218) ' nx           = ', nx
    write(21,218) ' n_eta        = ', n_eta
    write(21,218) ' nz           = ', nz
    write(21,218) ' nk           = ', nk
    write(21,218) ' nmu          = ', nmu
    write(21,'(a16, i5)') ' nt           = ', nt
    write(21,'(a16, i5)') ' tscrap       = ', t_scrap
    write(21,218) ' nx_factor    = ', nx_factor
    write(21,217) ' cover_k      = ', cover_k
    write(21,217) ' cover_mu     = ', cover_mu
217 format(a16, 2(f0.6,x))
    write(21,219) ' cmin         = ', cmin
219 format(a16, es8.2)
    write(21,218) ' No. threads  = ', OMP_get_max_threads()
    write(21,*)
    write(21,*) '----------------------------- Options ------------------------------'
    write(21,'(a20,l1)') ' ccv             =  ', ccv
    write(21,'(a20,l1)') ' def_contrib     =  ', def_contrib
    write(21,'(a20,l1)') ' collat_constr   =  ', collateral_constraint
    write(21,'(a20,l1)') ' kappa_in_01     =  ', kappa_in_01
    write(21,'(a20,l1)') ' exogenous_xgrid =  ', exogenous_xgrid
    write(21,'(a20,l1)') ' loms_in_logs    =  ', loms_in_logs
    write(21,'(a20,l1)') ' pooled_regress  =  ', pooled_regression
    write(21,'(a20,l1)') ' surv_rates      =  ', surv_rates
    write(21,'(a20,l1)') ' opt_zbren       =  ', opt_zbren
    write(21,'(a20,l1)') ' opt_zbrak       =  ', opt_zbrak
    write(21,216)        ' scale_AR        =  ', scale_AR
    write(21,216)        ' scale_IR        =  ', scale_IR
    write(21,*)
    write(21,*) '------------------------------ Guesses -----------------------------'
    write(21,216) ' r_ms_guess      =  ', r_ms_guess
216 format(a20, f0.4)
    write(21,216) ' mu_ms_guess     =  ', ms_guess%mu
    write(21,216) ' k_ms_guess      =  ', ms_guess%k
    write(21,216) ' apmax_factor    =  ', apmax_factor
    write(21,216) ' apmax_curv      =  ', apmax_curv
    write(21,*)
    write(21,*) '---------------------------- Tolerances ----------------------------'
    write(21,220) ' tol_asset_eul   =  ', tol_asset_eul
220 format(a20, es8.2)
    write(21,220) ' tol_coeffs      =  ', tol_coeffs
    write(21,*)
    write(21,*)
    write(21,*)
    write(21,*) '--------------------------------------------------------------------'
    write(21,*) '----------------------------- Comments -----------------------------'
    write(21,*) ' This file reports the effectively used values, e.g:'
    write(21,*) ' - option opt_initial_ms_guess selects how ms_guesses are chosen'
    write(21,*) ' - value of std(zeta) reported here is zet_sd = zeta_std*(1+scale_AR)'
    write(21,*) ' - value of collateral_constraint, which might not be available in'
    write(21,*) '     some calibration files and thus reverts to the default value.'
    write(21,*) ' This file helps as counter check for the various calibration files.'
    close(21)
end subroutine SaveParams


subroutine params_set_real(param_name, new_value)
! this is inefficient, should do one setter method for each
    character(len=*), intent(in) :: param_name
    real(dp), intent(in)     :: new_value
	    select case (param_name)
	    case ('scale_IR')
	       if (new_value < -1.0) then
	           print* , 'params_set: scale_IR new_value < -1.0, setting to -1.0'
	           scale_IR = -1.0
           else
	           scale_IR = new_value
           endif
        case ('scale_AR')
           if (new_value < -1.0) then
               print* , 'params_set: scale_AR new_value < -1.0, setting to -1.0'
               scale_AR = -1.0
           else
               scale_AR = new_value
           endif
        case ('tau')
            if (new_value < 0.0) then
                print* , 'params_set: tau new_value < 0.0, setting to 0.0'
                tau = 0.0
            elseif (new_value > 1.0) then
                print* , 'params_set: tau new_value > 1.0, setting to 1.0'
                tau = 1.0
            else
                tau= new_value
            endif
		case default
		    print '(a,a)', 'params_mod:params_set: Cannot set parameter ',param_name
		end select
end subroutine params_set_real

subroutine params_set_integer(param_name, new_value)
! this is inefficient, should do one setter method for each
    character(len=*), intent(in) :: param_name
    integer, intent(in)     :: new_value
        select case (param_name)
        case ('run_counter_start')
            run_counter_start = new_value
        case default
            print '(a,a)', 'params_mod:params_set: Cannot set parameter ',param_name
        end select
end subroutine params_set_integer

subroutine params_set_logical(param_name, new_value)
    character(len=*), intent(in) :: param_name
    logical, intent(in)     :: new_value
    select case (param_name)
        case ('partial_equilibrium')
            partial_equilibrium = new_value
    end select
end subroutine params_set_logical

end module params_mod
