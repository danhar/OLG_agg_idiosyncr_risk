!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

! main.f90 of program Equity Premium Social Security
program EPSS

    use ifport             ,only: system  ! Intel Fortran portability library
	use params_mod         ,only: SetDefaultValues,ReadCalibration, SetRemainingParams, CheckParams, cal_id, params_set, params_set_thisrun, welfare_decomposition, alt_insurance_calc, surv_rates, debugging,&
	                              n_end_params, run_n_times, run_counter_start, twosided_experiment, scale_AR, scale_IR, scale_AR_orig, scale_IR_orig, tau_experiment, tau, surv_rates, ccv, dp, calc_euler_errors, &
	                              tau_increment, tau_calib, tau_GE0, stockshare_fixed_orig
	use calibration_mod    ,only: calibrate
	use run_model_mod

    implicit none
	real(dp)                  :: secs, elapsed_time
	real(dp)     ,allocatable :: welfare(:,:), agg_cons(:,:), cev(:), agg_cons_ratio(:), cev_ins(:,:), welfare_ins(:,:,:), risk_scale(:)
	integer                   :: sys_error, rc, i
	integer(8)                :: start_time, end_time, count_rate, count_max ! to avoid maximum time, but relies on intel extension
	logical                   :: exit_main_loop
	character(:) ,allocatable :: projectname, calib_name, calib_name_base
	character(len=8)          :: runchar

    call system_clock(start_time)
!    sys_error = system('rm -fr model_output/*') ! Deletes existing output files

    call get_projectname(projectname)
    print*, ' '
   	print*, 'Starting program ', projectname
    print*, ' '

    do
        call SetDefaultValues
        call get_calibration_name(calib_name, exit_main_loop)
        if (exit_main_loop) exit

        ! copy the input to new, because that is the directory that will be worked with
        sys_error = system('cp -r model_input/last_results/'//cal_id(calib_name)//'/tau* model_input/last_results/'//cal_id(calib_name)//'/new/')
        calib_name_base = calib_name

	    print*, '- main: Reading calibration file '// calib_name
	    call ReadCalibration(trim(adjustl(calib_name)))
	    call SetRemainingParams(calib_name)

        ! First, calibrate the economy to targets (if applicable).
	    if (n_end_params > 0) then
	        print*, ' '
	        print*, '- main: Starting calibration routine'
            call params_set('tau', tau_calib)
	        call params_set_thisrun
            call CheckParams
            write(runchar,'(a4)') ',cal'
            calib_name = trim(calib_name)//trim(runchar)
            sys_error = system('mkdir model_output/'//cal_id(calib_name)) ! could create different folder with _cal attached?
	        call calibrate(projectname, calib_name)
	        call params_set('tau', tau_GE0)
        endif

        ! Second, set up experiments and run them.
        if (twosided_experiment .and. run_n_times>1) call params_set('run_counter_start', -1*run_n_times+2)
        if (allocated(welfare)) deallocate(welfare)
        allocate(welfare(run_counter_start:run_n_times,2))
        agg_cons = welfare
        if (allocated(risk_scale)) deallocate(risk_scale)
        allocate(risk_scale(run_counter_start:run_n_times))
        if (allocated(cev)) deallocate(cev)
        if (tau_experiment)  allocate(cev(run_counter_start:run_n_times))
        if (allocated(agg_cons_ratio)) deallocate(agg_cons_ratio)
        if (tau_experiment) agg_cons_ratio = cev
        if (alt_insurance_calc) then
            if (allocated(welfare_ins)) deallocate(welfare_ins)
            allocate(welfare_ins(lbound(welfare,1):ubound(welfare,1),lbound(welfare,2):ubound(welfare,2),5))
            if (allocated(cev_ins)) deallocate(cev_ins)
            if (tau_experiment)  allocate(cev_ins(lbound(cev,1):ubound(cev,1),lbound(welfare_ins,3):ubound(welfare_ins,3)))
        endif

        write (runchar, *) ' '
        do rc=run_counter_start,run_n_times ! always run one time without experiment
            if (run_counter_start/=run_n_times) then
                if (scale_IR_orig .ne. 0.0 .and. scale_IR_orig .ne. -1.0) then
                    call params_set('scale_IR', scale_IR_orig*real(rc-1,dp))
                    risk_scale(rc)= 1.0 + scale_IR
                    write(runchar,'(a3,f4.2)') ',IR', risk_scale(rc)
                elseif (scale_AR_orig .ne. 0.0 .and. scale_AR_orig .ne. -1.0) then
                    call params_set('scale_AR', scale_AR_orig*real(rc-1,dp))
                    risk_scale(rc)= 1.0 + scale_AR
                    write(runchar,'(a3,f4.2)') ',AR', risk_scale(rc)
                elseif (rc ==0) then
                    write(runchar,'(a4)') ',GE0'
                    call params_set('tau_experiment', .false.)
                    ! The following is necessary, because in this case we need to calibrate to the smaller tau so as to get consistent results
                    ! if (scale_AR == -1.0) call params_set('tau', tau+ tau_increment)
                elseif (rc ==1) then
                    write(runchar,'(a4)') ',GE1'
                    call params_set('tau_experiment', .true.)
                    call params_set('tau', tau- tau_increment) ! because we always calibrate to the higher tau
                    call params_set('partial_equilibrium', .false.) ! redundant since tau_experiment=.false. in previous GE0, but keep for safety.
                    call params_set('check_dynamic_efficiency', .true.)
                elseif (rc ==2) then ! the following are for the welfare decomposition
                    call params_set('surv_rates', .false.)
                    write(runchar,'(a7)') ',noSURV'
                    ! The current and following are partial equilibrium, so we don't calc Euler errors.
                    call params_set('check_dynamic_efficiency', .false.)
                    if (calc_euler_errors)   call params_set('calc_euler_errors', .false.)
                    call params_set('stockshare_fixed', .false.)
                elseif (rc ==3) then
                    call params_set('ccv', .false.)
                    write(runchar,'(a6)') ',noCCV'
                elseif (rc ==4) then
                    call params_set('scale_IR', -1.0_dp)
                    write(runchar,'(a5)') ',noIR'
                elseif (rc ==5) then
                    call params_set('scale_AR', -1.0_dp)
                    if (scale_IR_orig .ne. -1.0) call params_set('scale_IR', 0.0_dp)
                    write(runchar,'(a5)') ',noAR'
                elseif (rc ==6) then
                    call params_set('scale_AR', -1.0_dp)
                    call params_set('scale_IR', -1.0_dp)
                    write(runchar,'(a7)') ',norisk'
                else
                    print*, '-main: WARNING: rc=',rc,', but scale_IR =', scale_IR, ', scale_AR=', scale_AR
                    exit
                endif
                calib_name = calib_name_base//trim(runchar)
            endif

            print*, '- main: Setting experiment-specific parameters and saving params.txt'
            call params_set_thisrun
	        call CheckParams
	        sys_error = system('mkdir model_output/'//cal_id(calib_name))

            if (alt_insurance_calc) then
    	        call run_model(projectname, calib_name, welfare(rc,1), welfare_ins_o=welfare_ins(rc,1,:), agg_cons_o=agg_cons(rc,1))
            else
                call run_model(projectname, calib_name, welfare(rc,1), agg_cons_o=agg_cons(rc,1))
            endif

    	    if (tau_experiment) then
    	        call params_set('partial_equilibrium', .true.)
    	        call params_set('tau', tau+ tau_increment) ! better make function increase tau and put tau_increment in params
    	        write(runchar,'(a4,f4.3)') ',tau',tau
    	        calib_name = calib_name//runchar
    	        sys_error = system('mkdir model_output/'//cal_id(calib_name))
    	        if (alt_insurance_calc) then
    	            call run_model(projectname, calib_name, welfare(rc,2) ,welfare_ins_o=welfare_ins(rc,2,:), agg_cons_o=agg_cons(rc,2))
    	            cev_ins(rc,:) = welfare_ins(rc,2,:)/welfare_ins(rc,1,:) - 1.0
                else
    	            call run_model(projectname, calib_name, welfare(rc,2), agg_cons_o=agg_cons(rc,2))
                endif
                call params_set('tau', tau- tau_increment)
    	        cev(rc) = welfare(rc,2)/welfare(rc,1) - 1.0
    	        agg_cons_ratio(rc) = agg_cons(rc,1)/agg_cons(rc,2)
	        endif
	    enddo
	    if (tau_experiment .and. lbound(agg_cons,1)<1) agg_cons_ratio(lbound(agg_cons_ratio)) = agg_cons(1,1)/agg_cons(0,1)
	    if (welfare_decomposition) then
	        if (alt_insurance_calc) then
	            call write2file(welfare,cev,agg_cons_ratio,'welfare',cev_ins)
            else
                call write2file(welfare,cev,agg_cons_ratio,'welfare')
            endif
	    else
            if(size(welfare,1)>1 .or. tau_experiment) then
                if (alt_insurance_calc) then
                    call write2file(welfare,cev,agg_cons_ratio, 'welfare',cev_ins,risk_scale)
                else
                    call write2file(welfare,cev,agg_cons_ratio, 'welfare',scaling_o=risk_scale)
                endif
            endif
            if(size(welfare,1)>1) call plot('cev_regression')
        endif
    enddo

    call system_clock(end_time,count_rate,count_max) ! since the arguments are of kind integer(8), this relies on the intel fortran extension
    if (end_time < start_time) then
        elapsed_time = real(abs(end_time),dp) + real(count_max,dp) - real(start_time,dp)
    else
        elapsed_time= real(end_time - start_time,dp)
    endif
    secs= elapsed_time/real(count_rate,dp)
    print*, '*********** ..... Program ', projectname, ' completed ..... ***********'
	print '(a10, f0.2, a9, f0.2, a6)', 'CPU time: ', secs,' secs (= ', secs/60,' mins)'

contains
!-------------------------------------------------------------------------------
! Internal procedures in order:
! - subroutine get_projectname(projectname)
! - subroutine get_calibration_name(calib_name_out, exit_main_loop)
! - subroutine write2file(vector)
! - subroutine plot(mfile)
!-------------------------------------------------------------------------------
    subroutine get_projectname(projectname)
    ! Get the name of the project (parent directory)
        character(:), allocatable, intent(out) :: projectname
        character(len=100) :: path
        integer            :: i

        call get_command_argument(0,path) !returns the full path of the executable, which resides in the subdirectory of the current build config (e.g. Debug)
        i = index(trim(path), '/', .true.)
        sys_error=system('cd '//path(:i)//' && cp subdir.mk ../model_output') ! since the makefiles are in the subdir of the current build config
        projectname=trim(path(i+1:)) ! the last part is the executable name
    end subroutine get_projectname

!-------------------------------------------------------------------------------

    subroutine get_calibration_name(calib_name_out, exit_main_loop)
        character(:), allocatable, intent(out) :: calib_name_out
        logical, intent(out) :: exit_main_loop
        character(len=100) :: calib_name
        integer, save :: this_calibration_line = 1
        integer :: line, error_stat

        exit_main_loop = .false.

        open(unit=21, file='model_input/select_calibration_here.txt', status='OLD', action='READ', iostat=error_stat)

openif: if (error_stat == 0) then

	        ! fast forward to line before the current line
	        do line = 1, this_calibration_line - 1
	            read(21,'(a)')
	        enddo

	        ! read next calibration
stupid:     do ! this stupid do-loop is only here to allow for comments (preceded by an exclamation mark)
	            read(21,'(a)',iostat=error_stat) calib_name
	            if (error_stat==0) then
	               this_calibration_line = this_calibration_line + 1
		           if (index(calib_name,'!!!')>0) then
		               exit_main_loop = .true.
		           elseif (index(calib_name,'!')>0 .or. len_trim(calib_name)== 0) then
		              cycle
		           elseif (index(cal_id(calib_name),'cal')>0 .or. index(cal_id(calib_name),'GE0')>0 .or. index(cal_id(calib_name),'GE1')>0) then
		              print*, '- main: ERROR: calibration name must not contain cal, GE0, or GE1!'
		              print*, '        calib_name: '//calib_name
		              stop 'STOP'
	               endif
		        else
		           print '(a,i6)', '- main: ERROR reading calibration name in line ', line
		           exit_main_loop = .true.
	            endif

	            exit stupid

            enddo stupid

	    else openif
	        print '(a,i6)', '- main: ERROR opening select_calibration_here.txt ,IOSTAT=',error_stat
	        exit_main_loop = .true.
	    endif openif

	    close(unit=21)

	    if (error_stat > 0) then
		    open (unit=21, file='model_output/CRITICAL_ERROR.txt', status = 'replace')
	        write(21,'(a)') 'A critical error occurred in main while reading a calibration file.'
	        write(21,'(a)') 'Delete this file after reading to avoid confusion.'
	        close(21)
        endif

        calib_name_out = trim(adjustl(calib_name))

    end subroutine get_calibration_name
!-------------------------------------------------------------------------------

    subroutine write2file(welfare, cev, agg_cons_ratio, filename, cev_ins_o, scaling_o)
        real(dp) ,intent(in) :: welfare(0:,1:), cev(0:), agg_cons_ratio(0:)
        real(dp) ,intent(in) ,optional :: scaling_o(0:), cev_ins_o(0:,1:)
        character(len=*) ,intent(in) :: filename
        character(:) ,allocatable :: cal_id_temp
        real(dp) :: GE, PE, NR, IR, AR, LCI, CCV, SR, GE_INS, PE_INS, NR_INS, IR_INS, AR_INS, LCI_INS, CCV_INS, SR_INS, GE_MEAN, PE_MEAN, NR_MEAN, IR_MEAN, AR_MEAN, LCI_MEAN, CCV_MEAN, SR_MEAN
        integer  :: cev_ins_index, nr_cev

        cal_id_temp = cal_id(calib_name_base)    ! Could remove this line and put cal_id(calib_name) directly into open statement, but compiler bug.

        open(21,file='model_output/'//filename//'_'//cal_id_temp//'.txt')

        if (.not. present(scaling_o)) then

            GE = (welfare(0,1)/welfare(1,1) -1.0)*100.0
            PE =  cev(1)*100.0
            NR =  cev(6)*100.0
            IR = (cev(5)-cev(6))*100.0
            AR = (cev(4)-cev(6))*100.0
            LCI= (cev(3)*100.0-(cev(6)*100.0 + IR + AR))
            CCV= (cev(2) - cev(3))*100.0
            SR = (cev(1) - cev(2))*100.0

            GE_INS  = (welfare(0,1)/welfare(1,1) * agg_cons_ratio(0) -1.0)*100.0
            PE_INS  = ((cev(1)+1.0)*agg_cons_ratio(1)-1.0)*100.0
            NR_INS  = ((cev(6)+1.0)*agg_cons_ratio(6)-1.0)*100.0
            IR_INS  = ((cev(5)+1.0)*agg_cons_ratio(5)-1.0)*100.0
            AR_INS  = ((cev(4)+1.0)*agg_cons_ratio(4)-1.0)*100.0
            LCI_INS = ((cev(3)+1.0)*agg_cons_ratio(3)-1.0)*100.0 - IR_INS - AR_INS
            CCV_INS = ((cev(2)+1.0)*agg_cons_ratio(2) - (cev(3)+1.0)*agg_cons_ratio(3))*100.0
            SR_INS  = ((cev(1)+1.0)*agg_cons_ratio(1) - (cev(2)+1.0)*agg_cons_ratio(2))*100.0

            GE_MEAN  = GE - GE_INS
            PE_MEAN  = PE - PE_INS
            NR_MEAN  = NR - NR_INS
            IR_MEAN  = IR - IR_INS
            AR_MEAN  = AR - AR_INS
            LCI_MEAN  = LCI - LCI_INS
            CCV_MEAN  = CCV - CCV_INS
            SR_MEAN  = SR - SR_INS

            write(21,'(a)') 'Welfare changes, reported in % of consumption equivalent variation, CEV'
            write(21,'(a)') 'First line is total effect, second line insurance effect g_c^{ins}, third line mean effect g_c^{mean}'
            write(21,*)
            write(21,'(a)') '     GE     PE  CrowdOut '
            write(21,'(3(f7.2))') GE, PE, GE-PE
            write(21,'(3(f7.2))') GE_INS, PE_INS, GE_INS-PE_INS
            write(21,'(3(f7.2))') GE_MEAN, PE_MEAN, GE_MEAN-PE_MEAN
            write(21,*)
            if (surv_rates .or. debugging) then
                nr_cev = ubound(cev,1)
                write(21,'(a)') ' g_c(0,0)   g_c(0,IR)   g_c(AR,0)   g_c(AR,IR)   g_c(CCV)     g_c(SR)'
                write(21,'(<nr_cev>(3x,f6.2,3x))') (cev(nr_cev:1:-1))*100.0
                write(21,*)
                write(21,'(a)') ' g_c(0,0)    dg_c(IR)    dg_c(AR)   dg_c(LCI)   dg_c(CCV)    dg_c(SR) | dg_c(LCI)/dg_c(AR)   (dg_c(LCI)+dg_c(CCV))/PE'
                write(21,'(<nr_cev>(3x,f6.2,3x),2(10x,f6.2))') NR, IR, AR, LCI, CCV, SR, LCI/AR, (LCI + CCV)/PE
                write(21,'(<nr_cev>(3x,f6.2,3x),2(10x,f6.2))') NR_INS, IR_INS, AR_INS, LCI_INS, CCV_INS, SR_INS, LCI_INS/AR_INS, (LCI_INS + CCV_INS)/PE_INS
                write(21,'(<nr_cev>(3x,f6.2,3x),2(10x,f6.2))') NR_MEAN, IR_MEAN, AR_MEAN, LCI_MEAN, CCV_MEAN, SR_MEAN, LCI_MEAN/AR_MEAN, (LCI_MEAN + CCV_MEAN)/PE_MEAN

            else
                nr_cev = ubound(cev,1)-1
                write(21,'(a)') ' g_c(0,0)   g_c(0,IR)   g_c(AR,0)   g_c(AR,IR)   g_c(CCV)'
                write(21,'(<nr_cev>(3x,f6.2,3x))') (cev(nr_cev+1:2:-1))*100.0
                write(21,*)
                write(21,'(a)') ' g_c(0,0)    dg_c(IR)    dg_c(AR)   dg_c(LCI)   dg_c(CCV) | dg_c(LCI)/dg_c(AR)   (dg_c(LCI)+dg_c(CCV))/PE'
                write(21,'(<nr_cev>(3x,f6.2,3x),2(10x,f6.2))') NR, IR, AR, LCI, CCV, LCI/AR, (LCI + CCV)/PE
                write(21,'(<nr_cev>(3x,f6.2,3x),2(10x,f6.2))') NR_INS, IR_INS, AR_INS, LCI_INS, CCV_INS, LCI_INS/AR_INS, (LCI_INS + CCV_INS)/PE_INS
                write(21,'(<nr_cev>(3x,f6.2,3x),2(10x,f6.2))') NR_MEAN, IR_MEAN, AR_MEAN, LCI_MEAN, CCV_MEAN, LCI_MEAN/AR_MEAN, (LCI_MEAN + CCV_MEAN)/PE_MEAN
            endif

            cev_ins_index = 1 ! second run is GE without socsec, lbound is zero
        else
            if (size(welfare,1) > 1) then
                if (scale_IR_orig .ne. 0.0 .and. scale_IR_orig .ne. -1.0) then
                    write(21,'(a56)') ' IR      welfare,tau_0    welfare,tau_1              cev'
                elseif (scale_AR_orig .ne. 0.0 .and. scale_AR_orig .ne. -1.0) then
                    write(21,'(a56)') ' AR      welfare,tau_0    welfare,tau_1              cev'
                else
                    write(21,'(a)') ' WARNING: Experiments correctly set?'
                    write(21,'(a4,i2,a16,f5.2,a11,f5.2)') ' rc=',rc,', but scale_IR =', scale_IR, ', scale_AR=', scale_AR
                    write(21,'(a56)') ' risk?   welfare,tau_0    welfare,tau_1              cev'
                endif
            else
                if (scale_AR_orig == -1.0 .and. scale_IR_orig == -1.0) then
                    write(21,'(a56)') 'norisk   welfare,tau_0    welfare,tau_1              cev'
                else
                    write(21,'(a56)') ' risk    welfare,tau_0    welfare,tau_1              cev'
                endif
            endif
            do i=0,size(welfare,1)-1
                if (tau_experiment) then
                    write(21,'(f5.2,2x,3(es15.6,2x))') scaling_o(i), welfare(i,1), welfare(i,2), cev(i)
                else
                    write(21,'(f5.2,2x,es15.6,2x)') scaling_o(i), welfare(i,1)
                endif
            enddo
            cev_ins_index = 0 ! only one experiment
        endif


        if (present(cev_ins_o)) then
            write(21,*)
            write(21,*)
            write(21,*) 'New insurance calc, where the risk is removed by averaging the respective policy function.'
            write(21,*) 'Reported in % CEV, not mean adjusted!'
            write(21,*)
            write(21,'(a)') ' g_c(0,0)   g_c(0,IR)   g_c(AR,0)   g_c(AR,IR)   g_c(CCV)'
            write(21,'(<size(cev_ins_o,2)>(3x,f6.2,3x))') (cev_ins_o(cev_ins_index,5:1:-1))*100.0
            write(21,*)
        endif

        close(21)

    end subroutine write2file
!-------------------------------------------------------------------------------

    subroutine plot(mfile)
	    ! Call Matlab and plot
	    use ifport            ,only: system  ! Intel Fortran portability library
	    character(len=*), intent(in)   :: mfile
	    integer                        :: err_matl
	    character(:) ,allocatable :: cal_id_temp

	    cal_id_temp = cal_id(calib_name_base)    ! Could remove this line and put cal_id(calib_name) directly into open statement, but compiler bug.

	    err_matl = system('CALIBNAME='//cal_id_temp//' && export CALIBNAME && cd src_matlab && matlab -nodesktop -nosplash -r '//mfile//' > cl_output.txt')
	    if (err_matl ==-1) print*, 'Warning in main:plot: An error occured while plotting graphs'
    end subroutine plot

end program EPSS
