! main.f90 of program Equity Premium Social Security
! Copyleft 2012 Daniel Harenberg
program EPSS

    use ifport             ,only: system  ! Intel Fortran portability library
	use params_mod         ,only: SetDefaultValues,ReadCalibration, SetRemainingParams, CheckParams, SaveParams, cal_id, params_set, params_set_thisrun, &
	                              n_end_params, run_n_times, run_counter_start, twosided_experiment, scale_AR, scale_IR, scale_AR_orig, scale_IR_orig, tau_experiment, tau, dp
	use calibration_mod    ,only: calibrate
	use run_model_mod

    implicit none
	real(dp)                  :: secs
	real(dp)     ,allocatable :: welfare(:,:), cev(:), risk_scale(:)
	integer                   :: sys_error, rc, i, start_time, end_time, count_rate
	logical                   :: exit_main_loop
	character(:) ,allocatable :: projectname, calib_name, calib_name_base
	character(len=7)          :: runchar
    real(dp), parameter       :: tau_increment =0.02_dp

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

	    print*, '- main: Reading calibration file '// calib_name
	    call ReadCalibration(trim(adjustl(calib_name)))
	    call SetRemainingParams

	    if (n_end_params > 0) then
	        call params_set_thisrun
            call CheckParams
	        call calibrate(projectname, calib_name)
        endif

        if (twosided_experiment .and. run_n_times>1) call params_set('run_counter_start', -1*run_n_times+2)
        if (allocated(welfare)) deallocate(welfare)
        allocate(welfare(run_counter_start:run_n_times,2))
        if (allocated(risk_scale)) deallocate(risk_scale)
        allocate(risk_scale(run_counter_start:run_n_times))
        if (allocated(cev)) deallocate(cev)
        if (tau_experiment)  allocate(cev(run_counter_start:run_n_times))
        write (runchar, *) ' '
        calib_name_base = calib_name
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
	        call SaveParams(projectname, calib_name)

    	    call run_model(projectname, calib_name, welfare(rc,1))

    	    if (tau_experiment) then
    	        call params_set('partial_equilibrium', .true.)
    	        call params_set('tau', tau+ tau_increment) ! better make function increase tau and put tau_increment in params
    	        calib_name = calib_name//',tau'
    	        sys_error = system('mkdir model_output/'//cal_id(calib_name))
                call SaveParams(projectname, calib_name)
    	        call run_model(projectname, calib_name, welfare(rc,2))
    	        call params_set('tau', tau- tau_increment)
    	        cev(rc) = welfare(rc,2)/welfare(rc,1) - 1.0
	        endif
	    enddo
        if(size(welfare,1)>1 .or. tau_experiment) call write2file(welfare,cev,risk_scale,'welfare')
        if(size(welfare,1)>1) call plot('cev_regression')
    enddo

    call system_clock(end_time,count_rate)
    secs= real(end_time - start_time,dp)/real(count_rate,dp)
    print*, '*********** ..... Program ', projectname, ' completed ..... ***********'
	print '(a10, f8.2, a9, f6.2, a6)', 'CPU time: ', secs,' secs (= ', secs/60,' mins)'

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

    subroutine write2file(welfare, cev, scaling, filename)
        real(dp) ,intent(in) :: welfare(:,:), cev(:), scaling(:)
        character(len=*) ,intent(in) :: filename
        character(:) ,allocatable :: cal_id_temp

        cal_id_temp = cal_id(calib_name_base)    ! Could remove this line and put cal_id(calib_name) directly into open statement, but compiler bug.

        open(21,file='model_output/'//filename//'_'//cal_id_temp//'.txt')
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
        do i=1,size(welfare,1)
	        if (tau_experiment) then
	            write(21,'(f5.2,2x,3(es15.6,2x))') scaling(i), welfare(i,1), welfare(i,2), cev(i)
	        else
	            write(21,'(f5.2,2x,es15.6,2x)') scaling(i), welfare(i,1)
	        endif
        enddo
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
