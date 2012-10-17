module statistics
    ! This module is a bit of an overkill...
    use kinds ,only: dp

    implicit none
    private
    public tStats, Correlation, covariance, lorenz_calc, lorenz_err

    type tStats
        private
        real(dp) :: avg, std, auto, absmin, absmax, cv, & ! these include all realizations with an error
                    avg_exerr, min_exerr, max_exerr       ! these exclude all realizations with an error
        integer  :: namelength=10, digits2display_short=2, digits2display_long=8
        ! character(:), allocatable: name  ! instantiate with a user-written constructor. but overkill.
    contains
        procedure :: calc_stats => calculate_statistics
        procedure :: write => write_stats
        procedure :: set_number => set_number_stat
        procedure :: avg_exerr_ ! should make this without final subscript, and the variables with
        procedure :: avg_
        procedure :: max_exerr_
        procedure :: std_
        procedure :: get_namelength
        procedure :: get_digits2display
        procedure :: writing_format
    end type tStats

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure subroutine calculate_statistics(this, series, err_mu, err_K)
! - subroutine write_to_file(this, unit, format)
! - pure subroutine set_number(this, number)
! - pure real(dp) function avg_exerr_
! - pure function avg_(this) result(result)
! - pure real(dp) function max_exerr_
! - pure integer function get_namelength(this)
! - pure integer function get_digits2display(this, format)
! - pure integer function get_format(this, format)
! - pure real(dp)function Correlation(series1,stats1,series2,stats2)
!-------------------------------------------------------------------------------

    pure subroutine calculate_statistics(this, series, err_mu, err_K)
    ! This is a specific type bound procedure
        use params_mod ,only: t_scrap, nz, stat_dist_z
        class(tStats)          ,intent(out) :: this
        real(dp) ,dimension(:) ,intent(in)  :: series
        logical  ,dimension(:) ,intent(in)  :: err_mu, err_K
        integer :: lb, n ! lb = lower bound

        if (size(series)<= nz +2) then
            lb = 2 ! mean shock
        else
            lb = t_scrap+1
        endif

        n = size(series) -(lb-1)

        if (lb ==2) then ! mean shock
        ! A bit sloppy: _avg are mean shock values, _std deviations from this mean shock
        ! might be different for arithmetic average or weighted average with stat_dist_z
            this%avg = series(1)
            this%std = sqrt(sum(stat_dist_z*(series(lb:)- this%avg)**2))
            this%absmax = abs(series(1))
            this%absmin = abs(series(1))
            this%avg_exerr = this%avg
            this%max_exerr = this%absmax
            this%max_exerr = this%absmin
        else
            this%avg = sum(series(lb:))/n
            this%std = sqrt(sum((series(lb:) - this%avg)**2)/real(n-1,dp)) ! should I rather divide by n?
            this%absmax = maxval(abs(series(lb:)))
            this%absmin = minval(abs(series(lb:)))
            this%avg_exerr = sum(pack(series(lb:),err_mu(lb:) == .false. .and. err_K(lb:) == .false.))/count(err_mu(lb:) == .false. .and. err_K(lb:) == .false.)
            this%max_exerr = maxval(abs(pack(series(lb:),err_mu(lb:) == .false. .and. err_K(lb:) == .false.)))
            this%min_exerr = minval(abs(pack(series(lb:),err_mu(lb:) == .false. .and. err_K(lb:) == .false.)))
        endif
        if (this%avg == 0.0) then
            ! This happens in particular for pensions if tau = 0.0. Not assigning a NaN helps in debugging.
            this%cv = 0.0
        else
            this%cv   = this%std/this%avg
        endif

        if (this%std == 0.0) then
            this%auto =0.0
        else
            this%auto = sum((series(lb:size(series)-1)-this%avg)*(series(lb+1:)-this%avg))/(real((n-2),dp)*this%std**2)
        endif
    end subroutine calculate_statistics

    subroutine write_stats(this, unit, name_opt, format)
    ! This is a specific type bound procedure
        class(tStats) ,intent(in) :: this
        integer       ,intent(in) :: unit
        character(len=*), intent(in), optional :: name_opt, format ! name could be in type def, but then need constructor
        character(len=this%namelength) :: name
        character(:), allocatable :: fmt1
        integer :: nl, nd

        nl = this%namelength
        if (.not.present(format)) then
            nd =2
        else
            nd = this%get_digits2display(format)
        endif
        fmt1 = this%writing_format(nd,nl)

        if (present(name_opt)) then
            name = name_opt
            name = adjustl(name)
        else
            name = repeat(' ',nl)
        endif

        if (.not. present(format)) then
            write(unit,*) this%avg, this%std, this%cv, this%auto
        else
            select case (format)
            case ('short')
                write(unit,fmt1) ' '//name,this%avg,this%std,this%cv,this%auto
            case ('long')
                write(unit,fmt1) ' '//name,this%avg,this%std,this%cv,this%auto,this%absmax,this%absmin
            case ('short_max')
                write(unit,fmt1) ' '//name,this%absmax,this%max_exerr,this%avg,this%avg_exerr
            case ('long_max')
                write(unit,fmt1) ' '//name,this%absmax,this%max_exerr,this%avg,this%avg_exerr,this%absmin,this%min_exerr
            case default
                write(unit,fmt1) ' '//name,this%avg,this%std,this%cv,this%auto
            end select
        endif
    end subroutine write_stats

    pure subroutine set_number_stat(this, number)
        class(tStats) ,intent(out) :: this
        real(dp)      ,intent(in) :: number
        this%avg       = number
        this%std       = number
        this%auto      = number
        this%absmin    = number
        this%absmax    = number
        this%cv        = number
        this%avg_exerr = number
        this%min_exerr = number
        this%max_exerr = number
    end subroutine set_number_stat

    pure function avg_exerr_(this) result(result)
        class(tStats) ,intent(in)    :: this
        real(dp) :: result
            result = this%avg_exerr
    end function avg_exerr_

    pure function avg_(this) result(result)
        class(tStats) ,intent(in)    :: this
        real(dp) :: result
            result = this%avg
    end function avg_

    pure function max_exerr_(this) result(result)
        class(tStats) ,intent(in)    :: this
        real(dp) :: result
            result = this%max_exerr
    end function max_exerr_

    pure function std_(this) result(result)
        class(tStats) ,intent(in)    :: this
        real(dp) :: result
            result = this%std
    end function std_

    pure integer function get_namelength(this)
        class(tStats) ,intent(in) :: this
        get_namelength=this%namelength
    end function get_namelength

    pure integer function get_digits2display(this, format)
        class(tStats) ,intent(in) :: this
        character(len=*), intent(in) :: format
        select case (format)
        case ('short', 'short_max')
            get_digits2display=this%digits2display_short
        case ('long','long_max')
            get_digits2display=this%digits2display_long
        case default
            get_digits2display=this%digits2display_short
        end select
    end function get_digits2display

    pure function writing_format(this,nd,nl)
        character(:), allocatable :: writing_format
        class(tStats) ,intent(in) :: this
        integer, intent(in) :: nd, nl
        character(:), allocatable :: fp
        character(2) :: nd_char, nd2_char, nl_char

        write(nd_char,'(i2)') nd
        write(nd2_char,'(i2)') nd+7
        write(nl_char,'(i2)') nl+1
        fp ='es'//trim(adjustl(nd2_char))//'.'//trim(adjustl(nd_char))
        writing_format = '(a'//trim(nl_char)//','//fp//',5(3x,'//fp//'))'
    end function writing_format


    pure real(dp) function Correlation(series1,stats1,series2,stats2)
        use params_mod, only: t_scrap, nt
        real(dp), dimension(:), intent(in) :: series1, series2
        class(tStats), intent(in)          :: stats1, stats2
        integer :: lb, n ! lb = lower bound

        if (size(series1) .ne. size(series2)) then
            Correlation = -2.0
            return
        endif

        if (size(series1)< nt) then
            lb = 2 ! mean shock
        else
            lb = t_scrap+1
        endif
        n = size(series1) -(lb-1)

        Correlation = sum((series1(lb:) - stats1%avg)*(series2(lb:) - stats2%avg))/(real(n-1,dp)*stats1%std*stats2%std)

    end function Correlation

    pure real(dp) function covariance(series1,series2,err_mu,err_K)
        use params_mod, only: t_scrap, nt
        real(dp), dimension(:), intent(in) :: series1, series2
        logical, dimension(:), optional, intent(in) :: err_mu,err_K
        real(dp), dimension(:) ,allocatable :: ser1_ex, ser2_ex
        real(dp) :: mean1, mean2
        integer :: lb, n ! lb = lower bound

        if (size(series1)< nt) then
            lb = 2 ! mean shock
        else
            lb = t_scrap+1
        endif

        if (present(err_mu) .and. present(err_K)) then
	        ser1_ex = pack(series1(lb:),err_mu(lb:) == .false. .and. err_K(lb:) == .false.)
	        ser2_ex = pack(series2(lb:),err_mu(lb:) == .false. .and. err_K(lb:) == .false.)
	        n = count(err_mu(lb:) == .false. .and. err_K(lb:) == .false.)
        else
            ser1_ex= series1(lb:)
            ser2_ex= series2(lb:)
            n = size(series1) -(lb-1)
        endif
        mean1 = sum(ser1_ex)/n
        mean2 = sum(ser2_ex)/n

        covariance = sum((ser1_ex - mean1)*(ser2_ex - mean2))/(n-1)

    end function covariance

	subroutine lorenz_calc(f,x,fx,gini,error)
	! Compute Lorenz curve and Gini coefficient by sorting vector x and using density f
	! Note: to link dlasrt2 add mkl_scalapack_core.lib to project's additional dependencies

!	use mkl_scalapack
	    implicit none
	    real(dp), dimension(:), intent(in)    :: f
	    real(dp), dimension(:), intent(inout) :: x
	    real(dp), dimension(:), intent(out)   :: fx
	    integer               , intent(out)   :: error
	    real(dp), intent(out) :: gini
	    integer, dimension(size(x)) :: key
	    integer :: n,i

	    n=size(x)
	    if (size(f)/=n .or. size(fx)/=n) then
	       error = 1
	       fx = -1.0
	       gini = -1.0
	       return
        endif

	    key=[(i,i=1,n)]
!	    call dlasrt2('I',n,x,key,error)
	    fx=f(key)
	    x=x*fx
	    gini=x(1)*fx(1)
	    do i=2,n
	        x(i)=x(i)+x(i-1)
	        gini=gini+(x(i)+x(i-1))*fx(i)
	        fx(i)=fx(i)+fx(i-1)
	    end do
	    gini=1-gini/x(n)
	    x=x/x(n)
	end subroutine lorenz_calc

    subroutine lorenz_err(error)
        integer, intent(in) :: error
        if (error == 1) then
            print *, 'ERROR in statistics:lorenz_calc: x, f, and fx must be of the same size'
            print *, ' gini and lorenz were set to -1'
        elseif (error < 1) then
            print *, 'ERROR in statistics:lorenz_calc:dlasrt2: the value of x was illegal at index ', error
        elseif (error == 0) then
            print *, 'statistics:lorenz_calc: all went fine'
        else
            print *, 'statistics:lorenz_err: Error code illegal, must be integer <=1, error=', error
        endif
    end subroutine lorenz_err

end module statistics
