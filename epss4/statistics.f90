module statistics
    ! This module is a bit of an overkill...
    use kinds ,only: dp

    implicit none
    private
    public tStats, tStats_logical, corr, cov, lorenz_calc, lorenz_err

    type tStats
        private
        real(dp) :: avg, std, auto, absmin, absmax, cv, & ! these include all realizations with an error
                    avg_exerr, min_exerr, max_exerr       ! these exclude all realizations with an error
        real(dp), allocatable :: series(:)                ! contains the whole series (over all parallel runs), excluding the t_scrap
        integer  :: namelength=10, digits2display_short=2, digits2display_long=8
        character(:), public, allocatable:: name  ! could make private and instantiate with a user-written constructor below
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

    type tStats_logical
        integer :: count
        real(dp) :: percent
        logical(dp)  ,allocatable :: series(:)
        character(:) ,allocatable :: name
    contains
        procedure :: calc_stats => calculate_statistics_logical
    end type tStats_logical

! Compiler bug in Intel Fortran 12.1.2 (but not in 13.0.0)
!    interface tStats
!        module procedure constructor
!    end interface
!    interface tStats_logical
!        module procedure constructor_logical
!    end interface

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure function constructor(varname) result (new_stats)
! - pure function constructor_logical(varname) result (new_stats)
! - pure subroutine calculate_statistics(this, simvars)
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

    pure function constructor(varname) result (new_stats)
        character(len=*) ,intent(in)  :: varname
        type(tStats) :: new_stats
        new_stats%name = varname
    end function constructor

    pure function constructor_logical(varname) result (new_stats)
        character(len=*) ,intent(in)  :: varname
        type(tStats_logical) :: new_stats
        new_stats%name = varname
    end function constructor_logical

    pure subroutine calculate_statistics(this, simvars)
        use params_mod    ,only: t_scrap, stat_dist_z
        use simvars_class ,only: tSimvars
        class(tStats)          ,intent(inout) :: this
        type(tSimvars)         ,intent(in)  :: simvars(:)
        real(dp) ,allocatable :: seriest(:), seriesp(:)
        logical  ,allocatable :: err_k(:), err_mu(:)
        integer :: i, lb, n ! lb = lower bound

        if (size(simvars(1)%get(this%name))<= size(stat_dist_z) +2) then
            lb = 2 ! mean shock equilibrium
        else
            lb = t_scrap+1
        endif

        this%series = [(simvars(i)%get(this%name,lb) ,i=1, size(simvars))]
        err_k  = [(simvars(i)%err_k(lb:)      ,i=1, size(simvars))]
        err_mu = [(simvars(i)%err_mu(lb:)     ,i=1, size(simvars))]

        n = size(this%series)

        if (lb ==2) then ! mean shock equilibrium
        ! Mean shock stats are not important. _avg are mean shock values, _std deviations from this mean shock
        ! might be different for arithmetic average or weighted average with stat_dist_z
            seriesp  = simvars(1)%get(this%name,1,1) ! This is only temporary
            this%avg = seriesp(1)
            this%std = sqrt(sum(stat_dist_z*(this%series - this%avg)**2))
            this%absmax = abs(this%avg)
            this%absmin = abs(this%avg)
            this%avg_exerr = this%avg
            this%max_exerr = this%absmax
            this%max_exerr = this%absmin
        else
            this%avg = sum(this%series)/n
            this%std = sqrt(sum((this%series - this%avg)**2)/real(n-1,dp)) ! should I rather divide by n?
            this%absmax = maxval(abs(this%series))
            this%absmin = minval(abs(this%series))
            this%avg_exerr = sum(pack(this%series,err_mu == .false. .and. err_K == .false.))/count(err_mu == .false. .and. err_K == .false.)
            this%max_exerr = maxval(abs(pack(this%series,err_mu == .false. .and. err_K == .false.)))
            this%min_exerr = minval(abs(pack(this%series,err_mu == .false. .and. err_K == .false.)))
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
            seriest = [(simvars(i)%get(this%name,lb,size(simvars(1)%get(this%name))-1) ,i=1, size(simvars))]
            seriesp = [(simvars(i)%get(this%name,lb+1) ,i=1, size(simvars))]

            this%auto = sum((seriest-this%avg)*(seriesp-this%avg))/(real((n-2),dp)*this%std**2)
        endif
    end subroutine calculate_statistics

    pure subroutine calculate_statistics_logical(this, simvars)
        use params_mod    ,only: t_scrap, stat_dist_z
        use simvars_class ,only: tSimvars
        class(tStats_logical) ,intent(inout) :: this
        type(tSimvars)        ,intent(in)    :: simvars(:)
        integer :: i, lb ! lb = lower bound

        if (size(simvars(1)%get_logical(this%name))<= size(stat_dist_z) +2) then
            lb = 2 ! mean shock equilibrium
        else
            lb = t_scrap+1
        endif

        this%series  = [(simvars(i)%get_logical(this%name,lb) ,i=1, size(simvars))]
        this%count   = count(this%series) ! counts only .true.
        this%percent = real(this%count,dp)/real(size(this%series),dp) * 100.0

    end subroutine calculate_statistics_logical


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


    pure real(dp) function corr(var1,var2)
        ! correlation, including all realizations (also those where err_mu or err_K true)
        type(tStats) ,intent(in) :: var1,var2

        if (size(var1%series) .ne. size(var2%series)) then
            corr = -2.0
            return
        endif

        corr = sum((var1%series - var1%avg)*(var2%series - var2%avg))/(real(size(var1%series)-1,dp)*var1%std*var2%std)

    end function corr

    pure real(dp) function cov(var1,var2)
        ! covariance, including all realizations (also those where err_mu or err_K true)
        type(tStats) ,intent(in) :: var1,var2

        if (size(var1%series) .ne. size(var2%series)) then
            cov = -huge(1.0_dp)
            return
        endif

        cov = sum((var1%series - var1%avg)*(var2%series - var2%avg))/real((size(var1%series)-1),dp)
    end function cov

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
