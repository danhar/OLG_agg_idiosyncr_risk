module simvars_class
    ! This derived type collects the aggregate variables over the simulations.
    use kinds      ,only: dp
    implicit none
    private

    public tSimvars, read_unformatted, write_unformatted, print_error

    type tSimvars
        integer , dimension(:), allocatable :: z     ! realizations of aggregate shock
        real(dp), dimension(:), allocatable ::    &
                K, mu, output, stock, bonds, B, invest, C, Phi_1, Phi_nx, err_aggr, err_income, eul_err_max, eul_err_avg, &      ! mu, per capita: k, bonds, consumption
                r, rf, r_pf_median, r_pf_kappa_med, wage, pens, tau, welf, bequests ! prices
        logical, dimension(:), allocatable :: err_K, err_mu

    contains
        procedure :: allocate   => allocate_simvars
        procedure :: deallocate => deallocate_simvars
        !procedure :: read_unformatted ! cannot be typ-bound procedure because need to operate on array of type tSimvars
        !procedure :: write_unformatted
        procedure, private :: get_real
        procedure :: get_logical
        generic   :: get => get_real !, get_logical !compiler bug
        procedure :: cons_grow
        procedure :: zeta
        procedure :: delta
        procedure :: K_Y
        procedure :: I_Y
        procedure :: ex_ret
        procedure :: mpk
        procedure :: bequest_rate
    end type tSimvars

    ! The following procedures cannot be type-bound procedure because they need to operate on an array of type tSimvars.
    interface read_unformatted
        module procedure read_unformatted_array
    end interface read_unformatted

    interface write_unformatted
        module procedure write_unformatted_array
    end interface write_unformatted

    interface print_error
        module procedure print_error_msg
    end interface print_error

contains
    elemental subroutine allocate_simvars(this,t)
        class(tSimvars), intent(out)  :: this
        integer,    intent(in)      :: t
        call deallocate_simvars(this)
        allocate(this%z(t), this%Phi_1(t), this%Phi_nx(t), this%err_aggr(t), this%err_income(t), this%eul_err_max(t), this%eul_err_avg(t))
        allocate(this%K(t+1),this%mu(t), this%output(t), this%stock(t), this%bonds(t), this%invest(t), this%C(t), this%welf(t)) ! Recall that stock and bond in today's per capita terms, that is why only t, not t+1
        allocate(this%r(t),this%rf(t+1), this%r_pf_median(t), this%r_pf_kappa_med(t), this%wage(t), this%pens(t), this%tau(t), this%bequests(t))
        allocate(this%B(t), this%err_K(t), this%err_mu(t))

        this%err_K = .false.
        this%err_mu = .false.

    end subroutine allocate_simvars

    elemental subroutine deallocate_simvars(this)
        class(tSimvars), intent(inout)  :: this
        ! deallocating in reverse order to allocation for memory purposes
        if (allocated(this%err_mu)) deallocate(this%err_mu)
        if (allocated(this%err_K)) deallocate(this%err_K)
        if (allocated(this%B)) deallocate(this%B)
        if (allocated(this%bequests)) deallocate(this%bequests)
        if (allocated(this%tau)) deallocate(this%tau)
        if (allocated(this%pens)) deallocate(this%pens)
        if (allocated(this%wage)) deallocate(this%wage)
        if (allocated(this%r_pf_kappa_med)) deallocate(this%r_pf_kappa_med)
        if (allocated(this%r_pf_median)) deallocate(this%r_pf_median)
        if (allocated(this%rf)) deallocate(this%rf)
        if (allocated(this%r)) deallocate(this%r)
        if (allocated(this%welf)) deallocate(this%welf)
        if (allocated(this%C)) deallocate(this%C)
        if (allocated(this%invest)) deallocate(this%invest)
        if (allocated(this%bonds)) deallocate(this%bonds)
        if (allocated(this%stock)) deallocate(this%stock)
        if (allocated(this%output)) deallocate(this%output)
        if (allocated(this%mu)) deallocate(this%mu)
        if (allocated(this%K)) deallocate(this%K)
        if (allocated(this%eul_err_avg)) deallocate(this%eul_err_avg)
        if (allocated(this%eul_err_max)) deallocate(this%eul_err_max)
        if (allocated(this%err_income)) deallocate(this%err_income)
        if (allocated(this%err_aggr)) deallocate(this%err_aggr)
        if (allocated(this%Phi_nx)) deallocate(this%Phi_nx)
        if (allocated(this%Phi_1)) deallocate(this%Phi_1)
        if (allocated(this%z)) deallocate(this%z)
    end subroutine deallocate_simvars

    pure function get_real(this,varname, lb_o, ub_o) result(get)
        real(dp) ,allocatable         :: get(:)
        class(tSimvars)  ,intent(in)  :: this
        character(len=*) ,intent(in)  :: varname
        integer, intent(in) ,optional :: lb_o, ub_o ! lower bound, upper bound
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%z)
        endif

        select case (varname)
        case ('K','k')
            !if (.not.present(ub_o)) ub = size(this%k) !better to have all of equal size - more comparable, and e.g. for covariances
            get = this%k(lb:ub)
        case ('mu')
            get = this%mu(lb:ub)
        case ('output')
            get = this%output(lb:ub)
        case ('stock')
            get = this%stock(lb:ub)
        case ('bonds')
            get = this%bonds(lb:ub)
        case ('B')
            get = this%B(lb:ub)
        case ('invest')
            get = this%invest(lb:ub)
        case ('cons','C')
            get = this%C(lb:ub)
        case ('Phi_1')
            get = this%Phi_1(lb:ub)
        case ('Phi_nx')
            get = this%Phi_nx(lb:ub)
        case ('err_aggr')
            ! Att: absolute value
            get = abs(this%err_aggr(lb:ub))
        case ('err_inc')
            ! Att: absolute value
            get = abs(this%err_income(lb:ub))
        case ('eul_err_max')
            get = abs(this%eul_err_max(lb:ub))
        case ('eul_err_avg')
            get = abs(this%eul_err_avg(lb:ub))
        case ('r')
            get = this%r(lb:ub)
        case ('rf')
            ! if (.not.present(ub_o)) ub = size(this%rf) !better to have all of equal size - more comparable, and e.g. for covariances
            get = this%rf(lb:ub)
        case ('rpf_med')
            get = this%r_pf_median(lb:ub)
        case ('rpf_kapm')
            get = this%r_pf_kappa_med(lb:ub)
        case ('netwage')
            get = this%wage(lb:ub)
        case ('pension')
            get = this%pens(lb:ub)
        case ('tau')
            get = this%tau(lb:ub)
        case ('welfare')
            get = this%welf(lb:ub)
        case ('bequests')
            get = this%bequests(lb:ub)
        case ('ex_ret')
            get = this%ex_ret(lb,ub)
        case ('mpk')
            get = this%mpk(lb,ub)
        case ('K_Y')
            get = this%K_Y(lb,ub)
        case ('I_Y')
            get = this%I_Y(lb,ub)
        case ('bequests,%')
            get = this%bequest_rate(lb,ub)
        case ('cons_grow')
            get = this%cons_grow(lb,ub)
        case ('zeta')
            get = this%zeta(lb,ub)
        case ('delta')
            get = this%delta(lb,ub)
        end select
    end function get_real

    pure function get_logical(this,varname, lb_o, ub_o) result(get)
        logical ,allocatable,dimension(:)         :: get
        class(tSimvars)  ,intent(in)  :: this
        character(len=*) ,intent(in)  :: varname
        integer, intent(in) ,optional :: lb_o, ub_o ! lower bound, upper bound
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%err_K)
        endif

        select case(varname)
        case ('err_K','err_k')
            get = this%err_K(lb:ub)
        case ('err_mu')
            get = this%err_mu(lb:ub)
        end select

    end function get_logical

    pure function cons_grow(this,lb_o,ub_o)
        real(dp), allocatable :: cons_grow(:), temp(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%C)
        endif

        temp = [0.0_dp, this%C(2:)/this%C(:size(this%C)-1) -1.0]
        cons_grow = temp(lb:ub)

    end function cons_grow

    pure function K_Y(this,lb_o,ub_o)
        ! Capital-output-ratio
        use params_mod ,only: alpha
        real(dp), allocatable :: K_Y(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%output)
        endif

        K_Y = this%K(lb:ub)**(1.0-alpha)
    end function K_Y

    pure function I_Y(this,lb_o,ub_o)
        ! Investment-output-ratio
        real(dp), allocatable :: I_Y(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%output)
        endif

        I_Y = this%invest(lb:ub)/this%output(lb:ub)
    end function I_Y

    pure function ex_ret(this,lb_o,ub_o)
        ! Excess returns
        real(dp), allocatable :: ex_ret(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%r)
        endif

        ex_ret = this%r(lb:ub) - this%rf(lb:ub)
    end function ex_ret

    pure function mpk(this,lb_o,ub_o)
        ! marginal product of capital
        use income, only: f_net_mpk
        real(dp), allocatable :: mpk(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%r)
        endif

        mpk = f_net_mpk(this%K(lb:ub), this%zeta(lb,ub), this%delta(lb,ub))
    end function mpk

    pure function bequest_rate(this,lb_o,ub_o)
        use params_mod ,only: alpha
        real(dp), allocatable :: bequest_rate(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%r)
        endif

        bequest_rate = this%bequests(lb:ub)/this%K(lb:ub)**alpha
    end function bequest_rate

    pure function zeta(this,lb_o,ub_o)
        use params_mod, only: zetaval => zeta
        real(dp), allocatable :: zeta(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%z)
        endif

        if (size(this%z) <= size(zetaval)+2) then  !mean shock equilibrium
            if (lb == 1) then
                zeta  = [1.0_dp ,  zetaval(this%z(lb+1:))]
            else
                zeta  = zetaval(this%z(lb:))
            endif
        else
            zeta = zetaval(this%z(lb:ub))
        endif
    end function zeta

    pure function delta(this,lb_o,ub_o)
        use params_mod, only: deltaval => delta
        real(dp), allocatable :: delta(:)
        class(tSimvars) ,intent(in) :: this
        integer, intent(in), optional :: lb_o, ub_o
        integer :: lb, ub

        if (present(lb_o)) then
            lb = lb_o
        else
            lb = 1
        endif

        if (present(ub_o)) then
            ub = ub_o
        else
            ub = size(this%z)
        endif

        if (size(this%z) <= size(deltaval)+2) then !mean shock equilibrium
            if (lb == 1) then
                delta = [sum(deltaval)/size(deltaval) , deltaval(this%z(lb+1:))]
            else
                delta = deltaval(this%z(lb:))
            endif
        else
            delta = deltaval(this%z(lb:ub))
        endif
    end function delta

    subroutine read_unformatted_array(this, input_path)
        type(tSimvars) ,allocatable ,intent(out) :: this(:)
        character(*) ,intent(in)                 :: input_path
        integer :: array_size, nt, i, io_stat

        open(55,file=input_path//'/simvars_sizes.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
        read(55) array_size, nt
        close(55)

        if (io_stat == 0) then

            allocate(this(array_size))
            call this%allocate(nt)

            open(55,file=input_path//'/simvars_ge.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
            ! Could I use standard derived type IO?
            do i=1,size(this)
               read(55) this(i)%z, &    ! integer
                        this(i)%K, this(i)%mu, this(i)%output,this(i)%stock,this(i)%bonds, this(i)%B, this(i)%invest, this(i)%C, this(i)%Phi_1, this(i)%Phi_nx, this(i)%err_aggr, this(i)%err_income, &
                        ! this(i)%eul_err_max, this(i)%eul_err_avg
                        this(i)%r, this(i)%rf, this(i)%r_pf_median, this(i)%r_pf_kappa_med, this(i)%wage, this(i)%pens, this(i)%tau, this(i)%welf, this(i)%bequests, &
                        this(i)%err_K, this(i)%err_mu   !logical
            enddo
            close(55)
        endif

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR reading simvars from unformatted file'
            print*, 'Try initializing with hard-coded guess (estimate_from_simvars=.false.)'
            stop 'STOP in in simvars_class:read_unformatted_array'
        endif

    end subroutine read_unformatted_array

    subroutine write_unformatted_array(this, input_path)
        ! One can't use standard derived type IO because of the allocatable components.
        class(tSimvars) ,intent(in) :: this(:)
        character(*)    ,intent(in) :: input_path
        integer :: i, io_stat

        open(55,file=input_path//'/simvars_sizes.unformatted',form='unformatted',access='stream',iostat=io_stat,action='write')
        write(55) size(this), size(this(1)%z)
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR in in simvars_class:write_unformatted_array'
        endif

        open(55,file=input_path//'/simvars_ge.unformatted',form='unformatted',access='stream',iostat=io_stat,action='write')
        do i=1,size(this)
            write(55) this(i)%z, &    ! integer
                      this(i)%K, this(i)%mu, this(i)%output,this(i)%stock,this(i)%bonds, this(i)%B, this(i)%invest, this(i)%C, this(i)%Phi_1, this(i)%Phi_nx, &
                      this(i)%err_aggr, this(i)%err_income, this(i)%eul_err_max, this(i)%eul_err_avg, &
                      this(i)%r, this(i)%rf, this(i)%r_pf_median, this(i)%r_pf_kappa_med, this(i)%wage, this(i)%pens, this(i)%tau, this(i)%welf, this(i)%bequests, &
                      this(i)%err_K, this(i)%err_mu   !logical
        enddo
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR in in simvars_class:write_unformatted_array'
        endif

    end subroutine write_unformatted_array
    !-------------------------------------------------------------------------------

    subroutine print_error_msg(simvars)
        use kinds      ,only: dp
        use params_mod ,only: nt, t_scrap
        type(tSimvars) ,intent(in) :: simvars(:)
        integer                    :: count_err, i, n
        real(dp)                   :: perc_err

        n = size(simvars)
        if (any([(simvars(i)%err_K ,i=1,n)])) then
            count_err = count([(simvars(i)%err_K(t_scrap+1:) ,i=1,n)])
            perc_err  = count_err/real((nt-t_scrap)*n,dp)*100_dp
            print 214, ' Warning: simulate_economy: # K  not in grid =', count_err,'  (',perc_err,'%)'
        endif

        if (any([(simvars(i)%err_mu ,i=1,n)])) then
            count_err = count([(simvars(i)%err_mu(t_scrap+1:) ,i=1,n)])
            perc_err  = real(count_err,dp)/real((nt-t_scrap)*n,dp)*100_dp
            print 214, ' Warning: simulate_economy: # mu not in grid =', count_err,'  (',perc_err,'%)'
        endif

    214 format((a,i6,a3,f5.1,a2))

    end subroutine print_error_msg

end module simvars_class
