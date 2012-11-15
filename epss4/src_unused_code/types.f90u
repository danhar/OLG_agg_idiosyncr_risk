! Module defining derived types and their procedures
! Makes heavy use of Fortran 2003: allocatable component of derived type
module types
    use kinds      ,only: dp
    implicit none

    private allocate_simvars, deallocate_simvars, cons_grow, delta, zeta, K_Y, ex_ret, bequest_rate

    type tSimvars
        integer , dimension(:), allocatable :: z     ! realizations of aggregate shock
        real(dp), dimension(:), allocatable ::    &
                K, output, stock, bonds, B, invest, mu, C, Phi_1, Phi_nx, err_aggr, err_income, &      ! mu, per capita: k, bonds, consumption
                r, rf, r_pf_median, r_pf_kappa_med, wage, pens, tau, welf, bequests ! prices
        logical,  dimension(:), allocatable :: err_K, err_mu

    contains
        procedure, private :: get_real
        procedure :: get_logical
        generic   :: get => get_real !, get_logical !compiler bug
        procedure :: cons_grow
        procedure :: zeta
        procedure :: delta
        procedure :: K_Y
        procedure :: ex_ret
        procedure :: bequest_rate
    end type tSimvars

    type tLifecycle
    !?!? Why do I need them to be allocatable?
        real(dp), dimension(:), allocatable :: ap, kappa, cons, stock, cons_var, return, return_var   ! life-cycle profiles
!    contains
!        procedure :: set_number => set_number_lifecycle        ! generates compiler bug
!        procedure :: average ! doesnt work because not scalar input
    end type tLifecycle

	interface set_number
    	module procedure set_number_lifecycle
	end interface set_number

    interface AllocateType
        module procedure allocate_simvars, allocate_lifecycle
    end interface AllocateType

    interface DeallocateType
        module procedure deallocate_simvars, deallocate_lifecycle
    end interface DeallocateType

    interface average
        module procedure average_lifecycle
    end interface

contains
    elemental subroutine allocate_simvars(simvars,t)
        type(tSimvars), intent(out)  :: simvars
        integer,    intent(in)      :: t

        allocate(simvars%z(t), simvars%Phi_1(t), simvars%Phi_nx(t), simvars%err_aggr(t), simvars%err_income(t))
        allocate(simvars%K(t+1),simvars%mu(t), simvars%output(t), simvars%stock(t), simvars%bonds(t), simvars%invest(t), simvars%C(t), simvars%welf(t)) ! Recall that stock and bond in today's per capita terms, that is why only t, not t+1
        allocate(simvars%r(t),simvars%rf(t+1), simvars%r_pf_median(t), simvars%r_pf_kappa_med(t), simvars%wage(t), simvars%pens(t), simvars%tau(t), simvars%bequests(t))
        allocate(simvars%B(t), simvars%err_K(t), simvars%err_mu(t))

        simvars%err_K = .false.
        simvars%err_mu = .false.

    end subroutine allocate_simvars

    pure subroutine deallocate_simvars(simvars)
        type(tSimvars), intent(inout)  :: simvars

        if (allocated(simvars%z)) deallocate(simvars%z)
        if (allocated(simvars%K)) deallocate(simvars%K)
        if (allocated(simvars%mu)) deallocate(simvars%mu)
        if (allocated(simvars%output)) deallocate(simvars%output)
        if (allocated(simvars%stock)) deallocate(simvars%stock)
        if (allocated(simvars%bonds)) deallocate(simvars%bonds)
        if (allocated(simvars%invest)) deallocate(simvars%invest)
        if (allocated(simvars%C)) deallocate(simvars%C)
        if (allocated(simvars%Phi_1)) deallocate(simvars%Phi_1)
        if (allocated(simvars%Phi_nx)) deallocate(simvars%Phi_nx)
        if (allocated(simvars%err_aggr)) deallocate(simvars%err_aggr)
        if (allocated(simvars%B)) deallocate(simvars%B)
        if (allocated(simvars%err_income)) deallocate(simvars%err_income)
        if (allocated(simvars%r)) deallocate(simvars%r)
        if (allocated(simvars%rf)) deallocate(simvars%rf)
        if (allocated(simvars%r_pf_median)) deallocate(simvars%r_pf_median)
        if (allocated(simvars%r_pf_kappa_med)) deallocate(simvars%r_pf_kappa_med)
        if (allocated(simvars%wage)) deallocate(simvars%wage)
        if (allocated(simvars%pens)) deallocate(simvars%pens)
        if (allocated(simvars%tau)) deallocate(simvars%tau)
        if (allocated(simvars%welf)) deallocate(simvars%welf)
        if (allocated(simvars%bequests)) deallocate(simvars%bequests)
        if (allocated(simvars%err_K)) deallocate(simvars%err_K)
        if (allocated(simvars%err_mu)) deallocate(simvars%err_mu)

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
        case ('cons')
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
        case ('K_Y')
            get = this%K_Y(lb,ub)
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


    pure subroutine allocate_lifecycle(lc,nj)
        class(tLifecycle), intent(inout)  :: lc
        integer,    intent(in)      :: nj
        call deallocate_lifecycle(lc)
        allocate(lc%ap(nj), lc%kappa(nj), lc%cons(nj), lc%stock(nj), lc%cons_var(nj), lc%return(nj), lc%return_var(nj))
    end subroutine allocate_lifecycle

    pure subroutine deallocate_lifecycle(lc)
        class(tLifecycle), intent(inout)  :: lc
        if (allocated(lc%ap)) deallocate(lc%ap)
        if (allocated(lc%kappa)) deallocate(lc%kappa)
        if (allocated(lc%cons)) deallocate(lc%cons)
        if (allocated(lc%stock)) deallocate(lc%stock)
        if (allocated(lc%cons_var)) deallocate(lc%cons_var)
        if (allocated(lc%return)) deallocate(lc%return)
        if (allocated(lc%return_var)) deallocate(lc%return_var)
    end subroutine deallocate_lifecycle

    pure function average_lifecycle(lc) result(average)
        type(tLifecycle) ,intent(in) :: lc(:)
        type(tLifecycle) :: average
        integer :: i

        call allocate_lifecycle(average,size(lc(1)%ap))   !!!! TURN into type bound procedure!!
        average%ap         = 0.0
        average%kappa      = 0.0
        average%cons       = 0.0
        average%stock      = 0.0
        average%cons_var   = 0.0
        average%return     = 0.0
        average%return_var = 0.0

        do i=1,size(lc)
            average%ap         = average%ap + lc(i)%ap      /size(lc)
            average%kappa      = average%kappa + lc(i)%kappa   /size(lc)
            average%cons       = average%cons + lc(i)%cons    /size(lc)
            average%stock      = average%stock + lc(i)%stock   /size(lc)
            average%cons_var   = average%cons_var + lc(i)%cons_var/size(lc)
            average%return     = average%return + lc(i)%return  /size(lc)
            average%return_var = average%return_var + lc(i)%return_var/size(lc)
        enddo
    end function average_lifecycle

    pure subroutine set_number_lifecycle(this, number)
        ! turning this into a specific type-bound procedure generates a compiler bug
        use params_mod, only: nj
        class(tLifecycle) ,intent(out) :: this
        real(dp)          ,intent(in) :: number
        call allocate_lifecycle(this,nj)
        this%ap         = number
        this%kappa      = number
        this%cons       = number
        this%stock      = number
        this%cons_var   = number
        this%return     = number
        this%return_var = number
    end subroutine set_number_lifecycle

end module types
