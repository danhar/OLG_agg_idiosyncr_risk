! Module defining derived types and their procedures
! Makes heavy use of Fortran 2003: allocatable component of derived type
module types
    use kinds      ,only: dp
    implicit none

    private allocate_simvars, deallocate_simvars

    type tSimvars
        integer , dimension(:), allocatable :: z     ! realizations of aggregate shock
        real(dp), dimension(:), allocatable ::    &
                K, output, stock, bonds, B, invest, mu, C, Phi_1, Phi_nx, err_aggr, err_income, &      ! mu, per capita: k, bonds, consumption
                r, rf, r_pf_median, r_pf_kappa_med, wage, pens, tau, welf, bequests ! prices
        logical,  dimension(:), allocatable :: err_K, err_mu
    end type tSimvars

    type tLifecycle
        real(dp), dimension(:), allocatable :: ap, kappa, cons, stock, cons_var, return, return_var   ! life-cycle profiles
!    contains
!        procedure :: set_number => set_number_lifecycle        ! generates compiler bug
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

contains
    pure subroutine allocate_simvars(simvars,t)
        type(tSimvars), intent(inout)  :: simvars
        integer,    intent(in)      :: t

        call deallocate_simvars(simvars)

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
