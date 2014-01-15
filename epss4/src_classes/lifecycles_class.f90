module lifecycles_class
    ! This derived type collects the various lifecycle profiles.
    use kinds      ,only: dp
    implicit none
    private

    public tLifecycle, average

    type tLifecycle
        real(dp), dimension(:),     allocatable :: ap, kappa, cons, stock, cons_var, return, return_var, log_cons, var_log_cons
        real(dp), dimension(:,:,:), allocatable :: exp_value, xgrid
    contains
        procedure :: allocate => allocate_lifecycle
        procedure :: deallocate => deallocate_lifecycle
        procedure :: set_number => set_number_lifecycle
!        procedure :: average ! doesnt work because not scalar input
    end type tLifecycle

    interface average
        module procedure average_lifecycle
    end interface

contains
    elemental subroutine allocate_lifecycle(this, nj_o, nx_o)
        use params_mod, only: n_eta, nj_param => nj, nx_param => nx, nx_factor
        class(tLifecycle) ,intent(inout) :: this
        integer ,intent(in) ,optional  :: nj_o, nx_o
        integer :: nj, nx

        if (present(nj_o)) then
            nj = nj_o
        else
            nj = nj_param
        endif

        if (present(nx_o)) then
            nx = nx_o
        else
            nx = nx_param* nx_factor ! multiplying with nx_factor is specific to the use in simulation_mod
        endif

        call this%deallocate()

        allocate(this%ap(nj), this%kappa(nj), this%cons(nj), this%stock(nj), this%cons_var(nj), this%return(nj), this%return_var(nj), this%log_cons(nj), this%var_log_cons(nj))
        allocate(this%exp_value(nx,n_eta,nj), this%xgrid(nx,n_eta,nj))
    end subroutine allocate_lifecycle

    elemental subroutine deallocate_lifecycle(this)
        class(tLifecycle), intent(inout)  :: this
        if (allocated(this%ap)) deallocate(this%ap)
        if (allocated(this%kappa)) deallocate(this%kappa)
        if (allocated(this%cons)) deallocate(this%cons)
        if (allocated(this%stock)) deallocate(this%stock)
        if (allocated(this%cons_var)) deallocate(this%cons_var)
        if (allocated(this%return)) deallocate(this%return)
        if (allocated(this%return_var)) deallocate(this%return_var)
        if (allocated(this%log_cons)) deallocate(this%log_cons)
        if (allocated(this%var_log_cons)) deallocate(this%var_log_cons)
        if (allocated(this%exp_value)) deallocate(this%exp_value)
        if (allocated(this%xgrid)) deallocate(this%xgrid)
    end subroutine deallocate_lifecycle

    elemental subroutine set_number_lifecycle(this, number)
        class(tLifecycle) ,intent(inout) :: this
        real(dp)          ,intent(in) :: number
        if (.not. allocated(this%ap)) call this%allocate
        this%ap           = number
        this%kappa        = number
        this%cons         = number
        this%stock        = number
        this%cons_var     = number
        this%return       = number
        this%return_var   = number
        this%log_cons     = number
        this%var_log_cons = number
        this%exp_value    = number
        this%xgrid        = number
    end subroutine set_number_lifecycle

    pure function average_lifecycle(lc) result(average)
        ! Takes an array of type(tLifecycle) of one dimension and returns the average lifecylces.
        ! This is not a type-bound procedure, because it takes an array of tLifecycles as input.
        type(tLifecycle) ,intent(in) :: lc(:)
        type(tLifecycle) :: average
        integer :: i

        call average%allocate(size(lc(1)%ap), size(lc(1)%exp_value,1))
        call average%set_number(0.0_dp)

        do i=1,size(lc)
            average%ap           = average%ap           + lc(i)%ap          /size(lc)
            average%kappa        = average%kappa        + lc(i)%kappa       /size(lc)
            average%cons         = average%cons         + lc(i)%cons        /size(lc)
            average%stock        = average%stock        + lc(i)%stock       /size(lc)
            average%cons_var     = average%cons_var     + lc(i)%cons_var    /size(lc)
            average%return       = average%return       + lc(i)%return      /size(lc)
            average%return_var   = average%return_var   + lc(i)%return_var  /size(lc)
            average%log_cons     = average%log_cons     + lc(i)%log_cons    /size(lc)
            average%var_log_cons = average%var_log_cons + lc(i)%var_log_cons/size(lc)
            average%exp_value    = average%exp_value    + lc(i)%exp_value   /size(lc)
            average%xgrid        = average%xgrid        + lc(i)%xgrid       /size(lc)
        enddo
    end function average_lifecycle

end module lifecycles_class
