module lifecycles_class
    ! This derived type collects the various lifecycle profiles.
    use kinds      ,only: dp
    implicit none
    private

    public tLifecycle

    type tLifecycle
        real(dp), dimension(:), allocatable :: ap, kappa, cons, stock, cons_var, return, return_var
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
    elemental subroutine allocate_lifecycle(this, nj_o)
        use params_mod, only: nj_param => nj
        class(tLifecycle) ,intent(out) :: this
        integer ,intent(in) ,optional  :: nj_o
        integer :: nj

        if (present(nj_o)) then
            nj = nj_o
        else
            nj = nj_param
        endif

        allocate(this%ap(nj), this%kappa(nj), this%cons(nj), this%stock(nj), this%cons_var(nj), this%return(nj), this%return_var(nj))
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
    end subroutine deallocate_lifecycle

    elemental subroutine set_number_lifecycle(this, number)
        class(tLifecycle) ,intent(out) :: this
        real(dp)          ,intent(in) :: number
        if (.not. allocated(this%ap)) call this%allocate
        this%ap         = number
        this%kappa      = number
        this%cons       = number
        this%stock      = number
        this%cons_var   = number
        this%return     = number
        this%return_var = number
    end subroutine set_number_lifecycle

    pure function average_lifecycle(lc) result(average)
        ! Takes an array of type(tLifecycle) of one dimension and returns the average lifecylces.
        ! This is not a type-bound procedure, because it takes an array of tLifecycles as input.
        type(tLifecycle) ,intent(in) :: lc(:)
        type(tLifecycle) :: average
        integer :: i

        call average%allocate(size(lc(1)%ap))

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

end module lifecycles_class
