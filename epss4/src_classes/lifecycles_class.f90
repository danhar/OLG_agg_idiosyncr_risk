module lifecycles_class
    use kinds      ,only: dp
    implicit none

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
        module procedure allocate_lifecycle
    end interface AllocateType

    interface DeallocateType
        module procedure deallocate_lifecycle
    end interface DeallocateType

    interface average
        module procedure average_lifecycle
    end interface

contains

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

end module lifecycles_class
