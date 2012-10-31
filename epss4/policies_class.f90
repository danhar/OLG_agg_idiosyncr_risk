module policies_class
    use kinds
    implicit none
    private

    public tPolicies

    type tPolicies
        real(dp), allocatable, dimension(:,:,:,:,:,:) :: apgrid, kappa, stocks, xgrid ! policies /grids
    contains
        procedure :: allocate => allocate_policies
        procedure :: deallocate => deallocate_policies
        procedure :: calc_kappa
        procedure :: consumption    ! At the moment, I am not using this anywhere
        procedure :: mean           ! At the moment, I am not using this anywhere
        procedure :: interpolate    ! At the moment, I am not using this anywhere
    end type tPolicies

contains
!-------------------------------------------------------------------------------
! Module procedures in order:
! - pure subroutine allocate_policies(p,nk,nmu)
! - pure subroutine deallocate_policies(p)
! - pure subroutine calc_kappa(p)
! - pure function consumption(this)
! - pure function mean(this,dimension_o,weight_o) result(mean_policy)
! - pure function interpolate(this,dim_x, gridx, x) result(pol_int)
!-------------------------------------------------------------------------------

pure subroutine allocate_policies(this,nz,nk,nmu)
    use params_mod, only: nx, n_eta, nj
    class(tPolicies), intent(inout)  :: this
    integer,    intent(in)      :: nz,nk,nmu
    call deallocate_policies(this)
    allocate(this%apgrid(nx,n_eta,nz,nj,nk,nmu),this%kappa(nx,n_eta,nz,nj,nk,nmu),this%stocks(nx,n_eta,nz,nj,nk,nmu),this%xgrid(nx,n_eta,nz,nj,nk,nmu))
end subroutine allocate_policies
!-------------------------------------------------------------------------------

pure subroutine deallocate_policies(this)
    class(tPolicies), intent(inout)  :: this
    ! deallocating in reverse order to allocation for memory purposes
    if (allocated(this%xgrid)) deallocate(this%xgrid)
    if (allocated(this%stocks)) deallocate(this%stocks)
    if (allocated(this%kappa)) deallocate(this%kappa)
    if (allocated(this%apgrid)) deallocate(this%apgrid)
end subroutine deallocate_policies
!-------------------------------------------------------------------------------

pure subroutine calc_kappa(this)
    class(tPolicies), intent(inout)  :: this
    where (this%apgrid .ne. 0.0)
        this%kappa = this%stocks/this%apgrid
    elsewhere
        this%kappa = 0.0
    end where
end subroutine calc_kappa
!-------------------------------------------------------------------------------

pure function consumption(this)
    real(dp), dimension(:,:,:,:,:,:), allocatable :: consumption
    class(tPolicies), intent(in)  :: this

    consumption = this%xgrid - this%apgrid

end function consumption
!-------------------------------------------------------------------------------

pure function mean(this,dimension_o,weight_o) result(mean_policy)
    ! At the moment, I am not using this anywhere, because typically I want the result to have one dimension less, which is not the case here
    class(tPolicies), intent(in)  :: this
    type(tPolicies)               :: mean_policy
    integer, intent(in), optional:: dimension_o
    real(dp), intent(in), optional:: weight_o
    integer :: dimension, nd, dc
    real(dp) :: weight

    if (.not. present(dimension_o)) then
        dimension = 3 ! corresponds to aggregate shocks nz
    else
        dimension = dimension_o
    endif

    nd = size(this%apgrid,dimension)

    if (.not. present(weight_o)) then
        weight = 1.0/nd
    else
        weight = weight_o
    endif

    select case (dimension)
    case (3)
        call mean_policy%allocate(1,size(this%apgrid,5),size(this%apgrid,6)) ! policies for given z, K, and mu
        mean_policy%apgrid = 0.0
        mean_policy%stocks = 0.0
        mean_policy%xgrid  = 0.0
        do dc = 1,nd
            mean_policy%apgrid(:,:,1,:,:,:) = mean_policy%apgrid(:,:,1,:,:,:) + weight*this%apgrid(:,:,dc,:,:,:)
            mean_policy%stocks(:,:,1,:,:,:) = mean_policy%stocks(:,:,1,:,:,:) + weight*this%stocks(:,:,dc,:,:,:)
            mean_policy%xgrid (:,:,1,:,:,:) = mean_policy%xgrid (:,:,1,:,:,:) + weight*this%xgrid (:,:,dc,:,:,:)
        enddo
    case default ! same as case 3
        call mean_policy%allocate(1,size(this%apgrid,5),size(this%apgrid,6)) ! policies for given z, K, and mu
        mean_policy%apgrid = 0.0
        mean_policy%stocks = 0.0
        mean_policy%xgrid  = 0.0
        do dc = 1,nd
            mean_policy%apgrid(:,:,1,:,:,:) = mean_policy%apgrid(:,:,1,:,:,:) + weight*this%apgrid(:,:,dc,:,:,:)
            mean_policy%stocks(:,:,1,:,:,:) = mean_policy%stocks(:,:,1,:,:,:) + weight*this%stocks(:,:,dc,:,:,:)
            mean_policy%xgrid (:,:,1,:,:,:) = mean_policy%xgrid (:,:,1,:,:,:)  + weight*this%xgrid(:,:,dc,:,:,:)
        enddo
    end select

    call mean_policy%calc_kappa
end function mean
!-------------------------------------------------------------------------------

pure function interpolate(this,dim_x, gridx, x) result(pol_int)
    ! linear interpolation for one value in one dimension
    ! At the moment, I am not using this anywhere, because typically I want the result to have one dimension less, which is not the case here
    ! For this, would need another type tPolicies5 (or Policies_d5) and
    ! class(tPolicies),allocatable  :: pol_int ! allocate later
    ! or
    ! type(Policies_d5)  :: pol_int
    use fun_locate      ,only: f_locate

    class(tPolicies), intent(in)  :: this
    class(tPolicies),allocatable  :: pol_int
    integer, intent(in)           :: dim_x
    real(dp), intent(in) :: x, gridx(:)
    integer :: i
    real(dp) :: w

    i        = f_locate(gridx, x)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
    w        = (x - gridx(i)) / (gridx(i+1) - gridx(i))
    ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds

    select case (dim_x)
    case (5)
        call pol_int%allocate(size(this%apgrid,3),1,size(this%apgrid,6)) ! policies for given z, K, and mu
        pol_int%apgrid(:,:,:,:,1,:)= (1-w)*this%apgrid(:,:,:,:,i,:) +w*this%apgrid(:,:,:,:,i+1,:)
        pol_int%stocks(:,:,:,:,1,:)= (1-w)*this%stocks(:,:,:,:,i,:) +w*this%stocks(:,:,:,:,i+1,:)
        pol_int%xgrid (:,:,:,:,1,:) = (1-w)*this%xgrid (:,:,:,:,i,:) +w*this%xgrid (:,:,:,:,i+1,:)
    case default ! same as case 5
        call pol_int%allocate(size(this%apgrid,3),1,size(this%apgrid,6)) ! policies for given z, K, and mu
        pol_int%apgrid(:,:,:,:,1,:)= (1-w)*this%apgrid(:,:,:,:,i,:) +w*this%apgrid(:,:,:,:,i+1,:)
        pol_int%stocks(:,:,:,:,1,:)= (1-w)*this%stocks(:,:,:,:,i,:) +w*this%stocks(:,:,:,:,i+1,:)
        pol_int%xgrid (:,:,:,:,1,:) = (1-w)*this%xgrid (:,:,:,:,i,:) +w*this%xgrid (:,:,:,:,i+1,:)
    end select

    call pol_int%calc_kappa
end function interpolate
!-------------------------------------------------------------------------------

end module policies_class
