!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

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
        procedure :: mean
        procedure :: interpolate
        procedure :: read_unformatted_kappa
        procedure :: write_unformatted_kappa
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

pure subroutine allocate_policies(this,nx,nz,nk,nmu)
    use params_mod, only: n_eta, nj
    class(tPolicies), intent(inout)  :: this
    integer,    intent(in)      :: nx,nz,nk,nmu
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

subroutine read_unformatted_kappa(this, equilibrium_type, input_path)
    ! After having allocated all policies, read portfolio share kappa only!
    use params_mod, only: params_set
    class(tPolicies) ,intent(inout) :: this
    character(*)   ,intent(in)  :: equilibrium_type, input_path
    integer :: io_stat

    open(55,file=input_path//'/kappa_'//equilibrium_type//'.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
    if (io_stat .ne. 0) then
        print*, 'I/O ERROR opening unformatted file in policies_class:read_unformatted_kappa.'
        print*, 'Check path: '//input_path
        !stop 'STOP in in open policies_class:read_unformatted'
    endif

    read(55,iostat=io_stat) this%kappa
    if (io_stat .ne. 0) then
        print*, 'I/O ERROR reading kappa from unformatted file in policies_class:read_unformatted_kappa.'
        print*, 'Check path: '//input_path
        !stop 'STOP in in read policiess_class:read_unformatted'
    endif

    close(55)

end subroutine read_unformatted_kappa
!-------------------------------------------------------------------------------

subroutine write_unformatted_kappa(this, equilibrium_type, input_path)
    ! Write portfolio share kappa only!
    class(tPolicies) ,intent(in) :: this
    character(*)   ,intent(in) :: equilibrium_type, input_path
    integer :: io_stat

    open(55,file=input_path//'/kappa_'//equilibrium_type//'.unformatted',form='unformatted',access='stream',iostat=io_stat,action='write')
    if (io_stat .ne. 0) then
        print*, 'I/O ERROR opening unformatted file in policies_class:write_unformatted_kappa.'
        print*, 'Check path: '//input_path
        !stop 'STOP in in open policies_class:read_unformatted'
    endif

    write(55,iostat=io_stat) this%kappa
    if (io_stat .ne. 0) then
        print*, 'I/O ERROR writing kappa to unformatted file in policies_class:write_unformatted_kappa.'
        print*, 'Check path: '//input_path
        !stop 'STOP in in read policiess_class:read_unformatted'
    endif

    close(55)

end subroutine write_unformatted_kappa
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

pure function mean(this,dimension,weight_o) result(mean_policy)
    ! this returns the mean along a dimension. The result has same dimensionality as the original,
    ! but the same values in the respective dimension.
    class(tPolicies), intent(in)  :: this
    type(tPolicies)               :: mean_policy
    integer, intent(in)           :: dimension
    real(dp), intent(in), optional:: weight_o(:)
    real(dp) ,dimension(:,:,:,:,:) ,allocatable :: apgrid, stocks, xgrid
    integer :: nd, dc
    real(dp) ,allocatable:: weight(:)

    nd = size(this%apgrid,dimension)
    if (.not. present(weight_o)) then
        allocate(weight(nd))
        weight = 1.0/nd
    else
        weight = weight_o
    endif

    select case (dimension)
    case (1)
        ! The following causes a compiler bug
        ! apgrid = weight(1)*this%apgrid(1,:,:,:,:,:)
        !stocks = weight(1)*this%stocks(1,:,:,:,:,:)
        !xgrid  = weight(1)*this%xgrid (1,:,:,:,:,:)
        ! So instead do in two steps
        apgrid = this%apgrid(1,:,:,:,:,:)
        apgrid = weight(1)*apgrid
        stocks = this%stocks(1,:,:,:,:,:)
        stocks = weight(1)*stocks
        xgrid  = this%xgrid (1,:,:,:,:,:)
        xgrid  = weight(1)*xgrid

        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(dc,:,:,:,:,:)
            stocks = stocks + weight(dc)*this%stocks(dc,:,:,:,:,:)
            xgrid  = xgrid  + weight(dc)*this%xgrid (dc,:,:,:,:,:)
        enddo

    case (2)
        ! The following causes a compiler bug
        ! apgrid = weight(1)*this%apgrid(:,1,:,:,:,:)
        ! So instead do in two steps
        apgrid = this%apgrid(:,1,:,:,:,:)
        apgrid = weight(1)*apgrid
        stocks = this%stocks(:,1,:,:,:,:)
        stocks = weight(1)*stocks
        xgrid  = this%xgrid (:,1,:,:,:,:)
        xgrid  = weight(1)*xgrid
        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(:,dc,:,:,:,:)
            stocks = stocks + weight(dc)*this%stocks(:,dc,:,:,:,:)
            xgrid  = xgrid  + weight(dc)*this%xgrid (:,dc,:,:,:,:)
        enddo

    case (3)
        apgrid = this%apgrid(:,:,1,:,:,:)
        apgrid = weight(1)*apgrid
        stocks = this%stocks(:,:,1,:,:,:)
        stocks = weight(1)*stocks
        xgrid  = this%xgrid (:,:,1,:,:,:)
        xgrid  = weight(1)*xgrid
        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(:,:,dc,:,:,:)
            stocks = stocks + weight(dc)*this%stocks(:,:,dc,:,:,:)
            xgrid  = xgrid  + weight(dc)*this%xgrid (:,:,dc,:,:,:)
        enddo

    case (4)
        apgrid = this%apgrid(:,:,:,1,:,:)
        apgrid = weight(1)*apgrid
        stocks = this%stocks(:,:,:,1,:,:)
        stocks = weight(1)*stocks
        xgrid  = this%xgrid (:,:,:,1,:,:)
        xgrid  = weight(1)*xgrid
        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(:,:,:,dc,:,:)
            stocks = stocks + weight(dc)*this%stocks(:,:,:,dc,:,:)
            xgrid  = xgrid  + weight(dc)*this%xgrid (:,:,:,dc,:,:)
        enddo

    case (5)
        apgrid = this%apgrid(:,:,:,:,1,:)
        apgrid = weight(1)*apgrid
        stocks = this%stocks(:,:,:,:,1,:)
        stocks = weight(1)*stocks
        xgrid  = this%xgrid (:,:,:,:,1,:)
        xgrid  = weight(1)*xgrid

        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(:,:,:,:,dc,:)
            stocks = stocks + weight(dc)*this%stocks(:,:,:,:,dc,:)
            xgrid  = xgrid  + weight(dc)*this%xgrid (:,:,:,:,dc,:)
        enddo

    case (6)
        apgrid = this%apgrid(:,:,:,:,:,1)
        apgrid = weight(1)*apgrid
        stocks = this%stocks(:,:,:,:,:,1)
        stocks = weight(1)*stocks
        xgrid  = this%xgrid (:,:,:,:,:,1)
        xgrid  = weight(1)*xgrid

        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(:,:,:,:,:,dc)
            stocks = stocks + weight(dc)*this%stocks(:,:,:,:,:,dc)
            xgrid  = xgrid  + weight(dc)*this%xgrid (:,:,:,:,:,dc)
        enddo

    case default ! same as case 3
        apgrid = this%apgrid(:,:,1,:,:,:)
        stocks = this%stocks(:,:,1,:,:,:)
        xgrid  = this%xgrid (:,:,1,:,:,:)
        apgrid = weight(1)*apgrid
        stocks = weight(1)*stocks
        xgrid  = weight(1)*xgrid

        do dc = 2,nd
            apgrid = apgrid + weight(dc)*this%apgrid(:,:,dc,:,:,:)
            stocks = stocks + weight(dc)*this%stocks(:,:,dc,:,:,:)
            xgrid  = xgrid  + weight(dc)*this%xgrid (:,:,dc,:,:,:)
        enddo
    end select

    mean_policy%apgrid = spread(apgrid, dimension, nd)
    mean_policy%stocks = spread(stocks, dimension, nd)
    mean_policy%xgrid  = spread(xgrid , dimension, nd)

    call mean_policy%calc_kappa
end function mean
!-------------------------------------------------------------------------------

pure function interpolate(this,dim_x, gridx, x) result(pol_int)
    ! linear interpolation for one value in one dimension
    ! the interpolated dimension has size 1 when returned
    ! Due to compiler bugs with allocatable arrays its a bit more complicated than necessary
    use fun_locate      ,only: f_locate

    class(tPolicies), intent(in) :: this
    type (tPolicies)             :: pol_int
    integer, intent(in)          :: dim_x
    real(dp), intent(in) :: x, gridx(:)
    real(dp), dimension(:,:,:,:,:), allocatable :: apgrid, stocks, xgrid
    integer :: i, nd
    real(dp) :: w

    nd = size(this%apgrid,dim_x) ! need only because compiler bug

    i        = f_locate(gridx, x)   ! In 'default', returns iu-1 if x>xgrid(iu-1)
    w        = (x - gridx(i)) / (gridx(i+1) - gridx(i))
    ! If w>1 or w<0 we get linear extrapolation at upper or lower bounds

    select case (dim_x)
    case (1)
        apgrid= (1-w)*this%apgrid(i,:,:,:,:,:) +w*this%apgrid(i+1,:,:,:,:,:)
        stocks= (1-w)*this%stocks(i,:,:,:,:,:) +w*this%stocks(i+1,:,:,:,:,:)
        xgrid = (1-w)*this%xgrid (i,:,:,:,:,:) +w*this%xgrid (i+1,:,:,:,:,:)
    case (2)
        apgrid= (1-w)*this%apgrid(:,i,:,:,:,:) +w*this%apgrid(:,i+1,:,:,:,:)
        stocks= (1-w)*this%stocks(:,i,:,:,:,:) +w*this%stocks(:,i+1,:,:,:,:)
        xgrid = (1-w)*this%xgrid (:,i,:,:,:,:) +w*this%xgrid (:,i+1,:,:,:,:)
    case (3)
        apgrid= (1-w)*this%apgrid(:,:,i,:,:,:) +w*this%apgrid(:,:,i+1,:,:,:)
        stocks= (1-w)*this%stocks(:,:,i,:,:,:) +w*this%stocks(:,:,i+1,:,:,:)
        xgrid = (1-w)*this%xgrid (:,:,i,:,:,:) +w*this%xgrid (:,:,i+1,:,:,:)
    case (4)
        apgrid= (1-w)*this%apgrid(:,:,:,i,:,:) +w*this%apgrid(:,:,:,i+1,:,:)
        stocks= (1-w)*this%stocks(:,:,:,i,:,:) +w*this%stocks(:,:,:,i+1,:,:)
        xgrid = (1-w)*this%xgrid (:,:,:,i,:,:) +w*this%xgrid (:,:,:,i+1,:,:)
    case (5)
        ! The following line doesnt work because of compiler bug
        ! apgrid= (1-w)*this%apgrid(:,:,:,:,i,:) +w*this%apgrid(:,:,:,:,i+1,:)
        apgrid= this%apgrid(:,:,:,:,i,:) ! Because of compiler bug
        stocks= this%stocks(:,:,:,:,i,:)
        xgrid = this%xgrid (:,:,:,:,i,:)
        apgrid= (1-w)*apgrid +w*this%apgrid(:,:,:,:,i+1,:)
        stocks= (1-w)*stocks +w*this%stocks(:,:,:,:,i+1,:)
        xgrid = (1-w)*xgrid  +w*this%xgrid (:,:,:,:,i+1,:)
    case (6)
        apgrid= this%apgrid(:,:,:,:,:,i)
        stocks= this%stocks(:,:,:,:,:,i)
        xgrid = this%xgrid (:,:,:,:,:,i)
        apgrid= (1-w)*apgrid +w*this%apgrid(:,:,:,:,:,i+1)
        stocks= (1-w)*stocks +w*this%stocks(:,:,:,:,:,i+1)
        xgrid = (1-w)*xgrid  +w*this%xgrid (:,:,:,:,:,i+1)

    case default ! same as case 5
        apgrid= (1-w)*this%apgrid(:,:,:,:,i,:) +w*this%apgrid(:,:,:,:,i+1,:)
        stocks= (1-w)*this%stocks(:,:,:,:,i,:) +w*this%stocks(:,:,:,:,i+1,:)
        xgrid = (1-w)*this%xgrid (:,:,:,:,i,:) +w*this%xgrid (:,:,:,:,i+1,:)
    end select

    pol_int%apgrid = spread(apgrid, dim_x, 1)
    pol_int%stocks = spread(stocks, dim_x, 1)
    pol_int%xgrid  = spread(xgrid , dim_x, 1)

    call pol_int%calc_kappa
end function interpolate
!-------------------------------------------------------------------------------

end module policies_class
