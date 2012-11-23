module aggregate_grids_class
    use kinds,  only : dp
    implicit none
    private

    public tAggGrids

    type tAggGrids
        real(dp), allocatable, dimension(:) :: k, mu ! Aggregate grids
    contains
        procedure :: allocate => allocate_grids
        procedure :: deallocate => deallocate_grids
        procedure :: read_unformatted
        procedure :: write_unformatted
        procedure :: construct =>construct_aggr_grid
        procedure :: update => update_grid_with_stats
    end type tAggGrids

contains

    elemental subroutine allocate_grids(this,nk,nmu)
        class(tAggGrids), intent(out)  :: this
        integer,    intent(in)      :: nk,nmu
        allocate(this%k(nk),this%mu(nmu))
    end subroutine allocate_grids

    elemental subroutine deallocate_grids(this)
        class(tAggGrids), intent(out)  :: this
    end subroutine deallocate_grids

    subroutine read_unformatted(this,equilibrium_type)
        class(tAggGrids) ,intent(out) :: this
        character(len=*) ,intent(in)  :: equilibrium_type
        integer :: nk, nmu, io_stat

        open(55,file='model_input/last_results/aggr_grid_size_'//equilibrium_type//'.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
        read(55) nk, nmu
        close(55)

        if (io_stat == 0) then
            call this%allocate(nk,nmu)

            open(55,file='model_input/last_results/grids_'//equilibrium_type//'.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='read')
            read(55) this%k, this%mu
            close(55)
        endif

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR reading aggregate grids from unformatted file'
            stop 'STOP in in aggregate_grids_class:read_unformatted'
        endif

    end subroutine read_unformatted

    subroutine write_unformatted(this,equilibrium_type)
        class(tAggGrids) ,intent(in) :: this
        character(len=*) ,intent(in)  :: equilibrium_type
        integer :: io_stat

        open(55,file='model_input/last_results/aggr_grid_size_'//equilibrium_type//'.unformatted',form='unformatted',access='stream',iostat=io_stat, action='write')
        write(55) size(this%k), size(this%mu)
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR in aggregate_grids_class:write_unformatted'
        endif

        open(55,file='model_input/last_results/grids_'//equilibrium_type//'.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='write')
        write(55) this%k, this%mu
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR in aggregate_grids_class:write_unformatted'
        endif

    end subroutine write_unformatted

	elemental subroutine construct_aggr_grid(this, mean_guess, factor_k, factor_mu, cover_k, cover_mu, nk,nmu)
	! Create the grids for the aggregate variables k and mu
		use makegrid_mod

        class(tAggGrids) ,intent(out):: this
        type(tAggGrids) ,intent(in) :: mean_guess
	    real(dp)        ,intent(in) :: factor_k, factor_mu, cover_k, cover_mu
	    integer         ,intent(in) :: nk, nmu
	    real(dp)					:: lb, ub		! lower and upper bound

        call this%allocate(nk,nmu)

	    lb=factor_k*mean_guess%k(1)*(1.0-cover_k)
	    ub=factor_k*mean_guess%k(1)*(1.0+cover_k)
	    this%k=MakeGrid(lb,ub,nk)

	    lb=factor_mu*mean_guess%mu(1)*(1.0-cover_mu)
	    ub=factor_mu*mean_guess%mu(1)*(1.0+cover_mu)
	    this%mu=MakeGrid(lb,ub,nmu)
	end subroutine construct_aggr_grid

    elemental subroutine update_grid_with_stats(this, k_mean, mu_mean, k_std, mu_std)
    ! Update the grid using mean and variance of k and mu
    ! Wanted to use statistics, only: tStats and have type(tStats) ,intent(in) :: k, mu,
    ! but circular dependency statistics-> aggregate_grids_class -> params_mod -> statistics
        use makegrid_mod

        class(tAggGrids),intent(inout):: this
        real(dp) ,intent(in) :: k_mean, mu_mean, k_std, mu_std
        integer   :: n
        real(dp)                    :: lb, ub       ! lower and upper bound
        real(dp), parameter :: cover_k_l = 2.0, cover_k_u = 3.0,cover_mu_l = 2.0, cover_mu_u = 2.0, min_k=1.0_dp, min_mu = 0.01_dp

        lb=k_mean -cover_k_l*k_std
        if (lb < min_k) lb = min_k
        ub=k_mean +cover_k_u*k_std
        n = size(this%k)
        this%k=MakeGrid(lb,ub,n)

        lb=mu_mean -cover_mu_l*mu_std
        if (lb < min_mu) lb = min_mu
        ub=mu_mean +cover_mu_u*mu_std
        n = size(this%mu)
        this%mu=MakeGrid(lb,ub,n)
    end subroutine update_grid_with_stats

end module aggregate_grids_class
