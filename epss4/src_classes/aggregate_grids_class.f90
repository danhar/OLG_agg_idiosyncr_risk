! Copyright (C) 2016 Daniel Harenberg - All Rights Reserved
module aggregate_grids_class
    use kinds,  only : dp
    implicit none
    private

    public tAggGrids

    type tAggGrids
        real(dp), allocatable, dimension(:) :: k, mu ! Aggregate grids
        real(dp), private :: min_k=0.5_dp, max_k=16.0, min_mu = 0.0001_dp, max_mu = 0.12_dp, & ! min_mu = 0.01_dp
                             cover_k_l = 7.0, cover_k_u = 7.0, cover_mu_l = 4.0, cover_mu_u = 4.0, curv=1.75_dp
        logical :: fixed = .true.
    contains
        procedure :: allocate => allocate_grids
        procedure :: deallocate => deallocate_grids
        procedure :: set_params
        procedure :: read_unformatted
        procedure :: write_unformatted
        procedure :: construct =>construct_aggr_grid
        procedure :: update => update_grid_with_stats
    end type tAggGrids

contains

    elemental subroutine allocate_grids(this,nk,nmu)
        class(tAggGrids), intent(inout)  :: this
        integer,    intent(in)      :: nk,nmu
        call this%deallocate()
        allocate(this%k(nk),this%mu(nmu))
    end subroutine allocate_grids

    elemental subroutine deallocate_grids(this)
        class(tAggGrids), intent(inout)  :: this
        if (allocated(this%k)) deallocate(this%k)
        if (allocated(this%mu)) deallocate(this%mu)
    end subroutine deallocate_grids

    elemental subroutine set_params(this,k_min,k_max,mu_min,mu_max)
        class(tAggGrids), intent(inout)  :: this
        real(dp), intent(in) :: k_min,k_max,mu_min,mu_max
        this%min_k = k_min
        this%max_k = k_max
        this%min_mu= mu_min
        this%max_mu= mu_max
    end subroutine set_params

    subroutine read_unformatted(this,equilibrium_type, input_path)
        class(tAggGrids) ,intent(out) :: this
        character(len=*) ,intent(in)  :: equilibrium_type, input_path
        integer :: nk, nmu, io_stat

        open(55,file=input_path//'/aggr_grid_size_'//equilibrium_type//'.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
        read(55) nk, nmu
        close(55)

        if (io_stat == 0) then
            call this%allocate(nk,nmu)

            open(55,file=input_path//'/grids_'//equilibrium_type//'.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='read')
            read(55) this%k, this%mu
            close(55)
        endif

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR reading aggregate grids from unformatted file'
            stop 'STOP in in aggregate_grids_class:read_unformatted'
        endif

    end subroutine read_unformatted

    subroutine write_unformatted(this,equilibrium_type,input_path)
        class(tAggGrids) ,intent(in) :: this
        character(len=*) ,intent(in)  :: equilibrium_type, input_path
        integer :: io_stat

        open(55,file=input_path//'/aggr_grid_size_'//equilibrium_type//'.unformatted',form='unformatted',access='stream',iostat=io_stat, action='write')
        write(55) size(this%k), size(this%mu)
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR writing aggr_grid_size in aggregate_grids_class:write_unformatted'
        endif

        open(55,file=input_path//'/grids_'//equilibrium_type//'.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='write')
        write(55) this%k, this%mu
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR writing grids in aggregate_grids_class:write_unformatted'
        endif

    end subroutine write_unformatted

	elemental subroutine construct_aggr_grid(this, mean_guess, factor_k, factor_mu, cover_k, cover_mu, nk,nmu)
	! Create the grids for the aggregate variables k and mu
		use makegrid_mod

        class(tAggGrids) ,intent(inout):: this
        type(tAggGrids) ,intent(in) :: mean_guess
	    real(dp)        ,intent(in) :: factor_k, factor_mu, cover_k, cover_mu
	    integer         ,intent(in) :: nk, nmu
	    real(dp)					:: lb, ub		! lower and upper bound

        if (this%fixed) then
            this%k = MakeGrid(this%min_k ,this%max_k ,nk , this%curv)
            this%mu= MakeGrid(this%min_mu,this%max_mu,nmu, this%curv)
        else
            lb=factor_k*mean_guess%k(1)*(1.0-cover_k)
            if (lb < this%min_k) lb = this%min_k
            ub=factor_k*mean_guess%k(1)*(1.0+cover_k)
            this%k=MakeGrid(lb,ub,nk, this%curv)

            lb=factor_mu*mean_guess%mu(1)*(1.0-cover_mu)
            if (lb < this%min_mu) lb = this%min_mu
            ub=factor_mu*mean_guess%mu(1)*(1.0+cover_mu)
            this%mu=MakeGrid(lb,ub,nmu, this%curv)
        endif
	end subroutine construct_aggr_grid


    elemental subroutine update_grid_with_stats(this, k_min, k_max, mu_min, mu_max)
    ! Update the grid using min and max of k and mu
    ! Wanted to use statistics, only: tStats and have type(tStats) ,intent(in) :: k, mu,
    ! but circular dependency statistics-> aggregate_grids_class -> params_mod -> statistics.
        use makegrid_mod

        class(tAggGrids),intent(inout):: this
        real(dp) ,intent(in) :: k_min, k_max, mu_min, mu_max
        integer   :: n
        real(dp)                    :: lb, ub       ! lower and upper bound
        real(dp), parameter :: cover = 0.3_dp

        if (this%fixed) return  ! fixed grid, no update

        lb=k_min * (1.0 - cover)
        if (lb < this%min_k) lb = this%min_k
        ub=k_max * (1.0 + cover)
        n = size(this%k)
        this%k=MakeGrid(lb,ub,n, this%curv)

        lb=mu_min * (1.0 - cover)
        if (lb < this%min_mu) lb = this%min_mu
        ub=mu_max * (1.0 + cover)
        n = size(this%mu)
        this%mu=MakeGrid(lb,ub,n, this%curv)

        this%fixed = .true. ! only update once, in particular do not update when calc new GE after experiment.

    end subroutine update_grid_with_stats

    elemental subroutine update_grid_with_stats_old(this, k_mean, mu_mean, k_std, mu_std)
    ! Update the grid using mean and variance of k and mu
    ! Wanted to use statistics, only: tStats and have type(tStats) ,intent(in) :: k, mu,
    ! but circular dependency statistics-> aggregate_grids_class -> params_mod -> statistics.
        use makegrid_mod

        class(tAggGrids),intent(inout):: this
        real(dp) ,intent(in) :: k_mean, mu_mean, k_std, mu_std
        integer   :: n
        real(dp)                    :: lb, ub       ! lower and upper bound

        if (this%fixed) return  ! fixed grid, no update

        lb=k_mean - this%cover_k_l*k_std
        if (lb < this%min_k) lb = this%min_k
        ub=k_mean +this%cover_k_u*k_std
        n = size(this%k)
        this%k=MakeGrid(lb,ub,n, this%curv)

        lb=mu_mean -this%cover_mu_l*mu_std
        if (lb < this%min_mu) lb = this%min_mu
        ub=mu_mean +this%cover_mu_u*mu_std
        n = size(this%mu)
        this%mu=MakeGrid(lb,ub,n, this%curv)

        this%fixed = .true. ! only update once, in particular do not update when calc new GE after experiment.

    end subroutine update_grid_with_stats_old

end module aggregate_grids_class
