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
        integer :: nk, nmu

!        open(55,file='model_input/last_results/aggr_grid_size_'//equilibrium_type//'.unformatted',form='unformatted',action='read')
!        read(55) nk, nmu
!        close(55)
        if (equilibrium_type == 'ms') then ! DELETE LATER!
            nk = 1
            nmu = 1
        else
            nk = 10
            nmu = 8
        endif

        call this%allocate(nk,nmu)

        open(55,file='model_input/last_results/grids_'//equilibrium_type//'.unformatted'  ,form='unformatted',action='read')
        read(55) this%k, this%mu
        close(55)
    end subroutine read_unformatted

    subroutine write_unformatted(this,equilibrium_type)
        class(tAggGrids) ,intent(in) :: this
        character(len=*) ,intent(in)  :: equilibrium_type

        open(55,file='model_input/last_results/aggr_grid_size_'//equilibrium_type//'.unformatted',form='unformatted',action='write')
        write(55) size(this%k), size(this%mu)
        close(55)

        open(55,file='model_input/last_results/grids_'//equilibrium_type//'.unformatted'  ,form='unformatted',action='write')
        write(55) this%k, this%mu
        close(55)
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

end module aggregate_grids_class
