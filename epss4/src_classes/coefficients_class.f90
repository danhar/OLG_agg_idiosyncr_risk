module coefficients_class
    use kinds,  only : dp
    implicit none
    private

    public tCoeffs

    type tCoeffs
        real(dp), dimension(:,:), allocatable :: k, mu, r_squared  ! coefficients for laws of motion, and residual sum of squares
        logical           :: normalize
        real(dp), private :: norm_factor_small = 10.0, norm_factor_large = 100.0 ! Could normalize everytime with a different value (e.g. I could normalize all coeffs = 1.0 every time in the krusell-smith alg)
    contains
        procedure :: allocate
        procedure :: deallocate
        procedure :: read_unformatted
        procedure :: write_unformatted
        procedure :: maketype
        procedure :: makevector
    end type tCoeffs

contains

    elemental subroutine allocate(this,ncoeffs,nz)
        class(tCoeffs), intent(out)  :: this
        integer,    intent(in)      :: ncoeffs,nz
        allocate(this%k(ncoeffs,nz),this%mu(ncoeffs,nz))
        allocate(this%r_squared(2,nz))
    end subroutine allocate
!-------------------------------------------------------------------------------

    elemental subroutine deallocate(this)
        class(tCoeffs), intent(out)  :: this
    end subroutine deallocate
!-------------------------------------------------------------------------------

    subroutine read_unformatted(this)
        class(tCoeffs) ,intent(out) :: this
        integer :: ncoeffs,nz , io_stat

        open(55,file='model_input/last_results/coeffs_size.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
        read(55) ncoeffs, nz
        close(55)

        if (io_stat == 0) then
            call this%allocate(ncoeffs,nz)

            open(55,file='model_input/last_results/coeffs_ge.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='read')
            read(55) this%k, this%mu, this%r_squared
            close(55)
        endif

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR reading coefficients from unformatted file'
            stop 'STOP in in coefficients_class:read_unformatted'
        endif

    end subroutine read_unformatted
!-------------------------------------------------------------------------------

    subroutine write_unformatted(this)
        class(tCoeffs) ,intent(in) :: this
        integer :: io_stat

        open(55,file='model_input/last_results/coeffs_size.unformatted',form='unformatted',access='stream',iostat=io_stat, action='write')
        write(55) size(this%k)
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR in aggregate_grids_class:write_unformatted'
        endif

        open(55,file='model_input/last_results/coeffs_ge.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='write')
        write(55) this%k, this%mu, this%r_squared
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR in aggregate_grids_class:write_unformatted'
        endif

    end subroutine write_unformatted
!-------------------------------------------------------------------------------

    pure subroutine maketype(this, coeff_vec)
        use params_mod ,only: n_coeffs, pooled_regression, pi1_delta, nz
        class(tCoeffs)         ,intent(out):: this
        real(dp), dimension(:) ,intent(in) :: coeff_vec
        integer :: zc

        call this%allocate(n_coeffs,nz)
        this%r_squared = 0.0

        if (pooled_regression) then
            do zc=1, nz
                this%k(:,zc) = coeff_vec(1:n_coeffs)
                this%mu(:,zc) = coeff_vec(n_coeffs+1:)
            enddo
        elseif (pi1_delta==1.0) then ! Alternatively, could have nz_1=count(stat_dist_z==0.0)
            this%k(:,1) = coeff_vec(1:n_coeffs)
            this%k(:,2) = 0.0
            this%k(2,2) = 1.0
            this%k(:,3) = this%k(:,2)
            this%k(:,4) = coeff_vec(n_coeffs+1:2*n_coeffs)
            this%mu(:,1)= coeff_vec(2*n_coeffs+1:3*n_coeffs)
            this%mu(:,2) = 0.0
            this%mu(:,3) = this%mu(:,2)
            this%mu(:,4)= coeff_vec(3*n_coeffs+1:4*n_coeffs)
        else
            do zc=1, nz
                this%k(:,zc)  = coeff_vec(n_coeffs*(zc-1)+1:n_coeffs*zc)
                this%mu(:,zc) = coeff_vec(n_coeffs*(zc-1)+nz*n_coeffs+1:n_coeffs*zc+nz*n_coeffs)
            enddo
        endif

        if (this%normalize) then
            this%mu(1,:) = this%mu(1,:)/this%norm_factor_small
            this%mu(2,:) = this%mu(2,:)/this%norm_factor_large
            if (n_coeffs >= 3) then
                this%k(3,:) = this%k(3,:)/this%norm_factor_small
                this%mu(3,:) = this%mu(3,:)/this%norm_factor_large
            endif
        endif

    end subroutine maketype
!-------------------------------------------------------------------------------

    pure function makevector(this) result(coeff_vec)
        use params_mod ,only: pooled_regression, pi1_delta
        real(dp), dimension(:), allocatable :: coeff_vec
        class(tCoeffs) ,intent(in) :: this
        type(tCoeffs) :: coeffs     ! need this only because this has intent(in), and if I want to normalize
        integer :: zc, nz, n_coeffs

        n_coeffs = size(this%k,1)
        nz       = size(this%k,2)
        coeffs = this

        if (this%normalize) then
            coeffs%mu(1,:) = coeffs%mu(1,:)*this%norm_factor_small
            coeffs%mu(2,:) = coeffs%mu(2,:)*this%norm_factor_large
            if (n_coeffs >= 3) then
                coeffs%k(3,:) = coeffs%k(3,:)*this%norm_factor_small
                coeffs%mu(3,:) = coeffs%mu(3,:)*this%norm_factor_large
            endif
        endif

        if (pooled_regression) then
            allocate(coeff_vec(2*n_coeffs))
            coeff_vec(1:n_coeffs) = coeffs%k(:,1)
            coeff_vec(n_coeffs+1:)= coeffs%mu(:,1)
        elseif (pi1_delta==1.0) then ! Alternatively, could have nz_1=count(stat_dist_z==0.0)
            allocate(coeff_vec(2*n_coeffs*nz/2))
            coeff_vec(1:n_coeffs) = coeffs%k(:,1)
            coeff_vec(n_coeffs+1:2*n_coeffs)= coeffs%k(:,4)
            coeff_vec(2*n_coeffs+1:3*n_coeffs)= coeffs%mu(:,1)
            coeff_vec(3*n_coeffs+1:4*n_coeffs)= coeffs%mu(:,4)
        else
            allocate(coeff_vec(2*n_coeffs*nz))
            do zc=1, nz
                coeff_vec(n_coeffs*(zc-1)+1:n_coeffs*zc) = coeffs%k(:,zc)
                coeff_vec(n_coeffs*(zc-1)+nz*n_coeffs+1:n_coeffs*zc+nz*n_coeffs)= coeffs%mu(:,zc)
            enddo
        endif

    end function makevector
!-------------------------------------------------------------------------------

end module coefficients_class
