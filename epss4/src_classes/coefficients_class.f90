module coefficients_class
    use kinds,  only : dp
    implicit none
    private

    public tCoeffs

    type tCoeffs
        real(dp), dimension(:,:), allocatable :: k, mu, r_squared  ! coefficients for laws of motion, and residual sum of squares
        logical           :: normalize
        real(dp), dimension(:,:), allocatable, private :: k_initial, mu_initial
    contains
        procedure :: allocate
        procedure :: deallocate
        procedure :: read_unformatted
        procedure :: write_unformatted
        procedure :: write => write_results
        procedure :: save_initial_values
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
        write(55) size(this%k,1), size(this%k,2)
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

    subroutine write_results(this, unit, show_digits)
        use params_mod ,only: loms_in_logs, lom_k_version, lom_mu_version
        class(tCoeffs) ,intent(in) :: this
        integer        ,intent(in) :: unit, show_digits
        integer :: i
        character(:), allocatable :: fmt1, ind_var_1, ind_var_2 ! name of independent variables
        integer, parameter :: nl =11
        character(len=nl) :: dep_var ! name of dependent variables

        if (loms_in_logs) then
            dep_var   = " log(k')   "
            ind_var_1 = "log(k)"
            ind_var_2 = "log(k')"
        else
            dep_var   = "     k'    "
            ind_var_1 = "k"
            ind_var_2 = "k'"
        endif

        select case(lom_k_version)
        case(1)
            write(unit,166) dep_var, 'constant', ind_var_1                                                  ,'R^2'
        case(2)
            write(unit,166) dep_var, 'constant', ind_var_1, ind_var_1//'^2'                                 ,'R^2'
        case(3)
            write(unit,168) dep_var, 'constant', ind_var_2, ind_var_2//'^2', 'mu'                           ,'R^2'
        case(4)
            write(unit,166) dep_var, 'constant', ind_var_1, 'mu'                                            ,'R^2'
        case(5)
            write(unit,166) dep_var, 'constant', ind_var_1, 'mu', ind_var_1//'*mu'                          ,'R^2'
        case(6)
            write(unit,166) dep_var, 'constant', ind_var_1, 'mu', ind_var_1//'^2', 'mu^2'                   ,'R^2'
        case(7)
            write(unit,166) dep_var, 'constant', ind_var_1, 'mu', ind_var_1//'^2', 'mu^2', ind_var_1//'*mu' ,'R^2'
        case default
            write(unit,166) dep_var, 'constant', ind_var_1, ind_var_1//'^2'                                 ,' R^2  '
        end select
166     format(a<nl>,tr1,a<show_digits+6>,<size(this%k,1)-1>a<show_digits+10>,a<show_digits+8>)

        do i=1,size(this%k,2)
            write(unit,167) " z  =",i, this%k(:,i), this%r_squared(1,i)
        enddo
167     format(a5,i2,tr4,<size(this%k,1)>(es<show_digits+7>.<show_digits>,3x),'|',f<show_digits+4>.<show_digits+2>)
        write(unit,*)

        dep_var   = "    mu'  |"
        select case(lom_mu_version)
        case(1)
            write(unit,168) dep_var, 'constant', ind_var_2                                                  ,'R^2'
        case(2)
            write(unit,168) dep_var, 'constant', ind_var_2, ind_var_2//'^2'                                 ,'R^2'
        case(3)
            write(unit,168) dep_var, 'constant', ind_var_2, ind_var_2//'^2', 'mu'                           ,'R^2'
        case(4)
            write(unit,168) dep_var, 'constant', ind_var_2, 'mu'                                            ,'R^2'
        case(5)
            write(unit,168) dep_var, 'constant', ind_var_2, 'mu', ind_var_2//'*mu'                          ,'R^2'
        case(6)
            write(unit,168) dep_var, 'constant', ind_var_2, 'mu', ind_var_2//'^2', 'mu^2'                   ,'R^2'
        case(7)
            write(unit,168) dep_var, 'constant', ind_var_2, 'mu', ind_var_2//'^2', 'mu^2', ind_var_2//'*mu' ,'R^2'
        case default
            write(unit,168) dep_var, 'constant', ind_var_2, ind_var_2//'^2'                                 ,'R^2'
        end select
168     format(a<nl>,tr1,a<show_digits+6>,<size(this%mu,1)-1>a<show_digits+10>,a<show_digits+8>)

        do i=1,size(this%mu,2)
            write(unit,169) " z' =",i, this%mu(:,i), this%r_squared(2,i)
        enddo
169     format(a5,i2,tr4,<size(this%mu,1)>(es<show_digits+7>.<show_digits>,3x),'|',f<show_digits+4>.<show_digits+2>)
        write(unit,*)

    end subroutine write_results
!-------------------------------------------------------------------------------

    pure subroutine save_initial_values(this)
        class(tCoeffs)         ,intent(inout):: this
        this%k_initial = this%k
        this%mu_initial = this%mu
    end subroutine save_initial_values
!-------------------------------------------------------------------------------

    pure subroutine maketype(this, coeff_vec)
        use params_mod ,only: n_coeffs, pooled_regression, pi1_delta, nz
        class(tCoeffs)         ,intent(inout):: this
        real(dp), dimension(:) ,intent(in) :: coeff_vec
        integer :: zc

        if (.not. allocated(this%k)) call this%allocate(n_coeffs,nz)
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
            this%mu(:,2:3) = this%k(:,2:3)
            this%mu(:,4)= coeff_vec(3*n_coeffs+1:4*n_coeffs)
        else
            do zc=1, nz
                this%k(:,zc)  = coeff_vec(n_coeffs*(zc-1)+1:n_coeffs*zc)
                this%mu(:,zc) = coeff_vec(n_coeffs*(zc-1)+nz*n_coeffs+1:n_coeffs*zc+nz*n_coeffs)
            enddo
        endif

        if (this%normalize) then
            this%mu = this%mu*this%mu_initial
            this%k = this%k*this%k_initial
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
            coeffs%mu = this%mu/this%mu_initial
            coeffs%k = this%k/this%k_initial
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
