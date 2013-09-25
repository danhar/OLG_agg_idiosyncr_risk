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
        generic :: allocate => allocate_ncoeffs, allocate_version
        procedure ,private :: allocate_ncoeffs
        procedure ,private :: allocate_version
        procedure ,private :: set_default_values
        procedure :: deallocate
        procedure :: read_unformatted
        procedure :: write_unformatted
        procedure :: write => write_results
        procedure :: save_initial_values
        procedure :: maketype
        procedure :: makevector
    end type tCoeffs

contains

    elemental subroutine allocate_ncoeffs(this,ncoeffs_k, ncoeffs_mu_o ,nz)
        class(tCoeffs) ,intent(out)   :: this
        integer ,intent(in)           :: ncoeffs_k ,nz
        integer ,intent(in) ,optional :: ncoeffs_mu_o
        integer :: ncoeffs_mu

        if (present(ncoeffs_mu_o)) then
            ncoeffs_mu = ncoeffs_mu_o
        else
            ncoeffs_mu = ncoeffs_k
        endif

        allocate(this%k(ncoeffs_k,nz),this%mu(ncoeffs_mu,nz))
        allocate(this%r_squared(2,nz))

        call this%set_default_values()
    end subroutine allocate_ncoeffs

    elemental subroutine allocate_version(this)
        use params_mod ,only: lom_k_version, lom_mu_version, nz
        class(tCoeffs), intent(out)  :: this
        integer :: ncoeffs_k, ncoeffs_mu

        select case(lom_k_version)
        case(1)
            ncoeffs_k = 2
        case(2)
            ncoeffs_k = 3
        case(3)
            ncoeffs_k = 4
        case(4)
            ncoeffs_k = 3
        case(5)
            ncoeffs_k = 4
        case(6)
            ncoeffs_k = 5
        case(7)
            ncoeffs_k = 6
        case default
            ncoeffs_k = 3
        end select

        select case(lom_mu_version)
        case(1)
            ncoeffs_mu = 2
        case(2)
            ncoeffs_mu = 3
        case(3)
            ncoeffs_mu = 4
        case(4)
            ncoeffs_mu = 3
        case(5)
            ncoeffs_mu = 4
        case(6)
            ncoeffs_mu = 5
        case(7)
            ncoeffs_mu = 6
        case default
            ncoeffs_mu = 3
        end select

        call this%allocate_ncoeffs(ncoeffs_k, ncoeffs_mu,nz)

    end subroutine allocate_version
!-------------------------------------------------------------------------------

    elemental subroutine deallocate(this)
        class(tCoeffs), intent(out)  :: this
    end subroutine deallocate
!-------------------------------------------------------------------------------

    elemental subroutine set_default_values(this)
        class(tCoeffs) ,intent(inout)   :: this
        this%k=0.0
        this%mu=0.0
        this%r_squared=0.0
        this%normalize=.true.
    end subroutine set_default_values
!-------------------------------------------------------------------------------

    subroutine read_unformatted(this,input_path)
        use params_mod, only: params_set
        class(tCoeffs) ,intent(out) :: this
        character(*)   ,intent(in)  :: input_path
        integer :: ncoeffs_k, ncoeffs_mu,nz , io_stat, io_stat2, lom_k_version_read, lom_mu_version_read

        open(55,file=input_path//'/coeffs_size.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
        read(55,iostat=io_stat2) ncoeffs_k, ncoeffs_mu, nz, lom_k_version_read, lom_mu_version_read
        close(55)

        if (io_stat2 == 0) then
            call params_set('lom_k_version', lom_k_version_read)
            call params_set('lom_mu_version', lom_mu_version_read)
        else
            open(55,file=input_path//'/coeffs_size.unformatted',form='unformatted',access='stream',iostat=io_stat,action='read')
            read(55,iostat=io_stat2) ncoeffs_k, nz
            close(55)
            ncoeffs_mu = ncoeffs_k
            print*, 'WARNING: couldnt read ncoeffs_mu, lom_k_version.'
            print*, 'Setting ncoeffs_mu = ncoeffs_k and using lom_k_version from calib file.'
        endif


        if (io_stat == 0 .and. io_stat2 == 0) then
            call this%allocate(ncoeffs_k, ncoeffs_mu, nz)

            open(55,file=input_path//'/coeffs_ge.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='read')
            read(55) this%k, this%mu, this%r_squared
            close(55)
        endif

        if (io_stat .ne. 0 .or. io_stat2 .ne. 0) then
            print*, 'I/O ERROR reading coefficients from unformatted file'
            stop 'STOP in in coefficients_class:read_unformatted'
        endif

    end subroutine read_unformatted
!-------------------------------------------------------------------------------

    subroutine write_unformatted(this, input_path)
        use params_mod ,only: lom_k_version, lom_mu_version
        class(tCoeffs) ,intent(in) :: this
        character(*)   ,intent(in) :: input_path
        integer :: io_stat

        open(55,file=input_path//'/coeffs_size.unformatted',form='unformatted',access='stream',iostat=io_stat, action='write')
        write(55) size(this%k,1), size(this%mu,1), size(this%k,2), lom_k_version, lom_mu_version
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR writing coeffs_size in coefficients_class:write_unformatted'
        endif

        open(55,file=input_path//'/coeffs_ge.unformatted'  ,form='unformatted',access='stream',iostat=io_stat,action='write')
        write(55) this%k, this%mu, this%r_squared
        close(55)

        if (io_stat .ne. 0) then
            print*, 'I/O ERROR writing coeffs_ge in coefficients_class:write_unformatted'
        endif

    end subroutine write_unformatted
!-------------------------------------------------------------------------------

    subroutine write_results(this, unit, show_digits)
        use params_mod ,only: loms_in_logs, lom_k_version, lom_mu_version
        class(tCoeffs) ,intent(in) :: this
        integer        ,intent(in) :: unit, show_digits
        integer :: i
        character(:), allocatable :: ind_var_1, ind_var_2 ! name of independent variables
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
            write(unit,166) dep_var, 'constant', ind_var_1, ind_var_1//'^2'                                 ,'R^2'
        end select
166     format(a<nl>,tr1,a<show_digits+6>,<size(this%k,1)-1>a<show_digits+10>,a<show_digits+8>)

        do i=1,size(this%k,2)
            write(unit,167) " z  =",i, this%k(:,i), this%r_squared(1,i)
        enddo
167     format(a5,i2,tr4,<size(this%k,1)>(es<show_digits+7>.<show_digits>,3x),'|',f<show_digits+4>.<show_digits+2>)
        write(unit,*)

        dep_var   = "    mu'   "
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
        use params_mod ,only: pooled_regression, pi1_delta
        class(tCoeffs)         ,intent(inout):: this
        real(dp), dimension(:) ,intent(in) :: coeff_vec
        integer :: zc, n_coeffs_k, n_coeffs_mu, nz

        if (.not. allocated(this%k)) call this%allocate()
        ! this%r_squared = 0.0

        n_coeffs_k  = size(this%k ,1)
        n_coeffs_mu = size(this%mu,1)
        nz          = size(this%k ,2)

        if (pooled_regression) then
            do zc=1, nz
                this%k(:,zc) = coeff_vec(1:n_coeffs_k)
                this%mu(:,zc) = coeff_vec(n_coeffs_k+1:)
            enddo
        elseif (pi1_delta==1.0) then ! Alternatively, could have nz_1=count(stat_dist_z==0.0)
            this%k(:,1) = coeff_vec(1:n_coeffs_k)
            this%k(:,2) = 0.0
            this%k(2,2) = 1.0
            this%k(:,3) = this%k(:,2)
            this%k(:,4) = coeff_vec(n_coeffs_k+1:2*n_coeffs_k)
            this%mu(:,1)= coeff_vec(2*n_coeffs_k+1:2*n_coeffs_k + n_coeffs_mu)
            this%mu(:,2:3) = this%k(:,2:3)
            this%mu(:,4)= coeff_vec(2*n_coeffs_k + n_coeffs_mu+1:2*(n_coeffs_k + n_coeffs_mu))
        else
            do zc=1, nz
                this%k(:,zc)  = coeff_vec(n_coeffs_k*(zc-1)+1:n_coeffs_k*zc)
                this%mu(:,zc) = coeff_vec(n_coeffs_mu*(zc-1)+nz*n_coeffs_k+1:n_coeffs_mu*zc+nz*n_coeffs_k)
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
        integer :: zc, nz, n_coeffs_k, n_coeffs_mu

        n_coeffs_k  = size(this%k ,1)
        n_coeffs_mu = size(this%mu,1)
        nz          = size(this%k,2)
        coeffs      = this

        if (this%normalize) then
            coeffs%mu = this%mu/this%mu_initial
            coeffs%k = this%k/this%k_initial
        endif

        if (pooled_regression) then
            allocate(coeff_vec(n_coeffs_k + n_coeffs_mu))
            coeff_vec(1:n_coeffs_k) = coeffs%k(:,1)
            coeff_vec(n_coeffs_k+1:)= coeffs%mu(:,1)
        elseif (pi1_delta==1.0) then ! Alternatively, could have nz_1=count(stat_dist_z==0.0)
            allocate(coeff_vec((n_coeffs_k + n_coeffs_mu)*nz/2))
            coeff_vec(1:n_coeffs_k) = coeffs%k(:,1)
            coeff_vec(n_coeffs_k+1:2*n_coeffs_k)= coeffs%k(:,4)
            coeff_vec(2*n_coeffs_k+1:2*n_coeffs_k+n_coeffs_mu)= coeffs%mu(:,1)
            coeff_vec(2*n_coeffs_k+n_coeffs_mu+1:2*(n_coeffs_k+n_coeffs_mu))= coeffs%mu(:,4)
        else
            allocate(coeff_vec((n_coeffs_k+n_coeffs_mu)*nz))
            do zc=1, nz
                coeff_vec(n_coeffs_k *(zc-1)+1:n_coeffs_k*zc) = coeffs%k(:,zc)
                coeff_vec(n_coeffs_mu*(zc-1)+nz*n_coeffs_k+1:n_coeffs_mu*zc+nz*n_coeffs_k)= coeffs%mu(:,zc)
            enddo
        endif

    end function makevector
!-------------------------------------------------------------------------------

end module coefficients_class
