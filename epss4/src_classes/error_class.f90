!*******************************************************************************
! Copyright (c) 2016 Daniel Harenberg - All rights reserved.
!*******************************************************************************

module error_class
    use kinds
    implicit none
    private

    public tErrors

    type tErrors
	    logical(1) ,dimension(:,:,:,:,:,:)   ,allocatable :: asset
	    logical(1) ,dimension(:,:,:,:,:,:,:) ,allocatable :: cons    ! true if xp(= cah tomorrow) < xgrid(1,1,1,jc+1) in policyfunctions:solve_policyfunctions
        logical(1) ,dimension(:,:,:)         ,allocatable :: kp, mup, rfp ! true if forecast out of grid
        logical :: not_converged
    contains
        private
        procedure :: write2file_errs
        procedure :: write2file_not_converged
        procedure :: print_error_msg_tErr
        procedure :: print_error_msg_not_converged
        procedure :: print_errors
        procedure, public :: allocate => allocate_errs
        procedure, public :: deallocate => deallocate_errs
        generic, public :: write2file => write2file_errs, write2file_not_converged
        generic, public :: print2stderr => print_error_msg_tErr, print_error_msg_not_converged, print_errors
    end type tErrors

contains
!------------------------------------------------------------------------------------------------------------------------------
    pure subroutine allocate_errs(this,nk,nmu)
    ! Constructor (maybe should follow Rouson/Xia/Xu, p.37 austronaut example)
        use params_mod, only: nx,n_eta,nz,nj
        class(tErrors), intent(out)  :: this
        integer,    intent(in)      :: nk,nmu

        if (allocated(this%cons)) call this%deallocate ! Compiler bug, can remove later

        allocate(this%cons(nz,nx,n_eta,nz,nj,nk,nmu),this%asset(nx,n_eta,nz,nj,nk,nmu))
        allocate(this%kp(nz,nk,nmu), this%mup(nz,nk,nmu), this%rfp(nz,nk,nmu))
        this%cons  = .false.
        this%asset = .false.
        this%kp    = .false.
        this%mup   = .false.
        this%rfp   = .false.
        this%not_converged = .false.
    end subroutine allocate_errs

    pure subroutine deallocate_errs(this)
        class(tErrors), intent(inout)  :: this
        if (allocated(this%cons)) deallocate(this%cons)
        if (allocated(this%asset)) deallocate(this%asset)
        if (allocated(this%kp)) deallocate(this%kp)
        if (allocated(this%mup)) deallocate(this%mup)
        if (allocated(this%rfp)) deallocate(this%rfp)
    end subroutine deallocate_errs
!------------------------------------------------------------------------------------------------------------------------------

    subroutine print_errors(this, max_errs)
    ! Print errors to standard output
        class(tErrors), intent(in)    :: this
        integer, intent(in):: max_errs
        integer :: count_errs, muc, kc, jc, zc, ec, xc, zpc

        count_errs=0
        do muc=1, size(this%asset,6)
            do kc=1, size(this%asset,5)
                do zc = 1, size(this%asset,3)
	                if (this%kp(zc,kc, muc)) then
	                    count_errs = count_errs+1
	                    if (count_errs> max_errs) go to 999
	                    print '(a86,3i3)','ERROR: policyfunctions:calc_vars_tomorrow: forecast kp out of grid at (zc,kc,muc)=', zc,kc, muc
	                endif
	                if (this%mup(zc, kc, muc)) then
	                    count_errs = count_errs+1
	                    if (count_errs> max_errs) goto 999
	                    print '(a87,3i3)','ERROR: policyfunctions:calc_vars_tomorrow: forecast mup out of grid at (zc,kc,muc)=', zc, kc, muc
	                endif
	                if (this%rfp(zc, kc, muc)) then
                        count_errs = count_errs+1
                        if (count_errs> max_errs) goto 999
                        print '(a74,3i3)','ERROR: policyfunctions:calc_vars_tomorrow: rfp < rp(1) at (zc,kc,muc)=', zc, kc, muc
                    endif
                    do jc = 1, size(this%asset,4)
                        do ec =1, size(this%asset,2)
	                        do xc =1, size(this%asset,1)
		                        if (this%asset(xc,ec,zc,jc,kc,muc)) then
			                        count_errs = count_errs+1
			                        if (count_errs> max_errs) goto 999
			                        print '(a91,6i3)','ERROR: policyfunctions:asset_allocation: couldnt bracket kappa at (xc,ec,zc,jc,kc,muc)=', xc,ec,zc,jc, kc, muc
		                        endif
		                        do zpc = 1, size(this%cons,1)
			                        if (this%cons(zpc,xc,ec,zc,jc,kc,muc)) then
			                            count_errs = count_errs+1
			                            if (count_errs> max_errs) goto 999
			                            print '(a82,6i3)','ERROR: policyfunctions:consumption: cons_interp<=0.0 at (zpc,xc,ec,zc,jc,kc,muc)=', zpc, xc,ec,zc,jc, kc, muc
			                        endif
			                    enddo
		                    enddo
	                    enddo
                    enddo
                enddo
            enddo
        enddo

        return
999     print*, 'Showed the first ', max_errs,' errors, but there were more...'

    end subroutine print_errors
!------------------------------------------------------------------------------------------------------------------------------

    subroutine print_error_msg_tErr(this)
    ! Print error summary statistics
        class(tErrors), intent(in)    :: this
	    integer :: count_err
	    real(dp) :: perc_err

	    if (any(this%kp)) then
	        count_err = count(this%kp)
	        perc_err  = real(count_err,dp)/real(size(this%kp),dp)*100.0
	        print 214, ' Warning: policyfunctions: # kp  not in grid =', count_err,'  (',perc_err,'%)'
	    endif

        if (any(this%mup)) then
            count_err = count(this%mup)
            perc_err  = real(count_err,dp)/real(size(this%mup),dp)*100.0
            print 214, ' Warning: policyfunctions: # mup not in grid =', count_err,'  (',perc_err,'%)'
        endif

        if (any(this%cons)) then
            count_err = count(this%cons)
            perc_err  = real(count_err,dp)/real(size(this%cons),dp)*100.0
            print 214, ' Warning: policyfunctions: # cons_interp < 0 =', count_err,'  (',perc_err,'%)'
        endif

        if (any(this%asset)) then
            count_err = count(this%asset)
            perc_err  = real(count_err,dp)/real(size(this%asset),dp)*100.0
            print 214, ' Warning: policyfunctions: # kappa not found =', count_err,'  (',perc_err,'%)'
        endif

214     format((a,i6,a,f5.1,a2))

    end subroutine print_error_msg_tErr
!------------------------------------------------------------------------------------------------------------------------------

    subroutine write2file_errs(this, path)
    ! Write errors to file
        class(tErrors), intent(in)    :: this
        character(len=*), intent(in) :: path

        if (any(this%cons)) then
            open(unit=31, file=path//'/err_cons.txt', status = 'replace')
            write(31,*) this%cons
            close(31)
        endif

        if (any(this%asset)) then
            open(unit=32, file=path//'/err_asset.txt', status = 'replace')
            write(32,*) this%asset
            close(32)
        endif

        if (any(this%kp)) then
            open(unit=32, file=path//'/err_kp.txt', status = 'replace')
            write(32,*) this%kp
            close(32)
        endif

        if (any(this%mup)) then
            open(unit=32, file=path//'/err_mup.txt', status = 'replace')
            write(32,*) this%mup
            close(32)
        endif

        if (any(this%rfp)) then
            open(unit=32, file=path//'/err_rfp.txt', status = 'replace')
            write(32,*) this%rfp
            close(32)
        endif
43  format(a86)
    end subroutine write2file_errs
!------------------------------------------------------------------------------------------------------------------------------

    subroutine print_error_msg_not_converged(this, dir)
        class(tErrors), intent(in)    :: this
        character(len=*), intent(in) :: dir
        print*, 'CRITICAL ERROR: main: root finder did not converge in '//dir
    end subroutine print_error_msg_not_converged


    subroutine write2file_not_converged(this, fvals, path)
    ! Write error to file if not_converged
        class(tErrors), intent(in)    :: this
        real(dp), dimension(:), intent(in) :: fvals
        character(len=*), intent(in) :: path
        open(unit=301, file=path//'/err_rootfinder.txt', status = 'replace')
        write(301,*) 'this%not_converged = ', this%not_converged
        write(301,*) 'fvals = ', fvals
        close(301)
    end subroutine write2file_not_converged

end module error_class
