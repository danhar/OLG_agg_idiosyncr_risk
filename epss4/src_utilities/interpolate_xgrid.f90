module interpolate_xgrid
! This interpolation is specific, because it first creates the new xgrid, then acts on
! various objects. More elegant would be to have a function create xgrid, and have
! interpolation functions for each object.
! This module is USEd only in krusell_smith, the function InterpolateXgrid is only executed IF (exogenous_xgrid).
    implicit none
contains
    pure subroutine InterpolateXgrid(policies, value, polx, valx)
        use kinds      ,only: dp
        use policies_class
        use fun_lininterp

        type(tPolicies) ,intent(in)    :: policies
        real(dp)        ,intent(in)    :: value(:,:,:,:,:,:)
        type(tPolicies) ,intent(out)   :: polx
        real(dp) ,allocatable ,intent(out)   :: valx(:,:,:,:,:,:)
        integer                        :: muc, kc, jc, zc, ec, nz, nk, nmu

        nz = size(value,3); nk = size(value,5); nmu=size(value,6)
        call polx%allocate(nz, nk,nmu)
        valx = value
        do muc = 1, nmu
            ! mean of first and last mu seems good, coz small aggregate grid and policies linear in that dimension
            polx%xgrid(:,:,:,:,:,muc) = 0.5_dp * policies%xgrid(:,:,:,:,:,1) + 0.5_dp * policies%xgrid(:,:,:,:,:,nmu)
            do kc = 1, nk
                do jc=1,size(value,4)
                    do zc = 1,size(value,3)
                        do ec = 1,size(value,2)
	                        valx(:,ec,zc,jc,kc,muc)        = f_lininterp(policies%xgrid(:,ec,zc,jc,kc,muc), value(:,ec,zc,jc,kc,muc), polx%xgrid(:,ec,zc,jc,kc,muc))
	                        polx%apgrid(:,ec,zc,jc,kc,muc) = f_lininterp(policies%xgrid(:,ec,zc,jc,kc,muc), policies%apgrid(:,ec,zc,jc,kc,muc), polx%xgrid(:,ec,zc,jc,kc,muc))
	                        polx%stocks(:,ec,zc,jc,kc,muc) = f_lininterp(policies%xgrid(:,ec,zc,jc,kc,muc), policies%stocks(:,ec,zc,jc,kc,muc), polx%xgrid(:,ec,zc,jc,kc,muc))
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call polx%calc_kappa
    end subroutine InterpolateXgrid
end module interpolate_xgrid
