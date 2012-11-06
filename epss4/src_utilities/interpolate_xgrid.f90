module interpolate_xgrid
! This interpolation is specific, because it first creates the new xgrid, then acts on
! various objects. More elegant would be to have a function create xgrid, and have
! interpolation functions for each object.
! This module is USEd only in krusell_smith, the function InterpolateXgrid is only executed IF (exogenous_xgrid).
    implicit none
    private
    public InterpolateXgrid

    interface InterpolateXgrid
        module procedure InterpolateXgrid_policies, InterpolateXgrid_Phi
    end interface InterpolateXgrid

contains
    pure subroutine InterpolateXgrid_policies(nx_factor,policies, value, polx, valx)
        use kinds      ,only: dp
        use policies_class
        use makegrid_mod
        use fun_lininterp

        integer         ,intent(in)    :: nx_factor
        type(tPolicies) ,intent(in)    :: policies
        real(dp)        ,intent(in)    :: value(:,:,:,:,:,:)
        type(tPolicies) ,intent(out)   :: polx
        real(dp) ,allocatable ,intent(out)   :: valx(:,:,:,:,:,:)
        integer                        :: muc, kc, jc, zc, ec, nx, nz, nk, nmu, n_eta, nj
        real(dp)                       :: xmin, xmax
        real(dp) ,parameter            :: cmin = 1e-10

        nx = size(value,1); n_eta=size(value,2); nz = size(value,3); nj=size(value,4); nk = size(value,5); nmu=size(value,6)
        call polx%allocate(nx*nx_factor, nz, nk,nmu)
        allocate(valx(nx*nx_factor,n_eta,nz,nj,nk,nmu))
        do muc = 1, nmu
            do kc = 1, nk
                do jc=1,nj
                    do zc = 1,nz
                        do ec = 1,n_eta
                            ! mean of first and last mu seems good, coz small aggregate grid and policies linear in that dimension
                            xmin = 0.5_dp * policies%xgrid(1,ec,zc,jc,kc,1) + 0.5_dp * policies%xgrid(1,ec,zc,jc,kc,nmu)
                            xmax = 0.5_dp * policies%xgrid(nx,ec,zc,jc,kc,1) + 0.5_dp * policies%xgrid(nx,ec,zc,jc,kc,nmu)
                            polx%xgrid(:,ec,zc,jc,kc,muc) = MakeGrid(xmin,xmax,nx*nx_factor, 2.0_dp)

	                        valx(:,ec,zc,jc,kc,muc)        = f_lininterp(policies%xgrid(:,ec,zc,jc,kc,muc), value(:,ec,zc,jc,kc,muc), polx%xgrid(:,ec,zc,jc,kc,muc))
	                        polx%apgrid(:,ec,zc,jc,kc,muc) = f_lininterp(policies%xgrid(:,ec,zc,jc,kc,muc), policies%apgrid(:,ec,zc,jc,kc,muc), polx%xgrid(:,ec,zc,jc,kc,muc))
	                        polx%stocks(:,ec,zc,jc,kc,muc) = f_lininterp(policies%xgrid(:,ec,zc,jc,kc,muc), policies%stocks(:,ec,zc,jc,kc,muc), polx%xgrid(:,ec,zc,jc,kc,muc))

	                        where ((polx%xgrid(:,ec,zc,jc,kc,muc) - polx%apgrid(:,ec,zc,jc,kc,muc)) < cmin) polx%apgrid(:,ec,zc,jc,kc,muc) = polx%xgrid(:,ec,zc,jc,kc,muc) - cmin

                        enddo
                    enddo
                enddo
            enddo
        enddo
        call polx%calc_kappa

    end subroutine InterpolateXgrid_policies

    pure subroutine InterpolateXgrid_Phi(Phi, xgrid_old, xgrid_new)
        use kinds      ,only: dp
        use fun_lininterp

        real(dp) ,intent(inout) :: Phi(:,:,:)
        real(dp) ,intent(in)    :: xgrid_old(:,:,:), xgrid_new(:,:,:)
        integer                 :: jc, ec

        do jc=1,size(Phi,3)
            do ec=1,size(Phi,2)
                Phi(:,ec,jc) = f_lininterp(xgrid_old(:,ec,jc), Phi(:,ec,jc), xgrid_new(:,ec,jc))
            enddo
        enddo
    end subroutine InterpolateXgrid_Phi

end module interpolate_xgrid
