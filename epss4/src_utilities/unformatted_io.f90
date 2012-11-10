module unformatted_io
    use kinds           ,only: dp
    use types           ,only: tSimvars, AllocateType
    use aggregate_grids ,only: tAggGrids, AllocateType
    use laws_of_motion  ,only: tCoeffs
    implicit none
    private

    public ReadUnformatted, SaveUnformatted
    interface ReadUnformatted
        module procedure ReadUnformatted_msgrids, ReadUnformatted_msvars, ReadUnformatted_ksvars
    end interface ReadUnformatted

    interface SaveUnformatted
        module procedure SaveUnformatted_msvars, SaveUnformatted_ksvars
    end interface SaveUnformatted

contains

    !---------------------------------------------------------------------------
    ! Save/ load unformatted MS variables
    !---------------------------------------------------------------------------
    subroutine ReadUnformatted_msgrids(ms_grids)
        type(tAggGrids), intent(out) :: ms_grids
        call AllocateType(ms_grids,1,1)
        open(55,file='model_input/last_results/grids_ms.unformatted'  ,form='unformatted',action='read')
        read(55) ms_grids%k, ms_grids%mu
        close(55)
    end subroutine ReadUnformatted_msgrids

    subroutine ReadUnformatted_msvars(Phi, ms_grids)
    ! Not used so far!
        real(dp), dimension(:,:,:), intent(inout) :: Phi
        type(tAggGrids), intent(out) :: ms_grids
        call AllocateType(ms_grids,1,1)
        open(55,file='model_input/last_results/Phi_ms.unformatted',form='unformatted',action='read')
        read(55) Phi
        close(55)
        open(55,file='model_input/last_results/grids_ms.unformatted'  ,form='unformatted',action='read')
        read(55) ms_grids%k, ms_grids%mu
        close(55)
!        open(55,file='model_input/last_results/def_benefits.unformatted',form='unformatted',action='read')
!        read(55) pensions
!        close(55)
    end subroutine ReadUnformatted_msvars

    subroutine SaveUnformatted_msvars(Phi, ms_grids, pensions)
    ! Writes Phi_ms, K_ms, mu_ms unformatted to model_input! For later use in PE / GE.
        real(dp)        ,intent(in) :: Phi(:,:,:), pensions
        type(tAggGrids) ,intent(in) :: ms_grids
        open(55,file='model_input/last_results/Phi_ms.unformatted',form='unformatted')
        write(55) Phi
        close(55)
        open(55,file='model_input/last_results/grids_ms.unformatted',form='unformatted')
        write(55) ms_grids%k, ms_grids%mu
        close(55)
        open(55,file='model_input/last_results/def_benefits.unformatted',form='unformatted')
        write(55) pensions
        close(55)
    end subroutine SaveUnformatted_msvars

    !---------------------------------------------------------------------------
    ! Save /load ksvars unformatted for PE
    !---------------------------------------------------------------------------
    subroutine SaveUnformatted_ksvars(grids, coeffs, simvars)
        type(tAggGrids) ,intent(in) :: grids
        type(tCoeffs)   ,intent(in) :: coeffs
        type(tSimvars)  ,intent(in) :: simvars(:)
        integer :: i

        open(55,file='model_input/last_results/grids_ge.unformatted',form='unformatted')
        write(55) grids%k, grids%mu
        close(55)

        open(55,file='model_input/last_results/coeffs_ge.unformatted',form='unformatted')
        write(55) coeffs%k, coeffs%mu
        close(55)

        open(55,file='model_input/last_results/simvars_ge.unformatted',form='unformatted')
        ! Could I make this elemental instead of the do-loop? or use standard derived type IO
        do i=1,size(simvars)
            write(55) simvars(i)%z, simvars(i)%K, simvars(i)%mu, simvars(i)%B, simvars(i)%C, simvars(i)%Phi_1, simvars(i)%Phi_nx, simvars(i)%err_aggr, &
                      simvars(i)%err_income, simvars(i)%r, simvars(i)%rf, simvars(i)%wage, simvars(i)%pens, simvars(i)%tau, simvars(i)%welf, simvars(i)%bequests, simvars(i)%err_K, simvars(i)%err_mu
        enddo
        close(55)

        open(55,file='model_input/last_results/nt.unformatted',form='unformatted')
        write(55) size(simvars(1)%z)
        close(55)

        open(55,file='model_input/last_results/size_simvars.unformatted',form='unformatted')
        write(55) size(simvars)
        close(55)
    end subroutine SaveUnformatted_ksvars

    subroutine ReadUnformatted_ksvars(grids, coeffs, simvars)
        use params_mod  ,only: nk, nmu, n_coeffs, nz
        type(tAggGrids) ,intent(out) :: grids
        type(tCoeffs)   ,intent(out) :: coeffs
        type(tSimvars)  ,allocatable ,intent(out) :: simvars(:)
        integer :: i

        call AllocateType(grids,nk,nmu)
        open(55,file='model_input/last_results/grids_ge.unformatted',form='unformatted',action='read')
        read(55) grids%k, grids%mu
        close(55)

        allocate(coeffs%k(n_coeffs,nz), coeffs%mu(n_coeffs,nz), coeffs%r_squared(2,nz))
        open(55,file='model_input/last_results/coeffs_ge.unformatted',form='unformatted',action='read')
        read(55) coeffs%k, coeffs%mu
        close(55)

        open(55,file='model_input/last_results/size_simvars.unformatted',form='unformatted')
        read(55) i
        close(55)
        allocate(simvars(i))

        open(55,file='model_input/last_results/nt.unformatted',form='unformatted')
        read(55) i
        close(55)
        call AllocateType(simvars,i)

        open(55,file='model_input/last_results/simvars_ge.unformatted',form='unformatted',action='read')
        do i=1,size(simvars)
            read(55) simvars(i)%z, simvars(i)%K, simvars(i)%mu, simvars(i)%B, simvars(i)%C, simvars(i)%Phi_1, simvars(i)%Phi_nx, simvars(i)%err_aggr, &
                     simvars(i)%err_income, simvars(i)%r, simvars(i)%rf, simvars(i)%wage, simvars(i)%pens, simvars(i)%tau, simvars(i)%welf, simvars(i)%bequests, simvars(i)%err_K, simvars(i)%err_mu
        enddo
        close(55)
    end subroutine ReadUnformatted_ksvars

end module unformatted_io
