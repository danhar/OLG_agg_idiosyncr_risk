module fun_lininterp
    use kinds
    implicit none
    private
    public f_lininterp

    interface f_lininterp
        module procedure lininterp_vector, lininterp_scalar
    end interface

contains

pure function lininterp_vector(xgrid,fvals,x)
! linear interpolation
! compare to su_myinterp1.f90 of villaverde, maybe change!

    use fun_locate
    real(dp) ,dimension(:)       ,intent(in) :: xgrid,fvals   ! x-values and function values at those x-values
    real(dp) ,dimension(:)       ,intent(in) :: x             ! x values to interpolate
    real(dp) ,dimension(size(x))             :: lininterp_vector, w
    integer  ,dimension(size(x))             :: i

	i=f_locate(xgrid,x)		! In 'default', returns ju-1 if x>xgrid(ju-1) !

	! Note that this also obtains for linear extrapolation at upper and lower bounds
	w=(x-xgrid(i))/(xgrid(i+1)-xgrid(i))
	lininterp_vector=(1-w)*fvals(i)+w*fvals(i+1)

end function lininterp_vector

pure function lininterp_scalar(xgrid,fvals,x)
    real(dp)                           :: lininterp_scalar
    real(dp) ,dimension(:) ,intent(in) :: xgrid,fvals   ! x-values and function values at those x-values
    real(dp)               ,intent(in) :: x             ! scalar
    real(dp) ,dimension(1)             :: x_array, interp_array

    x_array(1)   = x
    interp_array = lininterp_vector(xgrid,fvals,x_array)
    lininterp_scalar = interp_array(1)

end function lininterp_scalar

end module fun_lininterp
