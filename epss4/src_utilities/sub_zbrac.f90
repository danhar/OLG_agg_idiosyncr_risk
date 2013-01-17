module sub_zbrac
    implicit none
    private
    public s_zbrac, s_zbrac_array

contains
!---------------------------------------------------------------------------
! Module procedures in order:
! - pure subroutine s_zbrac(func,x1,x2,success)
! - subroutine s_zbrac_array(func,x1,x2,success)
!---------------------------------------------------------------------------

pure subroutine s_zbrac(func,x1,x2,success)
! Attempts to bracket a root of func by expanding [x1,x2] by FACTOR=1.6
! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1183
! Attention: made pure, but no error return!
	use kinds
	real(dp), intent(inout) :: x1,x2
	logical, intent(out) :: success
	interface
        pure function func(x)
		use kinds
		implicit none
		real(dp), intent(in) :: x
		real(dp) :: func
		end function func
	end interface
	integer, parameter :: NTRY=50
	real(dp), parameter :: FACTOR=1.6_dp
	integer :: j
	real(dp) :: f1,f2
	if (x1 == x2) then
!	   print*,'sub_zbrac: ERROR, you have to guess an initial range'
	   success=.false.
	   return
	endif
	f1=func(x1)
	f2=func(x2)
	success=.true.
	do j=1,NTRY
		if ((f1 > 0.0 .and. f2 < 0.0) .or. (f1 < 0.0 .and. f2 > 0.0)) return
		if (abs(f1) < abs(f2)) then
			x1=x1+FACTOR*(x1-x2)
			f1=func(x1)
		else
			x2=x2+FACTOR*(x2-x1)
			f2=func(x2)
		end if
	end do
	success=.false.
end subroutine s_zbrac
!---------------------------------------------------------------------------

subroutine s_zbrac_array(func,x1,x2,success)
! Same as above, only that it func takes an array with a scalar element
! and is not pure. And changed to FACTOR=1.0_dp !
    use kinds
    real(dp), intent(inout) :: x1,x2
    logical, intent(out) :: success
    interface
        function func(x)
        use kinds
        implicit none
        real(dp), intent(in) :: x(:)
        real(dp) :: func(size(x))
        end function func
    end interface
    integer, parameter :: NTRY=50
    real(dp), parameter :: FACTOR=1.0_dp  ! ATTN: orig was FACTOR=1.6_dp
    integer :: j
    real(dp) :: f1,f2
    if (x1 == x2) then
       print*,'sub_zbrac: ERROR, you have to guess an initial range'
       success=.false.
       return
    endif
    f1=sum(func([x1])) ! sum only a hack to make it a scalar
    f2=sum(func([x2]))
    success=.true.
    do j=1,NTRY
        if ((f1 > 0.0 .and. f2 < 0.0) .or. (f1 < 0.0 .and. f2 > 0.0)) return
        if (abs(f1) < abs(f2)) then
            x1=x1+FACTOR*(x1-x2)
            f1=sum(func([x1])) ! sum only a hack to make it a scalar
        else
            x2=x2+FACTOR*(x2-x1)
            f2=sum(func([x2]))
        end if
    end do
    success=.false.
end subroutine s_zbrac_array

end module sub_zbrac
