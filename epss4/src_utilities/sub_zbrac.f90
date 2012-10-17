module sub_zbrac
    implicit none
contains

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
end module sub_zbrac
