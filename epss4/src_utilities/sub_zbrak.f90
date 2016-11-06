module sub_zbrak
    implicit none
contains

subroutine s_zbrak(func,x1,x2,n,xb1,xb2)
! Searches for root of fun by subdividing (x2-x1) into n equi-spaced segments, and returns nb pairs xb1(1:nb), xb2(1,nb)
! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1184
! Attention: made pure, but no error return!
	use kinds
	use numrec_utils, only: arth
	integer, intent(in) :: n
	real(dp), intent(in) :: x1,x2
	real(dp), dimension(:), allocatable, intent(out) :: xb1,xb2
	interface
		function func(x)
		use kinds
		implicit none
		real(dp), intent(in) :: x
		real(dp) :: func
		end function func
	end interface
	integer :: i
	real(dp) :: dx
	real(dp), dimension(0:n) :: f,x
	logical, dimension(1:n) :: mask

	dx=(x2-x1)/n
	x=x1+dx*arth(0.0_dp,1.0_dp,n+1)
	do i=0,n
		f(i)=func(x(i))
	end do
	mask=f(1:n)*f(0:n-1) <= 0.0
	xb1=pack(x(0:n-1),mask)
	xb2=pack(x(1:n),mask)

end subroutine s_zbrak
end module sub_zbrak
