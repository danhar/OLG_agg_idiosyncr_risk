module fun_fdjac
    implicit none
contains

function f_fdjac(fname,x,fvec)
! Computes forward-difference approximation to Jacobian.
! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1197
	use kinds

	real(dp) ,dimension(:)       ,intent(in) :: x		! point at which Jacobian is to be evaluated
	real(dp) ,dimension(size(x)) ,intent(in) :: fvec	! vector of function values at x
	real(dp) ,dimension(size(x),size(x)) 	 :: f_fdjac	! Jacobian
	real(dp) ,parameter 					 :: eps=1.0e-4_dp
	! eps is different from NR! NR: approximate square root of the machine precision: 1.0e-8_dp.
	! But we need larger value because we have numerical approximation of function value.
	integer 								 :: j,n
	real(dp), dimension(size(x)) 			 :: xtemp,xph,h
	interface
		function fname(x)
		use kinds
		implicit none
		real(dp), dimension(:), intent(in)	 :: x
		real(dp), dimension(size(x))		 :: fname
		end function fname
	end interface

	n		= size(x)
	xtemp	= x
	h		= eps*abs(xtemp)
	where (h == 0.0) h=eps
	xph		= xtemp+h							    ! Trick to reduce finite precision error.
	h		= xph-xtemp
	do j=1,n
		xtemp(j)	= xph(j)
		f_fdjac(:,j)= (fname(xtemp)-fvec(:))/h(j)	! Forward difference formula.
		xtemp(j)	= x(j)
	end do

end function f_fdjac
end module fun_fdjac
