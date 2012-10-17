module fun_zbrent
    implicit none
contains

pure function f_zbrent(func,x1,x2,tol_o)
! Uses Brent's method to find the root of func that is bracketed by [x1,x2]
! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1188
! Attention: call only if x1, x2 are known to bracket a root, because there is
! no error return if not, zbrent will just return no result (so that the function remains pure).
	use kinds
	real(dp), intent(in) :: x1,x2
	real(dp), intent(in), optional :: tol_o
	real(dp) :: f_zbrent
	interface
		pure function func(x)
		use kinds
		implicit none
		real(dp), intent(in) :: x
		real(dp) :: func
		end function func
	end interface
	integer, parameter :: ITMAX=1000
	real(dp), parameter :: EPS=epsilon(x1)
	integer :: iter
	real(dp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,tol, xm

	if (present(tol_o)) then
	   tol=tol_o
	else
	   tol=(x1+x2)/2.0_dp * 1.0e-8_dp
	endif

	a=x1
	b=x2
	fa=func(a)
	fb=func(b)
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
!		print*, 'f_zbrent: ERROR, initial values dont bracket root'
		return
	endif
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_dp*EPS*abs(b)+0.5_dp*tol
		xm=0.5_dp*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			f_zbrent=b
			return
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_dp*xm*s
				q=1.0_dp-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
				q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_dp*p  <  min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=func(b)
	end do
!	print*, 'f_zbrent: ERROR, exceeded maximum iterations'
	f_zbrent=b
end function f_zbrent
end module fun_zbrent
