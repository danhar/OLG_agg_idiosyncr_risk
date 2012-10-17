module sub_linesearch
    implicit none
contains

subroutine lnsrch(xold,fold,g,p,x,f,fvec,stpmax,check,fname,MaxLns_o)
! Line search: find new x with original Newton step p,scaled down by lambda
! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1195
	use kinds
	real(dp) ,dimension(:)          ,intent(in)    :: xold
	real(dp)                        ,intent(in)    :: fold		! value of distance function at xold
	real(dp) ,dimension(size(xold)) ,intent(in)    :: g		! gradient at xold
	real(dp) ,dimension(size(xold)) ,intent(inout) :: p		! Newton step \delta x
	real(dp) ,dimension(size(xold)) ,intent(out)   :: x
	real(dp) ,dimension(size(xold)) ,intent(inout) :: fvec		! output of fname
    real(dp)                        ,intent(in)    :: stpmax   ! limit length of steps
	real(dp)                        ,intent(out)   :: f		! new value of distance function
	logical                         ,intent(out)   :: check	! false on normal exit, true when x too close to xold. In a minimization algorithm, this usually signals convergence and can be ignored. However, in a zero-finding algorithm the calling program should check whether the convergence is spurious.
    integer  ,optional              ,intent(in)    :: MaxLns_o ! maximum number of line searches
	real(dp) ,parameter	:: ALF=1.0e-4_dp ,tolx=epsilon(x)	! ALF ensures sufficient decrease in function value; tolx is the convergence criterion
	real(dp) 			:: a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,tmplam
	integer				:: n, MaxLns, its
	interface
		function fname(x)
		use kinds
		implicit none
		real(dp), dimension(:), intent(in)	:: x
		real(dp), dimension(size(x))			:: fname
		end function fname
	end interface

    if (present(MaxLns_o)) then
        MaxLns = MaxLns_o
    else
        MaxLns = huge(MaxLns)
    endif

	n=size(xold)
	check=.false.

	! check whether p exeeds maximum step in any direction, normalize all.
	!This differs from NR, who use uses the L2 norm (sqrt(dot_product(p,p)) for check and normalizaton.
	if (maxval(abs(p))>stpmax) then
		p=p*stpmax/maxval(abs(p))
	endif

	slope=dot_product(g,p)
	if (slope >= 0.0) print*, 'roundoff problem in lnsrch'

	alamin=tolx/maxval(abs(p(:))/max(abs(xold(:)),1.0))			! Compute \lambda_{min}
	alam=1.0														! Always try full Newton step first.
	do its = 1,MaxLns
		x(:)=xold(:)+alam*p(:)
		fvec=fname(x)
		f=0.5_dp*dot_product(fvec,fvec)
		if (alam < alamin) then										! Convergence on x
			x(:)=xold(:)
			check=.true.
			return
		elseif (f <= fold+ALF*alam*slope) then						! Sufficient function decrease.
			return
		else														! Backtrack.
			if (alam == 1.0) then									! First time.
				tmplam=-slope/(2.0_dp*(f-fold-slope))
			else													! Subsequent backtracks.
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
				if (a == 0.0) then
					tmplam=-slope/(2.0_dp*b)
				else
					disc=b*b-3.0_dp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_dp*alam
					elseif (b <= 0.0) then
						tmplam=(-b+sqrt(disc))/(3.0_dp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_dp*alam)
	end do

end subroutine lnsrch
end module sub_linesearch
