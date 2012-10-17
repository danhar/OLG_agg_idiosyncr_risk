module sub_broyden
    implicit none
contains

subroutine s_broyden(fname,x,fvec,check,df_o,get_fd_jac_o,tolf_o,maxstp_o, maxlnsrch_o)
	! Find root by Broyden's method embedded in a globally convergent strategy.
	! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1199
    ! There are a lot of optional values, but you better know what you are
    ! doing if you deviate from the default values, which are hard-coded in NR
	use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite ! could rename
	use kinds
	use numrec_utils, only: get_diag,lower_triangle, outerprod,put_diag,unit_matrix, qrdcmp,qrupdt,rsolv
	use fun_fdjac
	use sub_linesearch

	real(dp) ,dimension(:) ,intent(inout) :: x			! initial guess, and final root
	real(dp) ,dimension(size(x))          &
						   ,intent(out)	  :: fvec			! values returned by function fname
	logical                ,intent(out)   :: check		! false: converged; true: Broyden's method can make no further progress (e.g. local minimum). Try different initial guess.
	real(dp) ,dimension(size(x),size(x))   &
			 ,optional     ,intent(inout) :: df_o			! Initial Jacobian (in), last jacobian (out)
	logical  ,optional     ,intent(in)    :: get_fd_jac_o	! true: compute Jacobian by first difference (fd)
	real(dp) ,optional     ,intent(in)    :: tolf_o		! convergence criterion on function values
	real(dp) ,optional     ,intent(in)    :: maxstp_o		! maximum step in line search algorithm <1
	integer  ,optional     ,intent(in)    :: maxlnsrch_o  ! maximum number of line searches
	real(dp) ,dimension(size(x),size(x))  :: df			! stands in for the optional df_o
	logical								  :: get_fd_jac	! stands in for the optional get_fd_jac_o
	real(dp)							  :: tolf			! stands in for the optional tolf_o
	real(dp)							  :: maxstp		! stands in for the optional maxstp_o
	integer  ,parameter					  :: max_it=200	! maximum number of iterations
	real(dp) ,parameter 				  :: eps=epsilon(x)
	real(dp) ,parameter                   :: tolmin=1.0e-12_dp,tolx=eps !criterion for convergence to a minimum of fmin; convergence criterion on x
	real(dp) 							  :: f,fold		! f: value of distance function; fold: value of distance function at xold
	integer 							  :: n,i,its,k
	real(dp) ,dimension(size(x)) 		  :: c,d,fvcold,g,p,s,t,w,xold
	real(dp) ,dimension(size(x),size(x))  :: r,qt
	logical 					          :: restart,sing, initial_df_present
	interface
		function fname(x)
			use kinds, only: dp
			implicit none
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(size(x))	   :: fname
		end function fname
	end interface

    n= size(x)
	if (present(maxstp_o)) then
		maxstp = maxstp_o
	else
		maxstp = 100.0_dp*max(sqrt(dot_product(x,x)),real(n,dp)) ! NR value
	endif
	if (present(tolf_o)) then
		tolf = tolf_o
	else
		tolf = 1e-8_dp    ! Value proposed by NR, p.1363, for double precision
	endif
	if (present(get_fd_jac_o)) then
		get_fd_jac=get_fd_jac_o
	elseif (.not.present(df_o)) then
		get_fd_jac=.true.		! if df_o is not supplied, make sure jacobian will be evaluated by first differences
	else
		get_fd_jac=.false.
	endif
	if (present(df_o)) then
		initial_df_present=.true.
	else
		initial_df_present=.false.
	endif

	restart=.true.										! so that Jacobian is computed
 	fvec=fname(x)
	f=0.5_dp*(dot_product(fvec,fvec))


    if (any(ieee_is_nan(fvec)) .or. any(.not. ieee_is_finite(fvec))) then
        print *, 'sub_alg_qn: WARNING NaN or Inf encountered'
        check =.true.
        return
    elseif (maxval(abs(fvec(:))) < tolf) then ! Test for initial guess being a root.
        check=.false.
        return
    endif

	do its=1,max_it												! Start of iteration loop.
		if (restart) then
			if (get_fd_jac .and..not. initial_df_present) then
				print*,'sub_broyden: computing finite difference Jacobian'
				df=f_fdjac(fname,x,fvec)
			else
				df=df_o
				initial_df_present=.false.		! this ensures that we enter the first branch next time if get_fd_jac is true
			endif
			r=df												! Initialize or reinitialize Jacobian in r.
			call qrdcmp(r,c,d,sing)								! QR decomposition of Jacobian.
			if (sing) then
				print*,'sub_broyden: singular Jacobian'
				return
			endif
			call unit_matrix(qt)								! Form QT explicitly.
			do k=1,n-1
				if (c(k) /= 0.0) then
					qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),matmul(r(k:n,k),qt(k:n,:)))/c(k)
				endif
			end do
			where (lower_triangle(n,n)) r(:,:)=0.0_dp
			call put_diag(d(:),r(:,:))							! Form R explicitly.
		else													! Carry out Broyden update.
			s(:)=x(:)-xold(:)
			do i=1,n
				t(i)=dot_product(r(i,i:n),s(i:n))
			end do
			w(:)=fvec(:)-fvcold(:)-matmul(t(:),qt(:,:))
			where (abs(w(:)) < eps*(abs(fvec(:))+abs(fvcold(:)))) w(:)=0.0_dp	! Don't update with noisy components of w
			if (any(w(:) /= 0.0)) then
				t(:)=matmul(qt(:,:),w(:))
				s(:)=s(:)/dot_product(s,s)
				call qrupdt(r,qt,t,s)							! Update R and QT .
				d(:)=get_diag(r(:,:))							! Diagonal of R stored in d.
				if (any(d(:) == 0.0)) print*, 'sub_broyden: r singular'
			endif
		endif
		p(:)=-matmul(qt(:,:),fvec(:))							! r.h.s. for linear equations
		do i=1,n												! Compute gradient for the line search
			g(i)=-dot_product(r(1:i,i),p(1:i))
		end do
		xold(:)=x(:)											! Store x, F, and f.
		fvcold(:)=fvec(:)
		fold=f
		call rsolv(r,d,p)										! Solve linear equations for step

		! lnsrch returns new x and f. It also calculates fvec at the new x
		call lnsrch(xold,fold,g,p,x,f,fvec,maxstp,check,fname,maxlnsrch_o)

        if (any(ieee_is_nan(fvec)) .or. any(.not. ieee_is_finite(fvec))) then
            print *, 'sub_alg_qn: WARNING NaN or Inf encountered'
            check =.true.
            return
        elseif (maxval(abs(fvec(:))) < tolf) then         ! Test for convergence on function values.
            check=.false.
            if (present(df_o)) df_o = df
            return
        endif

		if (check) then											! True if line search failed to find a new x
			! If restart is true we have failure: We have already tried reinitializing the Jacobian.
			! The other test is for gradient of f= zero, i.e., spurious convergence.
			if (restart .or. maxval(abs(g(:))*max(abs(x(:)),1.0)/max(f,0.5_dp*n)) < tolmin) then
				print*,'sub_broyden: cannot converge, try new starting point'
				if (present(df_o)) df_o = df
				return
			endif
			restart=.true. 										! Try reinitializing the Jacobian.
		else													! Successful step; will use Broyden update
			restart=.false.
			if ( maxval((abs(x(:)-xold(:)))/max(abs(x(:)),1.0)) < tolx) then			! Test for convergence on dx
				if (present(df_o)) df_o = df
!				check=.true.
				return
			endif
		endif
	end do
	check=.true.
	print*,'sub_broyden: max_it exceeded'
	if (present(df_o)) df_o = df

end subroutine s_broyden
end module sub_broyden
