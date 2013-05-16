module sub_alg_qn
    implicit none
contains

subroutine s_alg_qn(fname,fvec,x,n,QTmat,Rmat,intlj,reevalj,check,rstit0,MaxLns,max_it,maxstp,tol_f)

    use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite ! could rename
	use kinds
    use numrec_utils, only: get_diag,lower_triangle, outerprod,put_diag,unit_matrix, qrdcmp,qrupdt,rsolv
    use fun_fdjac
    use sub_linesearch

    real(dp) ,dimension(n)   ,intent(out)   :: fvec    ! output of fname
    real(dp) ,dimension(n)   ,intent(inout) :: x       ! input of fname
	integer                  ,intent(in) 	:: n, rstit0 	  ! size of x andfvec; number of iterations before Jacobian is re-initialized
	integer                  ,intent(in) 	:: MaxLns, max_it ! max no of line searches, max no of iterations
	logical                  ,intent(in)	:: reevalj ! true: reevaulate Jacobian if line-search fails
	real(dp)                 ,intent(in) 	:: tol_f, maxstp  ! tolerance for f(x); maxstp is passed to lnsrch (maximum Newton step)
	logical                  ,intent(inout)	:: intlj   ! true: Jacobian and Q-R-decomp are computed in iteration 1, false: input Q and R
    real(dp) ,dimension(n,n) ,intent(out)   :: QTmat   ! Q-R-decomposition of Jacobian
	real(dp) ,dimension(n,n) ,intent(inout)	:: Rmat	   ! Q-R-decomposition of Jacobian
	logical                  ,intent(out)   :: check   ! true if nonconvergent
	integer 								:: i,its,k,rstit
	real(dp) 								:: tol_x,f,fold
	real(dp) ,dimension(n) 				    :: c,d,fvcold,g,p,s,t,w,xold
	real(dp) ,dimension(n,n) 				:: rold,qtold
	logical 								:: restrt,sing
	real(dp)                 ,parameter     :: eps=epsilon(x)

    interface
        function fname(x)
            use kinds, only: dp
            implicit none
            real(dp), dimension(:), intent(in) :: x
            real(dp), dimension(size(x))       :: fname
        end function fname
    end interface
	
	tol_x=tol_f			! tolerance on X should not be lower than tolalg
	restrt=.true.

    fvec=fname(x)
    f=0.5_dp*(dot_product(fvec,fvec))    
	print '(a,es13.6)', ' sub_alg_qn: maxval(abs(fvec(:))) in iteration    1: ', maxval(abs(fvec(:)))
	print '(a,es13.6)', '             0.5*( fvec*fvec )                     : ', f

    if (any(ieee_is_nan(fvec)) .or. any(.not. ieee_is_finite(fvec))) then
        print *, 'sub_alg_qn: WARNING NaN or Inf encountered'
        check =.true.
        return
	elseif (maxval(abs(fvec(:))) < tol_f) then ! Test for initial guess being a root.
		check=.false.
		return
	endif
	
	! vabs=sqrt(dot_product(x(:),x(:)))
	! maxstp=STPMX*max(vabs,1.0*n)						! Calculate maxstp for line searches.

	rold=Rmat
	qtold=QTmat
	rstit=rstit0
	do its=1,max_it										! Start of iteration loop
		if (its==rstit) then
			rstit=rstit+rstit0
			restrt=.true.
		endif

		if (restrt) then

			if (intlj) then								! initialj=1 if initial Jacobi evaluation is required
				! initialize Jacobi
				Rmat=f_fdjac(fname,x,fvec)
				print*, 'sub_alg_qn: done with Jacobi initialization'
			elseif (reevalj) then
				intlj=.true.							! make sure to in next iter if still required
			else
				Rmat=rold								! make sure to always reinitialize with old Jacobi matrix
				QTmat=qtold
			endif

			call qrdcmp(Rmat,c,d,sing)										! QR decomposition of Jacobian.
			if (sing) then
				print*,'sub_alg_qn: singular Jacobian in broyden'
				return
			endif
			call unit_matrix(QTmat)										! Form QT explicitly.
			do k=1,n-1
				if (c(k) /= 0.0) then
					QTmat(k:n,:)=QTmat(k:n,:)-outerprod(Rmat(k:n,k),&
					matmul(Rmat(k:n,k),QTmat(k:n,:)))/c(k)
				endif
			end do
			where (lower_triangle(n,n)) Rmat(:,:)=0.0
			call put_diag(d(:),Rmat(:,:))									! Form R explicitly.

		else															! Carry out Broyden update.

			s(:)=x(:)-xold(:)
			do i=1,n
				t(i)=dot_product(Rmat(i,i:n),s(i:n))
			end do
			w(:)=fvec(:)-fvcold(:)-matmul(t(:),QTmat(:,:))

			! it is better to update in any case even though the following conditions do not hold
			! where (abs(w(:)) < eps*(abs(fvec(:))+abs(fvcold(:)))) &
			! 	w(:)=0.0												! Dont update with noisy components of w

			if (any(w(:) /= 0.0)) then
				t(:)=matmul(QTmat(:,:),w(:))
				s(:)=s(:)/dot_product(s,s)
				call qrupdt(Rmat,QTmat,t,s)									! Update R and QT .
				d(:)=get_diag(Rmat(:,:))									! Diagonal of R stored in d.
				if (any(d(:) == 0.0)) &
					print*, 'sub_alg_qn: Rmat singular in broyden'
			endif
		endif
		p(:)=-matmul(QTmat(:,:),fvec(:))									! r.h.s. for linear equations
		do i=1,n														! Compute gradient for the line search
			g(i)=-dot_product(Rmat(1:i,i),p(1:i))
		end do

		xold(:)=x(:)													! Store x, F, and f.
		fvcold(:)=fvec(:)
		fold=f
		call rsolv(Rmat,d,p)												! Solve linear equations for step

		! lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.
		call lnsrch(xold,fold,g,p,x,f,fvec,maxstp,check,fname,MaxLns)

		print '(a,i4,a2,es13.6)', ' sub_alg_qn: maxval(abs(fvec(:))) in iteration ',its+1, ': ', maxval(abs(fvec(:)))
        print '(a,es13.6)', '             0.5*( fvec*fvec )                     : ', f

	    if (any(ieee_is_nan(fvec)) .or. any(.not. ieee_is_finite(fvec))) then
            print *, 'sub_alg_qn: WARNING NaN or Inf encountered'
	        check =.true.
	        return
	    elseif (maxval(abs(fvec(:))) < tol_f) then         ! Test for convergence on function values.
			check=.false.
			return
		endif
		if (.not.(check)) then											! True if line search failed to find a new x
			restrt=.false.
			if ( maxval(abs(x(:)-xold(:))) < tol_x) then					! Test for convergence on dx: Here: absolute deviation
				! better not do anything here!
				! restrt=.true.
			endif
		else
			if (.not.(restrt)) then
				restrt=.true.       ! Try re-initializing Jacobian
			! ELSE: keep your fingers crossed and continue
			endif
		endif
	end do
	check=.true.
    print*,'sub_alg_qn: ERROR: max_it exceeded'

end subroutine s_alg_qn
end module sub_alg_qn
