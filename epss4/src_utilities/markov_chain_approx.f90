module markov_chain_approx
! This module is adapted from the toolbox by Iskander Karibzhanov
    use kinds
    implicit none
contains

	subroutine tauchen(rho,sigma,p,y,s,tauch_spread)
	! Tauchen method to approximate univariate AR(1) process by Markov chain
	!      y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
	! INPUTS: rho - serial correlation coefficient,
	!         sigma - coefficient of variation
	! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
	!         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
	!         s is an n-by-1 vector of stationary distribution
	    real(dp), intent(in) :: rho, sigma
	    real(dp), dimension(:,:), intent(out) :: p
	    real(dp), dimension(:), intent(out) :: y,s
	    real(dp), optional, intent(in) :: tauch_spread
	    integer :: n
	    real(dp) :: smin,emin
	    n=size(y)
	    if (size(p,dim=1)/=n .or. size(p,dim=2)/=n) then
	        print '(a,i3,a,i3)', 'tauchen: p must be a square matrix of size ',n,' x ',n
	        stop 'program terminated by tauchen'
	    end if
	    if (size(s)/=n) then
	        print '(a,i3)', 'tauchen: y and s must be vectors of the same size ',n
	        stop 'program terminated by tauchen'
	    end if
	    if (present(tauch_spread)) then
	        call tauch(tauch_spread)
	    else
	        emin=brent(1.0_dp,2.5_dp,4.0_dp,err,smin)
	!       PRINT '(a,f,4x,a,e)', 'smin=',smin,'emin=',emin
	        call tauch(smin)
	    end if
	contains
	subroutine tauch(m)
	    real(dp), intent(in) :: m
	    real(dp) :: ybar, dy, sr
	    real(dp), dimension(size(y)) :: yd
	    integer :: i
	    ybar=m*sigma
	    call grid(y,-ybar,ybar,1.0_dp)
	    dy=ybar/(n-1)
	    sr=sigma*sqrt(1-rho**2)
	    do i=1,n
	        yd=(y-rho*y(i)+dy)/sr
	        ! vdcdfnorm is a MKL procedure. vector double cdf normal distribution
	        call vdcdfnorm(n,yd,p(i,:))
	    end do
	    p(:,n)=1
	    do i=n,2,-1
	        p(:,i)=p(:,i)-p(:,i-1)
	    end do
	    call ergodic(p,s)
	end subroutine tauch
	real(dp) function err(m)
        real(dp), intent(in) :: m
        real(dp) :: rho_, sigma_

        call tauch(m)
        call markovtest(p,y,s,rho_,sigma_)
        err=(log(1-rho_)/log(1-rho)-1)**2+(sigma_/sigma-1)**2;
    end function err
	end subroutine tauchen

	subroutine markovtest(p,y,s,rho,sigma)
	    real(dp), dimension(:,:), intent(in) :: p
	    real(dp), dimension(:), intent(in) :: y,s
	    real(dp), intent(out) :: rho, sigma
	    real(dp), dimension(size(y)) :: py
	    real(dp) :: Eyy, Ey2, E2y
	    integer :: n
	    n=size(y)
	    if (size(p,dim=1)/=n .or. size(p,dim=2)/=n) then
	        print '(a,i3,a,i3)', 'markovtest: p must be a square matrix of size ',n,' x ',n
	        stop 'program terminated by markovtest'
	    end if
	    if (size(s)/=n) then
	        print '(a,i3)', 'markovtest: y and s must be vectors of the same size ',n
	        stop 'program terminated by markovtest'
	    end if
	    py = matmul(p,y)    ! DH: correct?
	    Eyy=sum(s*y*py)
	    Ey2=sum(s*y**2)
	    E2y=sum(s*y)**2;
	    rho=(Eyy-E2y)/(Ey2-E2y);
	    sigma=sqrt(Ey2-E2y);
	end subroutine markovtest

	subroutine rouwenhorst(rho,sigma_eps,p,y,s)
	! Rouwenhorst method to approximate univariate AR(1) process by Markov chain
	!      y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
	! INPUTS: rho - serial correlation coefficient,
	!         sigma - standard deviation of innovation epsilon
	! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
	!         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
	!         s is an n-by-1 vector of stationary distribution (binomial)
	    real(dp), intent(in) :: rho, sigma_eps
	    real(dp), dimension(:,:), intent(out) :: p
	    real(dp), dimension(:), intent(out) :: y,s
	    real(dp) :: ybar, q, sigma_z
	    integer :: n
	    n=size(y)
	    if (size(p,dim=1)/=n .or. size(p,dim=2)/=n) then
	        print '(a,i3,a,i3)', 'rouwenhorst: p must be a square matrix of size ',n,' x ',n
	        stop 'program terminated by rouwenhorst'
	    end if
	    if (size(s)/=n) then
	        print '(a,i3)', 'rouwenhorst: y and s must be vectors of the same size ',n
	        stop 'program terminated by rouwenhorst'
	    end if
	    sigma_z = sigma_eps/sqrt((1.0-rho**2))
	    ybar=sigma_z*sqrt(real(n-1,dp))
	    q=(1+rho)/2
	    call rhmat(p)
	    call grid(y,-ybar,ybar,1.0_dp)
	    call binom(s)
	contains
	recursive subroutine rhmat(p)
	    real(dp), dimension(:,:), intent(out) :: p
	    real(dp), dimension(size(p,dim=1)-1,size(p,dim=2)-1) :: p1
	    integer :: h
	    h=size(p,dim=1)
	    if (size(p,dim=2)/=h) stop 'P must be a square matrix'
	    if (h<2) stop 'P must be at least 2-by-2 matrix'
	    if (h==2) then
	        p=reshape((/q,1-q,1-q,q/),(/2,2/))
	    else
	        call rhmat(p1)
	        p=0.0
	        p(1:h-1,1:h-1)=q*p1
	        p(1:h-1,2:h)=(1-q)*p1+p(1:h-1,2:h)
	        p(2:h,1:h-1)=(1-q)*p1+p(2:h,1:h-1)
	        p(2:h,2:h)=q*p1+p(2:h,2:h)
	        p(2:h-1,:)=p(2:h-1,:)/2.0
	    end if
	end subroutine rhmat
	end subroutine rouwenhorst

subroutine ergodic(p,s)
! Purpose: Compute ergodic distribution s of Markov transition matrix p
    use mkl95_lapack, only: geev
    implicit none
    real(dp), dimension(:,:), intent(in) :: p
    real(dp), dimension(:), intent(out) :: s
    real(dp), dimension(size(s),size(s)) :: ip,vl
    real(dp), dimension(size(s)) :: wr,wi,dw
    integer :: m,info,uw,w1(1)
    real(dp) :: ds

    m=size(s)
    if (size(p,dim=1)/=m .or. size(p,dim=2)/=m) then
        print '(a,i3,a,i3)', 'sd: p must be a square matrix of size ',m,' x ',m
        stop 'program terminated by sd'
    end if
    ip=p
    call geev(ip,wr,wi,vl=vl,info=info)
    call check('geev',info)
    dw=abs(sqrt(wr*wr+wi*wi)-1)
    w1=minloc(dw)
    uw=count(dw<1000*epsilon(dw))
    if (uw<1) print '(a)', 'Warning: No unitary eigenvalue is found. Stationary distribution of Markov chain does not exist.'
    if (uw>1) print '(a)', 'Warning: More than one unitary eigenvalue is found. Stationary distribution of Markov chain is not unique.'
    if (uw<1 .or. uw>1) print '(a,f20.15,a,f20.15)', 'Using eigenvalue ',wr(w1(1)),'+i',wi(w1(1))
    s=vl(:,w1(1))/sum(vl(:,w1(1)))
    if (any(s<0)) then
        print '(a)', 'The stationary distribution of Markov chain has negative values. Rebalancing...'
        ds=sum(s,mask=s<0)/count(s>=0)
        where(s<0)
            s=0
        elsewhere
            s=s+ds
        end where
    end if
end subroutine ergodic

subroutine grid(x,xmin,xmax,s)
! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
! s=1       linear spacing
! s>1       left skewed grid spacing with power s
! 0<s<1     right skewed grid spacing with power s
! s<0       geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
! s=-1      logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
! s=0       logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
    implicit none
    real(dp), dimension(:), intent(out) :: x
    real(dp), intent(in) :: xmin,xmax,s
    real(dp) :: c ! growth rate of grid subintervals for logarithmic spacing
    integer :: n,i
    n=size(x)
    forall(i=1:n) x(i)=(i-1)/real(n-1,dp)
    if (s>0.0_dp) then
        x=x**s*(xmax-xmin)+xmin
        if (s==1.0_dp) then
!           PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
        else
!           PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
        end if
    else
        if (s==-1.0_dp) then
            c=xmax-xmin+1
!       ELSEIF (s==0.0_dp) THEN
!           IF (xmin>0.0_dp) THEN
!               c=xmax/xmin
!           ELSE
!               STOP 'grid: can not use logarithmic spacing for nonpositive values'
!           END IF
        else
            c=-s
        end if
!       PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
        x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
    end if
end subroutine grid

subroutine check(func,info)
    implicit none
    character(len=*), intent(in) :: func
    integer, intent(in) :: info
    if (info<0) then
        print '(a,a,a,i2,a)', 'Error in ',func,': the ',-info,'-th parameter had an illegal value'
    end if
    if (info>0) then
        select case(func)
        case('geev')
            print '(a,a,a)', 'Error in ',func,': the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed.'
        case('getrf')
            print '(a,a,a)', 'Error in ',func,': the factorization has been completed, but U is exactly singular.'
        case default
            print '(a,a,a,i2)', 'Error in ',func,': info=',info
        end select
    end if
    if (info/=0) stop 'program terminated by check'
end subroutine check

subroutine binom(f)
! Binomial probability mass function with p=1/2
    implicit none
    real(dp), dimension(:), intent(out) :: f
    integer :: n,k
    n=size(f)
    f(1)=2.0_dp**(1-n)
    do k=1,n-1
        f(k+1)=f(k)*(n-k)/k
    end do
end subroutine binom

real(dp) function brent(ax,bx,cx,func,xmin)
    real(dp), intent(in) :: ax,bx,cx
    real(dp), intent(out) :: xmin
    interface
        real(dp) function func(x)
        use kinds
        implicit none
        real(dp), intent(in) :: x
        end function func
    end interface
    integer, parameter :: ITMAX=100
    real(dp), parameter :: TOL=sqrt(epsilon(ax)),CGOLD=0.381966011250105_dp,ZEPS=1.0e-3_dp*epsilon(ax)
    integer :: iter
    real(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX
        xm=0.5_dp*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0_dp*tol1
        if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
            xmin=x
            brent=fx
            return
        end if
        if (abs(e) > tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0_dp*(q-r)
            if (q > 0.0) p=-p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                p <= q*(a-x) .or. p >= q*(b-x)) then
                e=merge(a-x,b-x, x >= xm )
                d=CGOLD*e
            else
                d=p/q
                u=x+d
                if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
            end if
        else
            e=merge(a-x,b-x, x >= xm )
            d=CGOLD*e
        end if
        u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
        fu=func(u)
        if (fu <= fx) then
            if (u >= x) then
                a=x
            else
                b=x
            end if
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
        else
            if (u < x) then
                a=u
            else
                b=u
            end if
            if (fu <= fw .or. w == x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
                v=u
                fv=fu
            end if
        end if
    end do
    stop 'brent: exceed maximum iterations'
end function brent

end module markov_chain_approx
