module quadrature_mod

use kinds

implicit none
private
public quadrature, gauss_hermite_quad, gauss_legendre_quad, romberg_quad, &
       romberg_quad_open, trapezoid_quad, change_of_standard_normal_rv

interface arth
    module procedure arth_d, arth_i
end interface

contains
!-------------------------------------------------------------------------------

pure subroutine quadrature(nodes, weights, n_o, quadrature_method_o, factor_o)
! Alternatively, nodes and weights could be type components.
use global_constants   ,only: pi
! romberg_quad, romberg_quad_open take function as input.
real(dp), dimension(:) ,allocatable ,intent(out) :: nodes, weights
integer, intent(in) ,optional       :: n_o, factor_o
character(*), intent(in), optional  :: quadrature_method_o
integer                             :: flag, factor, n
character(:), allocatable           :: quadrature_method
character(*)           ,parameter   :: procedure_name    = 'quadrature'

if (present(n_o)) then
    n = n_o
else
    n = 3
endif
if (present(quadrature_method_o)) then
    quadrature_method = quadrature_method_o
else
    quadrature_method = 'gauss_hermite'
endif
if (present(factor_o)) then
    factor = factor_o
else
    factor = 1
endif

allocate(nodes(n*factor), weights(n*factor))

select case (quadrature_method)
case('gauss_hermite','gh')
    ! For standard-normal random variable:
    call gauss_hermite_quad(nodes, weights, flag)
    ! For expectations of standard-normal:
    weights = 1.0_dp/sqrt(pi) * weights
case('gauss_legendre','gl')
    call gauss_legendre_quad(-2._dp,2._dp,nodes,weights, flag)
    ! Assume uniform distribution and puts directly into weights:
    weights= weights/real(size(weights),dp)
case('uniform')
    ! These have to sum to one, else problembs with RS or EZ prefs:
    weights= 1.0_dp/real(size(weights),dp)
end select

end subroutine quadrature
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine change_of_standard_normal_rv(mean_o, standard_deviation, lognormal_o, nodes)
real(dp), intent(in), optional :: mean_o
real(dp), intent(in) :: standard_deviation
logical, intent(in), optional :: lognormal_o
real(dp), dimension(:) ,intent(inout) :: nodes
real(dp) :: mean
logical :: lognormal
character(*),parameter::procedure_name='change_of_standard_normal_variable'

if (present(lognormal_o)) then ! This check needs to come first.
    lognormal = lognormal_o
else
    lognormal = .false.
endif
if (present(mean_o)) then
    mean = mean_o
else
    if (lognormal) then
        mean = -(standard_deviation**2.0_dp/2.0_dp)
    else
        mean = 0.0
    endif
endif

nodes = mean + sqrt(2.0_dp) * standard_deviation *nodes
if (lognormal) nodes = exp(nodes)

end subroutine change_of_standard_normal_rv
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure subroutine gauss_hermite_quad(x,w, exitflag)
! Gauss Hermite-Quadrature, adapted from Press, et al. (1992), p. 1062.
! Returns abscissas and weights of N-point Gauss-Hermite quadrature formula.
! In contrast to Press et al., abscissas x returned in ascending order.
! exitflag: 1 if correct, -1 if too many iterations (its+1==maxits),
! -2 if  (size(x) .ne. size(w)).
    use global_constants, only: pi
    real(dp), dimension(:), intent(out) :: x,w
    integer               , intent(out) :: exitflag
    ! Relative precision in dp, see adjustment in Press, et al. (1992), p.
    real(dp), parameter :: eps=1.0e-14_dp, &
                           pim4=0.7511255444649425_dp ! 1/pi^{1/4}
    integer :: its,j,m,n
    integer, parameter :: maxit=10
    real(dp) :: anu
    real(dp), parameter :: c1=9.084064e-01_dp,c2=5.214976e-02_dp,&
        c3=2.579930e-03_dp,c4=3.986126e-03_dp
    real(dp), dimension((size(x)+1)/2) :: rhs,r2,r3,theta
    real(dp), dimension((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    logical, dimension((size(x)+1)/2) :: unfinished

    exitflag = 1
    if (size(x) .ne. size(w)) exitflag = -2

    n=size(x)

    m=(n+1)/2
    anu=2.0_dp*n+1.0_dp
    rhs=arth(3,4,m)*pi/anu
    r3=rhs**(1.0_dp/3.0_dp)
    r2=r3**2
    theta=r3*(c1+r2*(c2+r2*(c3+r2*c4)))
    z=sqrt(anu)*cos(theta)
    unfinished=.true.
    do its=1,maxit
        where (unfinished)
            p1=pim4
            p2=0.0
        end where
        do j=1,n
            where (unfinished)
                p3=p2
                p2=p1
                p1=z*sqrt(2.0_dp/j)*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
            end where
        end do
        where (unfinished)
            pp=sqrt(2.0_dp*n)*p2
            z1=z
            z=z1-p1/pp
            unfinished=(abs(z-z1) > eps)
        end where
        if (.not. any(unfinished)) exit
    end do
    if (its == maxit+1) exitflag = -1
    x(1:m)=-z
    x(n:n-m+1:-1)=z
    w(1:m)=2.0_dp/pp**2
    w(n:n-m+1:-1)=w(1:m)
end subroutine gauss_hermite_quad
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure subroutine gauss_legendre_quad(x1,x2,x,w,exitflag)
! Gauss-Legendre Quadrature, adapted from Press, et al. (1998), p.145 & p.1059.
! Returns abscissas and weights of N-point Gauss-Legendre quadrature formula.
! exitflag: 1 if correct, -1 if too many iterations (its+1==maxits),
! -2 if  (size(x) .ne. size(w)).
    use global_constants, only: pi
    real(dp), intent(in) :: x1,x2 ! Lower and upper limits of integration
    real(dp), dimension(:), intent(out) :: x,w ! abscissas and weights
    integer,                intent(out) :: exitflag
    real(dp), parameter :: eps=1.0e-14_dp ! relative precision
    integer :: its,j,m,n
    integer, parameter :: maxit=10
    real(dp) :: xl,xm
    real(dp), dimension((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    logical, dimension((size(x)+1)/2) :: unfinished

    exitflag = 1
    if (size(x) .ne. size(w)) exitflag = -2

    n=size(x)
    m=(n+1)/2
    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)
    z=cos(pi*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
    unfinished=.true.
    do its=1,maxit
        where (unfinished)
            p1=1.0
            p2=0.0
        end where
        do j=1,n
            where (unfinished)
                p3=p2
                p2=p1
                p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
            end where
        end do
        where (unfinished)
            pp=n*(z*p1-p2)/(z*z-1.0_dp)
            z1=z
            z=z1-p1/pp
            unfinished=(abs(z-z1) > eps)
        end where
        if (.not. any(unfinished)) exit
    end do
    if (its == maxit+1) exitflag = -1
    x(1:m)=xm-xl*z
    x(n:n-m+1:-1)=xm+xl*z
    w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
    w(n:n-m+1:-1)=w(1:m)
end subroutine gauss_legendre_quad
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure function arth_d(first,increment,n)
! From Numerical Recipes nrutils.f90, see Press et al. 1992, p. 1371.
! Routine relating to polynomials and recurrences.
! Array function returning arithmetic progression.
! Used, e.g., by gauss_hermite_quad.
    real(dp), intent(in) :: first,increment
    integer, intent(in) :: n
    real(dp), dimension(n) :: arth_d
    integer :: k,k2
    real(dp) :: temp
    integer ,parameter :: NPAR_ARTH=16,NPAR2_ARTH=8

    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth_d(k)=arth_d(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth_d(k)=arth_d(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
end function arth_d

pure function arth_i(first,increment,n)
    integer, intent(in) :: first,increment,n
    integer, dimension(n) :: arth_i
    integer :: k,k2,temp
    integer ,parameter :: NPAR_ARTH=16,NPAR2_ARTH=8

    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth_i(k)=arth_i(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
end function arth_i
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure real(dp) function romberg_quad(func,a,b)
! Romberg integration of order 2K, where K=2 is Simpson's rule.
! Adapted from Numerical Recipes, Vol. 2, p. 1054.
    real(dp), intent(in) :: a,b
    interface
        pure function func(x)
        use kinds ,only: dp
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x)) :: func
        end function func
    end interface
    ! jmax is max number of steps, k is number of points used in extrapolation.
    integer, parameter :: jmax=20,jmaxp=jmax+1,k=5,km=k-1
    ! eps is fractional accuracy, determined by extrapolation error estimate.
    real(dp), parameter :: eps=1.0e-6_dp
    ! h, s store successive trapezoidal approximations and relative step sizes.
    real(dp), dimension(jmaxp) :: h,s
    real(dp) :: dqromb
    integer :: j, exitflag
    h(1)=1.0
    do j=1,jmax
        call trapezoid_quad(func,a,b,s(j),j)
        if (j >= k) then
            call polint(h(j-km:j),s(j-km:j),0.0_dp,romberg_quad,dqromb,exitflag)
            if (abs(dqromb) <= eps*abs(romberg_quad)) return
        end if
        s(j+1)=s(j)
        h(j+1)=0.25_dp*h(j)
    end do
    ! call nrerror('qromb: too many steps')
end function romberg_quad
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure subroutine polint(xa,ya,x,y,dy,exitflag)
! From Numerical Recipes, Vol 2, see Press et al. 1992, p. 1043.
! Routine for polynomial interpolation or extrapolation.
! Used, e.g., by romberg_quad
! exitflag: 1 if correct, -1 calculation failure, -2 if (size(xa).ne.size(ya)).
    real(dp), dimension(:), intent(in) :: xa,ya
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y,dy
    integer , intent(out) :: exitflag
    integer :: m,n,ns, imin(1)
    real(dp), dimension(size(xa)) :: c,d,den,ho

    exitflag = 1 ! no error
    if (size(xa) .ne. size(ya)) exitflag = -2

    n=size(xa)
    c=ya
    d=ya
    ho=xa-x
    imin=minloc(abs(x-xa))
    ns=imin(1)
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
        den(1:n-m)=ho(1:n-m)-ho(1+m:n)
        if (any(den(1:n-m) == 0.0)) exitflag = -1 ! polint: calculation failure
        den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
        d(1:n-m)=ho(1+m:n)*den(1:n-m)
        c(1:n-m)=ho(1:n-m)*den(1:n-m)
        if (2*ns < n-m) then
            dy=c(ns+1)
        else
            dy=d(ns)
            ns=ns-1
        end if
        y=y+dy
    end do
end subroutine polint
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure subroutine trapezoid_quad(func,a,b,s,n)
! Quadrature using Trapezoidal rule.
! Adapted from Numerical Recipes, Vol. 2, p. 1052.
! Used, e.g., by romberg_quad
    real(dp), intent(in) :: a,b
    real(dp), intent(inout) :: s
    integer, intent(in) :: n
    interface
        pure function func(x)
        use kinds ,only: dp
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x)) :: func
        end function func
    end interface
    real(dp) :: del,fsum
    integer :: it
    if (n == 1) then
        s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
    else
        it=2**(n-2)
        del=(b-a)/it
        fsum=sum(func(arth(a+0.5_dp*del,del,it)))
        s=0.5_dp*(s+del*fsum)
    end if
end subroutine trapezoid_quad
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure real(dp) function romberg_quad_open(func,a,b,choose)
! Romberg quadrature for improper integrals, e.g., with open interval.
! choose is an integrating function, at the moment can take only midpnt
! Adapted from Numerical Recipes, Vol. 1+2, p. 137
real(dp), intent(in) :: a,b
interface
    pure function func(x)
    use kinds ,only: dp
    implicit none
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: func
    end function func

    pure subroutine choose(funk,aa,bb,s,n)
    use kinds ,only: dp
    implicit none
    real(dp), intent(in) :: aa,bb
    real(dp), intent(inout) :: s
    integer, intent(in) :: n
    interface
        pure function funk(x)
        use kinds ,only: dp
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x)) :: funk
        end function funk
    end interface
    end subroutine choose
end interface
integer, parameter :: jmax=14,jmaxp=jmax+1,k=5,km=k-1
real(dp), parameter :: eps=1.0e-6
real(dp), dimension(jmaxp) :: h,s
real(dp) :: dqromo
integer :: j, exitflag
h(1)=1.0
do j=1,jmax
    call choose(func,a,b,s(j),j)
    if (j >= k) then
        call polint(h(j-km:j),s(j-km:j),0.0_dp,romberg_quad_open,dqromo,exitflag)
        if (abs(dqromo) <= eps*abs(romberg_quad_open)) return
    end if
    s(j+1)=s(j)
    h(j+1)=h(j)/9.0_dp
end do
! call nrerror('qromo: too many steps')
end function romberg_quad_open
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

pure subroutine midpnt(func,a,b,s,n)
! Quadrature using midpoint rule
! can be passed as argument to romberg_quad_open
! Adapted from Numerical Recipes
    real(dp), intent(in) :: a,b
    real(dp), intent(inout) :: s
    integer, intent(in) :: n
    interface
        pure function func(x)
        use kinds ,only: dp
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x)) :: func
        end function func
    end interface
    real(dp) :: del
    integer :: it
    real(dp), dimension(2*3**(n-2)) :: x
    if (n == 1) then
        s=(b-a)*sum(func( (/0.5_dp*(a+b)/) ))
    else
        it=3**(n-2)
        del=(b-a)/(3.0_dp*it)
        x(1:2*it-1:2)=arth(a+0.5_dp*del,3.0_dp*del,it)
        x(2:2*it:2)=x(1:2*it-1:2)+2.0_dp*del
        s=s/3.0_dp+del*sum(func(x))
    end if
end subroutine midpnt
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

end module quadrature_mod
