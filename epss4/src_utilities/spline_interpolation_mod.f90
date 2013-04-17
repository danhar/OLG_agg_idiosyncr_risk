module spline_interpolation_mod
    use kinds
    implicit none

    interface tridag
        module procedure tridag_par
    end interface

contains

pure subroutine spline(x,y,yp1,ypn,y2,err)
    implicit none
    real(dp), dimension(:), intent(in) :: x,y
    real(dp), intent(in) :: yp1,ypn
    real(dp), dimension(:), intent(out) :: y2
    logical ,intent(out) :: err
    integer :: n
    real(dp), dimension(size(x)) :: a,b,c,r

    err = .false.
    if (size(x)==size(y) .and. size(y) == size(y2)) then
        n = size(x)
    else
        err = .true.
        return
    endif
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30_dp) then
        r(1)=0.0
        c(1)=0.0
    else
        r(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        c(1)=0.5
    end if
    if (ypn > 0.99e30_dp) then
        r(n)=0.0
        a(n)=0.0
    else
        r(n)=(-3.0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
        a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
end subroutine spline

pure function splint(xa,ya,y2a,x)
    use fun_locate
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), intent(in) :: x
    real(dp) :: splint
    integer :: khi,klo,n
    real(dp) :: a,b,h
    if (size(xa)==size(ya) .and. size(ya)==size(y2a)) then
        n=size(xa)
    else
        return
    endif
    klo=max(min(f_locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) return !call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
end function splint

pure recursive subroutine tridag_par(a,b,c,r,u)
    real(dp), dimension(:), intent(in) :: a,b,c,r
    real(dp), dimension(:), intent(out) :: u
    integer, parameter :: NPAR_TRIDAG=4
    integer :: n,n2,nm,nx
    real(dp), dimension(size(b)/2) :: y,q,piva
    real(dp), dimension(size(b)/2-1) :: x,z
    real(dp), dimension(size(a)/2) :: pivc
    n=size(a)+1
    if (n < NPAR_TRIDAG) then
        call tridag_ser(a,b,c,r,u)
    else
!        if (maxval(abs(b(1:n))) == 0.0) call nrerror('tridag_par: possible singular matrix')
        n2=size(y)
        nm=size(pivc)
        nx=size(x)
        piva = a(1:n-1:2)/b(1:n-1:2)
        pivc = c(2:n-1:2)/b(3:n:2)
        y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
        q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
        if (nm < n2) then
            y(n2) = b(n)-piva(n2)*c(n-1)
            q(n2) = r(n)-piva(n2)*r(n-1)
        end if
        x = -piva(2:n2)*a(2:n-2:2)
        z = -pivc(1:nx)*c(3:n-1:2)
        call tridag_par(x,y,z,q,u(2:n:2))
        u(1) = (r(1)-c(1)*u(2))/b(1)
        u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
            -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
        if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    end if
end subroutine tridag_par

pure subroutine tridag_ser(a,b,c,r,u)
    implicit none
    real(dp), dimension(:), intent(in) :: a,b,c,r
    real(dp), dimension(:), intent(out) :: u
    real(dp), dimension(size(b)) :: gam
    integer :: n,j
    real(dp) :: bet
    n=size(a)+1
    bet=b(1)
    !if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
     !   if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 2')
        u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    end do
end subroutine tridag_ser

end module spline_interpolation_mod
