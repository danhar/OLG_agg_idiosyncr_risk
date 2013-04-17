module spline_interpolation_mod
use params_mod, only: dp
contains

pure subroutine spline(x,y,yp1,ypn,y2,err)

!    use nr, only : tridag
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
!    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
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

end module spline_interpolation_mod
