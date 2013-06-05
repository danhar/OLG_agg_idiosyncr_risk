module alg_gs_mod
    implicit none
contains

subroutine alg_gs(fname,fvec,x,n,check,gswght,maxits,tolf)
! fixed point iteration
use kinds, only: dp

implicit none

integer :: i,its,n,maxits
real(dp) :: tolf,f,fold,gswght
real(dp), dimension(n) :: x,fvec,xold,p,xn
real(dp), parameter :: eps=epsilon(x),tolx=1.0e-8
logical :: check
integer,parameter::maxlns=5
real(dp),parameter::stpfac=10.0_dp

interface
    function fname(x)
        use kinds, only: dp
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x))       :: fname
    end function fname
end interface

its=0
fvec=fname(x)

! Check convergence
f=maxval(abs(fvec(:)))
print*,'Error in iteration # 1: ', f


! Test for initial guess being a root.
if (f < tolf) then			
	check=.false.
	return
endif

do its=1,maxits										! Start of iteration loop
	xold(:)=x(:)													! Store x, F, and f.
	fold=f
    
    p=-gswght*fvec
    call lnsrch_simple(xold,fold,p,x,n,f,fvec,stpfac,check,fname,maxlns)
    
    print*,'Error in iteration # ',its+1, ':', f

	if (f < tolf) then							! Test for convergence on function values.
		check=.false.
		return
	elseif ( maxval(abs(x(:)-xold(:))) < tolx) then					! Test for convergence on dx: Here: absolute deviation 
		check=.false.
		return				
	endif
	
end do
check=.true.
print*,'MAXITS exceeded in ALG_GS'


contains

! ---------------------------------------------
subroutine lnsrch_simple(xold,fold,p,x,n,f,fvec,stpfac,check,fname,maxlns)

use kinds, only: dp
implicit none
integer,intent(in):: n,maxlns
real(dp),intent(in)::stpfac
real(dp),intent(inout)::xold(n),fold,p(n),x(n),f,fvec(n)
logical,intent(inout)::check
real(dp)::alam
real(dp),parameter::alamin=0.01_dp
integer::its
interface
    function fname(x)
        use kinds, only: dp
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x))       :: fname
    end function fname
end interface


check=.false.

alam=1.0_dp
do its=1,maxlns
	if (its>1) print*,'line search iteration # ',its-1

    x(:)=xold(:)+alam*p(:)
    fvec=fname(x)
    f=maxval(abs(fvec(:)))
    ! print*, 'current alam,f,fold,maxlns:', alam,f,fold,maxlns
    
	if (alam < alamin) then										! Convergence on x
		check=.true.
		return
	elseif (f < fold) then						! Sufficient function decrease.
		return
    else														! Backtrack.
	    alam=alam/stpfac
    endif
end do 

end subroutine lnsrch_simple
! ---------------------------------------------

end subroutine alg_gs

end module alg_gs_mod
