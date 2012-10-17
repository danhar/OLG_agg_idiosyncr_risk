module fun_locate
! overloaded f_locate: one for vectors only, one for scalars only.
    use kinds
    implicit none

    interface f_locate
        module procedure f_locate_vec, f_locate_scalar
    end interface f_locate

contains

pure function f_locate_vec(xgrid,x)
! Search ordered table, return index of lower bound. Returns ju-1 if x>xgrid(ju-1)!
! Works only for vectors!

	real(dp) ,dimension(:)       ,intent(in) :: xgrid		! grid
	real(dp) ,dimension(:)       ,intent(in) :: x			! values to be searched
	integer  ,dimension(size(x)) 			 :: f_locate_vec! indeces
	integer  ,dimension(size(x))			 :: jl,ju,jm

    jl=1
    ju=size(xgrid)
    do
        if (all(ju-jl<=1)) exit
        where (ju-jl>1)
	        jm=(ju+jl)/2
	        where (x>=xgrid(jm))
	            jl=jm
	        elsewhere
	            ju=jm
	        end where
	    end where
    end do
	f_locate_vec=jl

end function f_locate_vec


pure function f_locate_scalar(xgrid,x,method_o)
! Search ordered table, return index of lower bound. In 'default', returns ju-1 if x>xgrid(ju-1)!
! This version works only for scalar x!
! Just to make sure also included original NR version, method='numrec'

	real(dp), dimension(:), intent(in) 	 :: xgrid			! grid
	real(dp), intent(in) 				 :: x				! scalar(!) value to be searched
	character(len=*),intent(in),optional :: method_o
	integer						 		 :: f_locate_scalar	! index
	character(len=7)                     :: method
	integer								 :: jl, ju, jm, n
	logical :: ascnd

    if (.not.present(method_o)) then
        method = 'default'
    else
        method = method_o
    endif

	if (method=='default')	then 		! Returns ju-1 if x>xgrid(ju)!
	    jl=1
	    ju=size(xgrid)
	    do
	        if (ju-jl<=1) exit
	        jm=(ju+jl)/2
	        if (x>=xgrid(jm)) then
	            jl=jm
	        else
	            ju=jm
	        endif
	    end do
	f_locate_scalar=jl

	else ! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1045. Returns ju iff x>xgrid(ju)!
		if (method /= 'numrec') return
		n=size(xgrid)
		ascnd = (xgrid(n) >= xgrid(1))
		jl=0
		ju=n+1
		do
			if (ju-jl <= 1) exit
			jm=(ju+jl)/2
			if (ascnd .eqv. (x >= xgrid(jm))) then
				jl=jm
			else
				ju=jm
			end if
		end do
		if (x == xgrid(1)) then
			f_locate_scalar=1
		else if (x == xgrid(n)) then
			f_locate_scalar=n-1
		else
			f_locate_scalar=jl
		end if
	endif

end function f_locate_scalar

end module fun_locate
