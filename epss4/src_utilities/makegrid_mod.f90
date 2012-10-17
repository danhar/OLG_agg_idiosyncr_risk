module makegrid_mod
    use kinds

    implicit none
    private
    public MakeGrid

contains

    pure function MakeGrid(lb,ub,n,curv,method)
        real(dp)     ,intent(in)           :: lb,ub       ! lower and upper bound
        integer      ,intent(in)           :: n
        real(dp)     ,intent(in) ,optional :: curv        !'curvature'
        character(len=*) ,intent(in) ,optional :: method
        real(dp) ,dimension(n)             :: MakeGrid

    ! The following allows to supply curv and chebyshev, curv is then not used.
    if (present(method)) then
        if(method=='chebyshev') then
            MakeGrid = makegrid_chebyshev(n,lb,ub)
	    else
            MakeGrid = makegrid_exponential(lb,ub,n,curv) ! optional curv cascades
	    endif
    else
         MakeGrid = makegrid_exponential(lb,ub,n,curv)
    endif

    end function MakeGrid
!-------------------------------------------------------------------------------

	pure function makegrid_exponential(x1,x2,n,curv) result(grid)
	! function for grid generation
	! (optional) curv>1 -> more points close to x1; 0<curv<1 more points close to x2; default curv=1.0
		real(dp) ,intent(in)			 :: x1,x2		! lower and upper bound
		integer  ,intent(in)			 :: n
		real(dp) ,intent(in)   ,optional :: curv		!'curvature'
		real(dp) ,dimension(n)           :: grid
		integer							 :: i
		real(dp)						 :: scalefact, c

	    if (.not.present(curv)) then
	        c = 1.0
	    else
	        c = curv
	    endif

		scalefact = x2-x1
		grid(n)   = x2
		grid(1)   = x1		! this way, if n=1, we return only lower bound
		do i=2,n-1
		    grid(i) = x1 + scalefact*((i-1.0)/(n-1.0))**c
		end do

	end function makegrid_exponential
!-------------------------------------------------------------------------------

	pure function makegrid_chebyshev(m,lb,ub) result(cheb_extrema)
		! Computes the chebyshev nodes  (extrema) on an interval [lb,ub]. The input variables
		! lb and ub are optional. If they are not present, then the standard nodes on [-1,1]
		! are computed.
		integer, intent(in) :: m
		real(dp), intent(in), optional :: lb, ub
		real(dp) :: cheb_extrema(m)
		integer ::i             ! loop index
		real(dp), parameter     :: pi = 3.141592653589793238462643383279502884197_dp    ! number pi

		! Compute the 'm' Chebyshev extrema in [-1,1]
		do i = 1, m
		    cheb_extrema(i) = - COS( pi*(i-1.0) /real(m-1,dp) )
		end do

		! Check if need to move nodes to [lb,ub]
		if ( PRESENT( lb ) .and. PRESENT( ub ) ) then
		    cheb_extrema = ( cheb_extrema + 1.0 ) * 0.5_dp * ( ub - lb ) + lb
		end if

	end function makegrid_chebyshev

end module makegrid_mod
