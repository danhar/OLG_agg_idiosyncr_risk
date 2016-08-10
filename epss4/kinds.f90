! Copyright (C) 2016 Daniel Harenberg - All Rights Reserved
module kinds
! Type kind parameters
    implicit none
	integer, parameter ::								&
		! for integers, the argument n specifies the range +/-10**n
		! the naming i? corresponds to the Intel Fortran Compiler integer kinds: ? is bytes used
		! The comment gives the precise range resulting from the byte size
    	i1kind = selected_int_kind(2),                  &  ! 2**(8*1-1)=+/-128
    	i2kind = selected_int_kind(4),					&  ! 2**(8*2-1)=+/-32768
    	i4kind = selected_int_kind(8),					&  ! 2**(8*4-1)=+/-2.1e9
    	i8kind = selected_int_kind(16),					&  ! 2**(8*8-1)=+/-9.2e18


		! comments are approx machine precision and exponent range for the Intel Fortran compiler 11.1
    	sp = kind(1.0),									&	! 1e-8,  +/-38
    	dp = selected_real_kind(2*precision(1.0_sp)), 	&	! 1e-16, +/-308
    	qp = selected_real_kind(2*precision(1.0_dp)),	&	! 1e-40, +/-4932
        pr = dp,                                        &   ! preferred real precision

    	! minimum logical size. Intel Fortran Compiler has kind=4 as standard, using 4 bytes
    	sl = kind(.true.)	! 'single prec logical', equiv to logical(1). note that default logical is recommended in the manual

    ! literal constants (pi, etc)

end module kinds
