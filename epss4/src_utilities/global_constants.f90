!*******************************************************************************
! Copyright (c) 2010-2017 Daniel Harenberg - All rights reserved.
!*******************************************************************************

module global_constants
    use kinds, only: dp
    implicit none

    ! parameters influencing code execution
    logical ,parameter :: debugging = .false.


    character(*) ,parameter :: &
    ! relative paths
                   path_to_previous_results = 'model_input/previous_results/' ,&

    ! scalar formats for writing to file
                            format_real     = '(es13.6,1x)'   , &
                            format_real_E3  = '(es13.6E3,1x)' , &
                            format_logi     = '(l1,1x)'       , &
                            fmt_char_real   = '(a,x,f7.4)'    , &
                            fmt_char_es     = '(a,x,es13.6)'  , &
                            fmt_char_esE3   = '(a,x,es13.5E3)', &
                            fmt_char_int    = '(a,x,i3)'      , &
                            fmt_char_logi   = '(a,2x,l1)'     , &
                            fmt_char_char   = '(a,2x,a)'      , &
                            fmt_tau_char    = '(f4.3)'

    ! literal constants (pi, etc)
    real(dp) ,parameter :: euler = 0.5772156649015328606065120900824024310422_dp
    real(dp) ,parameter :: pi    = 3.141592653589793238462643383279502884197_dp

end module global_constants
