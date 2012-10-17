module markovchain_mod
    implicit none
contains

function MarkovChain(T,nt,start_value, seed) result(chain)
    use kinds
    use fun_locate

    real(dp) ,dimension(:,:) ,intent(in) :: T     ! Transition matrix
    integer                  ,intent(in) :: nt	  ! number of simulation periods
    integer  ,optional       ,intent(in) :: start_value  ! initial state
    integer  ,dimension(:) ,optional ,intent(in) :: seed
    integer  ,dimension(nt)              :: chain ! markov chain containing indeces of states
    real(dp) ,dimension(nt)		         :: randn ! holds real random numbers
    real(dp) ,dimension(size(T,1), size(T,2)+1) :: cumul ! cumulative probabilities of T
    integer  ,dimension(size(T,1))       :: states
    integer                              :: i, nz
    integer(kind(nt))                    :: tc

!    if(present(reset_seed)) then
!        if (reset_seed) call random_seed
!    endif

    if (present(seed)) call random_seed(put=seed)

    if (present(start_value)) then
        chain(1) = start_value
    else
        chain(1) = 1
    endif

    nz=size(T,1)
    states=(/(i,i=1,nz)/)
    cumul(:,1)=0.0
    do i =1, nz
        cumul(:,i+1)= cumul(:,i)+T(:,i)
    enddo

	call random_number(randn)

    do tc=2,nt
        chain(tc) = f_locate(cumul(chain(tc-1),:),randn(tc))
    enddo

end function MarkovChain
end module markovchain_mod
