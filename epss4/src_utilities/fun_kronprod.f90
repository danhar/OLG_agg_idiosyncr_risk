module fun_kronprod
    implicit none
contains
! kronecker product
pure function f_kronprod(A,B)
    use kinds
    real(dp), dimension(:,:), intent(in) :: A,B
	real(dp), dimension(size(A,1)*size(B,1),size(A,2)*size(B,2)) :: f_kronprod
	integer:: n,m,p,q,i,j,k,l

	m=size(A,1)
	n=size(A,2)
	p=size(B,1)
	q=size(B,2)
	f_kronprod(:,:)=0.0_dp
	do i=1,m						! go down rows of A
		k=(i-1)*p+1					! fill in rows of B
		do j=1,n					! go down columns of A
			l=(j-1)*q+1				! fill in columns of B
			f_kronprod(k:k+p-1,l:l+q-1)=A(i,j)*B
		end do
	end do

end function f_kronprod
end module fun_kronprod
