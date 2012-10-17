module numrec_utils
! Adapted from Press et al. (2002), Numerical Recipes in F90 (NR), p. 1364ff

	use kinds
	use fun_kronprod
	implicit none

    interface put_diag
        module procedure put_diag_vec, put_diag_scalar
    end interface

contains
	! ---------------------------------------------
	pure function get_diag(mat)
	! input: matrix of arbitrary size
	! output: diagonal elements of matrix

		real(dp),dimension(:,:),intent(in) :: mat
		real(dp), dimension(min(size(mat,1),size(mat,2))) :: get_diag
		integer j

		do j=1,min(size(mat,1),size(mat,2))
			get_diag(j)=mat(j,j)
		end do

	end function

	! ---------------------------------------------
	pure function lower_triangle(j,k ,extra)
	! (returns a lower triangular mask)

		integer, intent(in) :: j,k
		integer, optional, intent(in) :: extra
		logical, dimension(j,k) :: lower_triangle

		integer :: n,jj,kk
		n=0
		if (present(extra)) n=extra
		do jj=1,j
			do kk=1,k
				lower_triangle(jj,kk)= (kk-jj < n)
			end do
		end do

	end function lower_triangle

	! ---------------------------------------------
	pure function outerprod(a,b)
	! takes outer product of two vectors a and b
		real(dp), dimension(:), intent(in) :: a,b
		real(dp), dimension(size(a),size(b)) :: outerprod

		outerprod = spread(a,dim=2,ncopies=size(b)) * &
			spread(b,dim=1,ncopies=size(a))

	end function outerprod

	! ---------------------------------------------
	pure function ifirstloc(mask)
	! Returns the index (subscript value) of the first location, in a one-dimensional
	! logical mask, that has the value .TRUE., or returns size(mask)+1 if all
	! components of mask are .FALSE.
	! Note that while the reference implementation uses a do-loop, the function is
	! parallelized in nrutil by instead using the merge and maxloc intrinsics.
	! Reference implementation:

		logical, dimension(:), intent(in) :: mask
		integer:: ifirstloc
		integer:: i
		do i=1,size(mask)
			if (mask(i)) then
				ifirstloc=i
				return
			end if
		end do
		ifirstloc=i
	end function ifirstloc

	! ---------------------------------------------
	subroutine put_diag_vec(diag,mat)
	! Sets the diagonal of matrix mat equal to the argument diag,
	! a vector whose size must be the smaller of the two dimensions of matrix mat.
		real(dp), dimension(:), intent(in) :: diag
		real(dp), dimension(:,:), intent(inout) :: mat
		integer:: j,n,m
		n=size(diag)
		m=min(size(mat,1),size(mat,2))
		if (n==m) then
			do j=1,n
				mat(j,j)=diag(j)
			end do
		else
			print*, 'error in put_diag: n not equal to m'
			return
		endif
	end subroutine put_diag_vec

	! ---------------------------------------------
    pure subroutine put_diag_scalar(diag,mat)
    ! Sets the diagonal of matrix mat equal to the argument diag, a scalar
        real(dp), intent(in) :: diag
        real(dp), dimension(:,:), intent(inout) :: mat
        integer:: j,m
        m=min(size(mat,1),size(mat,2))
        forall (j=1:m) mat(j,j)=diag
    end subroutine put_diag_scalar

	! ---------------------------------------------
	pure subroutine unit_matrix(mat)
    ! Sets the diagonal components of mat to unity, all other components to zero.
    ! When mat is square, this will be the unit matrix; otherwise, a unit matrix
    ! with appended rows or columns of zeros.
		real(dp), dimension(:,:), intent(out) :: mat
		integer:: i,n
		n=min(size(mat,1),size(mat,2))
		mat(:,:)=0.0_dp
		do i=1,n
			mat(i,i)=1.0_dp
		end do
	end subroutine unit_matrix

	! ---------------------------------------------
	pure subroutine qrdcmp(a,c,d,sing)
	! Constructs the QR decomposition of the n-times-n matrix a. The upper triangular matrix R is
	! returned in the upper triangle of a, except for the diagonal elements of R, which are returned
	! in the n-dimensional vector d. The orthogonal matrix Q is represented as a product of n-1
	! Householder matrices Q1 . . .Qn-1, where Qj = 1 - uj . uj/cj. The ith component of uj
	! is zero for i = 1, . . . , j - 1 while the nonzero components are returned in a(i,j) for
	! i = j, . . . ,n. sing returns as true if singularity is encountered during the decomposition,
	! but the decomposition is still completed in this case.

		real(dp), dimension(:,:), intent(inout) :: a
		real(dp), dimension(:), intent(out) :: c,d
		real(dp)::vabs
		logical, intent(out) :: sing
		integer:: k,n
		real(dp) :: scal,sigma

		sing=.false.
		n=size(a,1)
		do k=1,n-1
			scal=maxval(abs(a(k:n,k)))
			if (scal == 0.0) then				! Singular case.
				sing=.true.
				c(k)=0.0
				d(k)=0.0
			else								! Form Qk and Qk times A.
				a(k:n,k)=a(k:n,k)/scal
				vabs=sqrt(dot_product(a(k:n,k),a(k:n,k)))
				sigma=sign(vabs,a(k,k))
				a(k,k)=a(k,k)+sigma
				c(k)=sigma*a(k,k)
				d(k)=-scal*sigma
				a(k:n,k+1:n)=a(k:n,k+1:n)-outerprod(a(k:n,k),matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
			end if
		end do
		d(n)=a(n,n)
		if (d(n) == 0.0) sing=.true.

	end subroutine qrdcmp

	! ---------------------------------------------
	pure subroutine qrupdt(r,qt,u,v)
	! Given the QR decomposition of some nxn matrix, calculates the QR decomposition of
	! the matrix Q (R + kron(u,v)). Here r and qt are nxn matrices, u and v are n-dimensional
	! vectors. Note that QT is input and returned in qt.
		real(dp), dimension(:,:), intent(inout) :: r,qt
		real(dp), dimension(:), intent(inout) :: u
		real(dp), dimension(:), intent(in) :: v
		integer:: i,k,n

		n=size(r,1)
		k=n+1-ifirstloc(u(n:1:-1) /= 0.0)		! Find largest k such that u(k) .= 0.
		if (k < 1) k=1

		do i=k-1,1,-1						! Transform R + kron(u,v) to upper Hessenberg.
			call rotate(r,qt,i,u(i),-u(i+1))
			u(i)=pythag(u(i),u(i+1))
		end do
		r(1,:)=r(1,:)+u(1)*v
		do i=1,k-1							! Transform upper Hessenberg matrix to upper triangular.
			call rotate(r,qt,i,r(i,i),-r(i+1,i))
		end do

	end subroutine qrupdt

	! ---------------------------------------------
	pure subroutine rotate(r,qt,i,a,b)
	! Given nxn matrices r and qt, carry out a Jacobi rotation on rows i and i+1 of each matrix.
	! a and b are the parameters of the rotation
		real(dp), dimension(:,:), target, intent(inout) :: r,qt
		integer, intent(in) :: i
		real(dp), intent(in) :: a,b
		real(dp), dimension(size(r,1)) :: temp
		integer:: n
		real(dp) :: c,fact,s

		n=size(r,1)
		if (a == 0.0) then							! Avoid unnecessary overflow or underflow.
			c=0.0
			s=sign(1.0,b)
		elseif (abs(a) > abs(b)) then
			fact=b/a
			c=sign(1.0/sqrt(1.0+fact**2),a)
			s=fact*c
		else
			fact=a/b
			s=sign(1.0/sqrt(1.0+fact**2),b)
			c=fact*s
		endif
		temp(i:n)=r(i,i:n)							! Premultiply r by Jacobi rotation.
		r(i,i:n)=c*temp(i:n)-s*r(i+1,i:n)
		r(i+1,i:n)=s*temp(i:n)+c*r(i+1,i:n)
		temp=qt(i,:)								! Premultiply qt by Jacobi rotation.
		qt(i,:)=c*temp-s*qt(i+1,:)
		qt(i+1,:)=s*temp+c*qt(i+1,:)

	end subroutine rotate

	! ---------------------------------------------
	pure function pythag(a,b)
		real(dp), intent(in) :: a,b
		real(dp) :: pythag
		real(dp) :: absa,absb
		absa=abs(a)
		absb=abs(b)
		if (absa > absb) then
			pythag=absa*sqrt(1.0+(absb/absa)**2)
		else
			if (absb == 0.0) then
				pythag=0.0
			else
				pythag=absb*sqrt(1.0+(absa/absb)**2)
			endif
		endif
	end function pythag

	! ---------------------------------------------
	pure subroutine rsolv(a,d,b)
	! Solves the set of n linear equations R � x = b, where R is an upper triangular matrix stored
	! in a and d. The n�n matrix a and the vector d of length n are input as the output of the
	! routine qrdcmp and are not modified. b is input as the right-hand-side vector of length n,
	! and is overwritten with the solution vector on output.
		real(dp), dimension(:,:), intent(in) :: a
		real(dp), dimension(:), intent(in) :: d
		real(dp), dimension(:), intent(inout) :: b
		integer:: i,n
		n=size(a,1)
		b(n)=b(n)/d(n)
		do i=n-1,1,-1
			b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i)
		end do
	end subroutine rsolv

	! ---------------------------------------------
    pure function arth(first,increment,n)
        real(dp), intent(in) :: first,increment
        integer, intent(in) :: n
        real(dp), dimension(n) :: arth
        integer :: k,k2
        real(dp) :: temp, NPAR_ARTH, NPAR2_ARTH
        NPAR_ARTH=16
        NPAR2_ARTH=8
        if (n > 0) arth(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
                arth(k)=arth(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
                arth(k)=arth(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                if (k >= n) exit
                k2=k+k
                arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                temp=temp+temp
                k=k2
            end do
        end if
    end function arth

end module numrec_utils
