!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for applying the Laplacian
!
!    \Delta = (d/dr)^2 + 1/r (d/dr) + 1/r^2 (d/dt)^2                                      (1)
!
!  to functions on an annulus.  Functions are represented via their values at a 
!  two-dimensional discretization grid which is the tensor product of two 
!  one-dimensional discretization grids. One, in the t variable, discretizes 
!  the exponential functions
!
!    { exp(ikt)  : -n <= k <= n }
!
!  and a second discretization grid --- this one in the r variable 0 --- discretizes
!  the polynomials
!
!    r^j                       for j= 0,1,...,m
!
!  The first is of length 2n+1 and the second of length m+1.  The fast Fourier 
!  transform, which is applied using Swarztrauber's fftpack library, is used to 
!  accelerate the operations.
!
!  This code is somewhat inefficient as the initialization data for the FFT
!  and other data needed to apply the operators are generated repeatedly rather 
!  than being generated once upon initialization and stored.  This is done to
!  make parallelization easier.
!
!  The following subroutines should be regarded as publicly callable:
!
!    helmop_init - return the tensor product discretization grid given n and m
!
!    helmop_delta - apply the operator (1) a user-specified function supplied
!      by its values at the two-dimensional discretization grid
!
!    helmop_helmholtz - apply the operator 
!
!           \Delta + lambda^2 ( 1 + q(r) ) 
!
!      with q(r) a radially symmetric function supplied by the user 
!      
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module helmop

use utils

contains


subroutine helmop_init(R1,R2,n,m,ts,rs)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)    :: ts(:)
double precision, allocatable, intent(out)    :: rs(:)
!
!  Construct the tensor product discretization grid used to represent
!  functions and perform any necessary initializations for the FFT.
!
!  Input parameters:
!    R1 - the inner radius of the annulus
!    R2 - the outer radius of the annulus
!    n - a parameter which determines the size of the discretization grid in the
!      t variable
!    m - a parameter which determines the size of the discretization grid in the
!      r variable
!
!  Output parameters:
!    ts - an array of length (2n+1) containing the discretization nodes in the
!      t variable
!    rs - an array of length (m+1) containing the discretization nodes in the
!      r variable
!

double precision, allocatable :: ts0(:)

pi = acos(-1.0d0)

allocate(ts(2*n+1),ts0(m+1),rs(m+1))

do i=1,2*n+1
ts(i) = (i-1)*2*pi/(2*n+1.0d0)
end do


a   = R1
b   = R2

do i=1,m+1
ts0(i) = (m+1-i)*pi/(m+0.0d0)
end do

! Chebyshev nodes
rs = (b-a)/2 * cos(ts0) +  (b+a)/2


end subroutine


subroutine helmop_helmholtz(dlambda,R1,R2,n,m,ts,rs,qs,x,y)
implicit double precision (a-h,o-z)
double precision :: ts(n),rs(m),qs(m)
double complex   :: x(2*n+1,m+1),y(2*n+1,m+1)
!
!  Apply the operator \Delta + q(r) I to a function represented by its values
!  at the discretization nodes.
!
!  Input parameters:
!    dlambda - the wavenumber of the operator
!    R1 - the inner radius of the annulus
!    R2 - the outer radius of the annulus
!    n - the parameter controlling the number of discretization nodes in the t
!      variable
!    m - the parameter controlling the number of discretizatoon nodes in the r
!      variables
!    ts - the discretization nodes in the t variable as returned by helmop_init
!    rs - the discretization nodes in the r variable as returned by helmop_init
!    qs - an array giving the value of q at the discretization nodes in the r
!       variable
!    x - an array dimensions (2n+1,m+1) whose (i,j) entry gives the value
!      of the input function at the point ( ts(i), rs(j) )
!
!  Ouput parameters:
!
!    y - an array dimensions (2n+1,m+1) whose (i,j) entry gives the value
!      of the output function at the point ( ts(i), rs(j) )
!


call helmop_delta(R1,R2,n,m,ts,rs,x,y)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,qval)
!$OMP DO
do j=1,m+1
qval    = qs(j)
!qval    = 0         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y(:,j)  = y(:,j) + x(:,j)*dlambda**2*(1+qval)
end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine


subroutine helmop_delta(R1,R2,n,m,ts,rs,x,y)
implicit double precision (a-h,o-z)
double precision :: ts(n),rs(m)
double complex   :: x(2*n+1,m+1),y(2*n+1,m+1)
!
!  Apply the Laplacian to a user-supplied function specified via its values
!  at the tensor product discretization grid.
!
!  Input parameters:
!    R1 - the inner radius of the annulus
!    R2 - the outer radius of the annulus
!    n - the parameter controlling the number of discretization nodes in the t
!      variable
!    m - the parameter controlling the number of discretizatoon nodes in the r
!      variables
!    ts - the discretization nodes in the t variable as returned by helmop_init
!    rs - the discretization nodes in the r variable as returned by helmop_init
!
!    x - an array dimensions (2n+1,m+1) whose (i,j) entry gives the value
!      of the input function at the point ( ts(i), rs(j) )
!
!  Ouput parameters:
!
!    y - an array dimensions (2n+1,m+1) whose (i,j) entry gives the value
!      of the output function at the point ( ts(i), rs(j) )
!

double complex, allocatable :: y1(:,:),y2(:,:),y3(:,:)
double complex, allocatable :: xx(:),yy(:)


allocate(y1(2*n+1,m+1))
allocate(y2(2*n+1,m+1))
allocate(y3(2*n+1,m+1))
allocate(xx(m+1),yy(m+1))


!
!  Apply 1/r^2 (d/dt)^2 and place the result in y1
!

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,r)
!$OMP DO 
do j=1,m+1
r = rs(j)
call helmop_apply_ddt(R1,R2,n,ts,x(:,j),y1(:,j))
y1(:,j) = y1(:,j) / r**2
end do
!$OMP END DO 
!$OMP END PARALLEL

!
!  Apply 1/r (d/dr) and place the result in y2
!


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,xx,yy)
!$OMP DO
do i=1,2*n+1
xx = x(i,:)
call helmop_apply_dr(R1,R2,m,rs,xx,yy)
y2(i,:) = yy
end do
!$OMP END DO
!$OMP END PARALLEL



!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,r)
!$OMP DO
do j=1,m+1
r       = rs(j)
y2(:,j) = y2(:,j)/r
end do
!$OMP END DO
!$OMP END PARALLEL


!
!  Apply (d/dr)^2 and place the result in y3
!

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,xx,yy)
!$OMP DO
do i=1,2*n+1
xx = x(i,:)
call helmop_apply_dr(R1,R2,m,rs,xx,yy)
call helmop_apply_dr(R1,R2,m,rs,yy,xx)
y3(i,:) = xx
end do
!$OMP END DO
!$OMP END PARALLEL

!
!  Sum the results of applying the three operators
!

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j)
!$OMP DO
do j=1,m+1
y(:,j) = y1(:,j) + y2(:,j) + y3(:,j)
end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine


subroutine helmop_apply_ddt(R1,R2,n,ts,x,y)
implicit double precision (a-h,o-z)
double precision       :: ts(:)
double complex         :: x(:),y(:)
!
!  Apply the discretization of the operator (d/dt)^2 to a periodic input function
!  specified via its values at the (2n+1) discretization points in the t
!  variable.
!
!  Input parameters:
!    R1 - the inner radius of the annulus
!    R2 - the outer radius of the annulus
!    n - the parameter controlling the number of discretization nodes in the
!      t variable
!    ts - the discretization nodes in the t variable
!    x - the input vector of length 2n+1 which gives the values of
!      the input function at the discretization nodes in the 
!      variable
!
!  Output parameters:
!    y - the output vector of length 2n+1 giving the values of
!      the second derivative of the input function
!
!

double complex, allocatable       :: wsave(:), z(:)
!double precision, allocatable     :: wsave(:)

allocate(wsave(15*n+1000),z(-n:n))

y = x
call zffti(2*n+1,wsave)
call zfftf(2*n+1,y,wsave)

! reorder and scale the coefficients
y    = y*1.0d0/(2*n+1.0d0)


do i=1,n
z(i-1) = y(i)
end do

do i=-1,-n,-1
z(i) = y(2*n+2+i)
end do



! apply the differential operator
do i=-n,n
z(i) = -z(i) * (i**2)
end do


! reorder the coefficients

y = 0

do i=1,n
y(i) = z(i-1) 
end do

do i=-1,-n,-1
 y(2*n+2+i) = z(i)
end do



! perform the inverse fourier transform
call zfftb(2*n+1,y,wsave)


deallocate(wsave)

end subroutine



subroutine helmop_apply_dr(R1,R2,m,rs,x,y)
implicit double precision (a-h,o-z)
double precision :: rs(m+1)
double complex   :: x(m+1),y(m+1)

!
!  Apply the discretization of the operator d/dr to an input function
!  specified via its values at the (2n+1) discretization points in the r
!  variable.
!
!  Input parameters:
!    R1 - the inner radius of the annulus
!    R2 - the outer radius of the annulus
!    m - the parameter controlling the number of discretization nodes in the
!      r variable
!    rs - the discretization nodes in the r variable
!
!    x - the input vector of length m+1 which gives the values of
!      the input function at the discretization nodes in the r
!      variable
!
!  Output parameters:
!    y - the output vector of length m+1 giving the values of
!      the derivative of the input function at the discretization nodes
!

double precision, allocatable :: xr(:),yr(:)
double precision, allocatable :: xi(:),yi(:)
double complex                :: ima

ima = (0.0d0,1.0d0)

allocate(xr(m+1),yr(m+1))
allocate(xi(m+1),yi(m+1))

xr = real(x)
xi = imag(x)

call helmop_apply_dr0(R1,R2,m,rs,xr,yr)
call helmop_apply_dr0(R1,R2,m,rs,xi,yi)

y = (yr + ima * yi)

end subroutine



subroutine helmop_apply_dr0(R1,R2,m,rs,x,y)
implicit double precision (a-h,o-z)
double precision :: x(:),y(:),rs(:)

!
!  Apply the discretization of the operator d/dr to an input function
!  specified via its values at the (2n+1) discretization points in the r
!  variable.
!
!  Input parameters:
!    m - the parameter controlling the number of discretization nodes in the
!      r variable
!    rs - the discretization nodes in the r variable
!
!    x - the input vector of length m+1 which gives the values of
!      the input function at the discretization nodes in the r
!      variable
!
!  Output parameters:
!    y - the output vector of length m+1 giving the values of
!      the derivative of the input function at the discretization nodes
!
double precision, allocatable :: wsave2(:),wsave3(:),ts0(:)
double precision, allocatable :: coefs(:),coefs2(:)

pi = acos(-1.0d0)

allocate(wsave2(6*m+100))
allocate(wsave3(6*m+100))
allocate(ts0(m+1))
allocate(coefs(m+1),coefs2(m-1))

a = R1
b = R2

call dcosti(m+1,wsave2)
call dsinti(m-1,wsave3)

do i=1,m+1
ts0(i) = (i-1)*pi/(m+0.0d0)
end do



coefs = x
call dcost(m+1,coefs,wsave2)

coefs    = coefs * 1.0d0/(m+0.0d0)
coefs(1) = coefs(1)/2


do j=1,m-1
coefs2(j) = -(j) * coefs(j+1)
end do


call dsint(m-1,coefs2,wsave3)
coefs2    = coefs2 / 2.0d0

y(2:m)   = -coefs2 / sin(ts0(2:m))

y(1)       = 0
y(m+1)     = 0

dsign     = -1

do j=0,m
dd        = j**2*coefs(j+1)
y(1)      = y(1)   + dd
y(m+1)    = y(m+1) + dsign*dd
dsign     = -dsign
end do

y = y * 2/(b-a)

y = -y
end subroutine




end module
