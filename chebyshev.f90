!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for representing functions via piecewise Chebyshev 
!  expansions.  More accurately, functions on an interval [a,b] can be represented either 
!  via their values  at the Chebyshev grids on a collection of subintervals [a,b] or
!  via their Chebyshev expansion coefficients on said subintervals.
!
!  The following routines should be regarded as public:
!    
!    chebexps - construct the n-point Clenshaw-Curtis quadrature on the interval [-1,1],
!      the matrix which takes the values of a function to the Chebyshev coefficients of the 
!      polynomial interpolating the function at the nodes of the quadrature rule, as well
!      as the "left" and "right" spectral integration matrices.
!
!    chebs - evaluate the Chebyshev polynomials of orders 0 through n at a specified
!      point in the interval [-1,1]
!
!    chebeval - evaluate a polynomial of degree n-1 whose values are given 
!      on the n-point Chebyshev grid on an interval [a,b] using the well-known 
!      barycentric interpolation formula
!
!    chebeval2 - use the Clenshaw algorithm in order to evaluate an n-term Chebyshev
!      expansion on the interval [a,b] at a specified point t in [a,b]
!
!    chebevalder - evaluate an n-term Chebyshev expansion on the interval [a,b]
!      and its derivative at a specified point x in [a,b]
!
!    chebpw_eval - evaluate a function specified by its values at k-point Chebyshev
!      grids on a collection of subintervals of [a,b] at an arbitrary point x in
!      [a,b]
!
!    chebpw_eval2 - evaluate a piecewise Chebyshev expansion given on a collection
!      of subintervals of [a,b] at a specified point x in [a,b]
!
!    chebpw_merge - merge a list of intervals into a second list of intervals
!
!    chebpw_aint - apply the matrix interpolating a function from a piecewise
!      Chebyshev grid to a user-specified collection of points to a vector
!
!    chebpw_aintt - apply the transpose of the matrix interpolating a function from 
!      a piecewise Chebyshev grid to a user-specified collection of points to a vector
!
!    chebpw_aint_c - apply the matrix interpolating a function from a piecewise
!      Chebyshev grid to a user-specified collection of points to a complex-valued
!      vector
!
!    chebpw_aintt_c - apply the transpose of the matrix interpolating a function from 
!      a piecewise Chebyshev grid to a user-specified collection of points to a complex-
!      valued vector
!
!    chebpw_plot - produce a PDF file showing the plot of a function represented
!      via its values on a  piecewise Chebyshev grid
!
!    chebroots - find the real roots of a Chebyshev expansion on the interval (a,b)
!      using the companion matrix method 
!
!       WARNING -- THIS ROUTINE RELIES ON LAPACK ROUTINES
!
!    chebpw_roots - find the real roots of a piecewise Chebyshev expansion using
!      the companion matrix method 
!
!       WARNING THIS ROUTINE RELIES ON LAPACK ROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module chebyshev

use utils
use iso_c_binding


double precision, private, parameter         :: dtail = 0.25d0

interface

subroutine chebadapfun(t,val,userptr)
import c_ptr
double precision, intent(in)      :: t
double precision, intent(out)     :: val
type(c_ptr)                       :: userptr
end subroutine



subroutine chebadapfun2(t,val1,val2,userptr)
import c_ptr
double precision, intent(in)      :: t
double precision, intent(out)     :: val1,val2
type(c_ptr)                       :: userptr
end subroutine


end interface

!
!  The chebexps_data structure contains the following elements:
!
!    xs - an array containing the n quadrature nodes
!    whts - an array containing the n quadrature weights
!    u - the (n,n) matrix which takes the values of an n-term Chebyshev expansion
!      at the n quadrature nodes to the n expansion coefficients
!    v - the (n,n) matrix which takes the coefficients of an nterm-Chebyshev
!      expansion to its values at the n quadratue nodes
!    aintl - the "left" spectral integration matrix which takes the values
!      of a function f(t) on the Chebyshev nodes to the value of the function g(t) 
!      defined via the formula
!
!                     t
!          g(t) = \int  f(u) du
!                     a
!
!    aintr - the "right" spectral integration matrix which takes the values
!      of a function f(t) on the Chebyshev nodes to the value of the function g(t) 
!      defined via the formula
!
!                     t
!          g(t) = \int_ f(u) du
!                     b
!    aintr2, aintr3, aintl2, aintl3 - the higher order left and right spectral 
!      integration matrices
!

type      chebexps_data
integer                       :: k
double precision, allocatable :: xs(:),whts(:),u(:,:),v(:,:)

! left and right spectral integration matrices 
double precision, allocatable :: aintr(:,:),aintr2(:,:),aintr3(:,:)
double precision, allocatable :: aintl(:,:),aintl2(:,:),aintl3(:,:)

! spectral differentiation matrix
double precision, allocatable :: adiff(:,:)
end type  chebexps_data


contains



subroutine chebexps(n,chebdata)
implicit double precision (a-h,o-z)
integer                               :: n
type(chebexps_data), intent(out)      :: chebdata

!
!  Populate a chebexps structure.
!
!  Input parameters:
!    n - the number of points in the chebyshev grid
!
!  Output parameters:
!   chebdata - the chebexps_data structure containing the elements described
!    above
!   
!

double precision, allocatable :: pols(:),c(:,:),d(:,:),xx(:,:)
data pi /3.14159265358979323846264338327950288d0/

allocate(chebdata%xs(n),chebdata%whts(n),chebdata%u(n,n),chebdata%v(n,n))
allocate(pols(n+1),c(1,n),d(1,n))

chebdata%k = n
!
!  Construct the nodes
!

h = pi/(n-1)
do i=1,n
chebdata%xs(n-i+1) = cos(h*(i-1))
end do

!
!  Construct the matrix u which takes values to coefficients
!

do i=1,n
x = chebdata%xs(i)
call chebs(x,n-1,pols)
do j=1,n
chebdata%u(j,i) = pols(j)
chebdata%v(i,j) = pols(j)
end do
end do


chebdata%u(1,:) = chebdata%u(1,:)/2
chebdata%u(n,:) = chebdata%u(n,:)/2
chebdata%u(:,1) = chebdata%u(:,1)/2
chebdata%u(:,n) = chebdata%u(:,n)/2
chebdata%u = chebdata%u*2.0d0/(n-1)

!
!  Construct the weights by multiplying u^t on the left by the
!  integrals of the Chebyshev polynomials.
!

c=0
c(1,1) = 2.0d0
do i=2,n-1,2
c(1,i+1) = 1.0d0/(i+1)-1.0d0/(i-1)
end do

d = matmul(c,chebdata%u)
chebdata%whts = d(1:n,1)

!
!  Form the matrix which takes the values of a function f(t) to the values of
!
!              t
!    g(t) =  \int     f(u) du
!              a
!  

allocate(xx(n,n),chebdata%aintr(n,n),chebdata%aintl(n,n))

do i=1,n
call chebs(chebdata%xs(i),n,pols)
xx(i,1) = chebdata%xs(i)
xx(i,2) = chebdata%xs(i)**2/2.0d0
do j=3,n
xx(i,j) = 0.5d0 * (pols(j+1)/j-pols(j-1)/(j-2))
end do
end do

do i=2,n
xx(i,:) = xx(i,:) - xx(1,:)
end do
xx(1,:) = 0

chebdata%aintl = matmul(xx,chebdata%u)

!
!  Form the matrix which takes the values of a function f(t) to the values of
!
!              t
!    g(t) =  \int     f(u) du
!              b
!  

xx = 0

do i=1,n
call chebs(chebdata%xs(i),n,pols)
xx(i,1) = chebdata%xs(i)
xx(i,2) = chebdata%xs(i)**2/2.0d0
do j=3,n
xx(i,j) = 0.5d0 * (pols(j+1)/j-pols(j-1)/(j-2))
end do
end do

do i=1,n-1
xx(i,:) = xx(i,:) - xx(n,:)
end do
xx(n,:) = 0

chebdata%aintr = matmul(xx,chebdata%u)

!
!  Form the higher order spectral integration matrices
!

allocate(chebdata%aintr2(n,n),chebdata%aintr3(n,n))
chebdata%aintr2 = matmul(chebdata%aintr,chebdata%aintr)
chebdata%aintr3 = matmul(chebdata%aintr,chebdata%aintr2)

allocate(chebdata%aintl2(n,n),chebdata%aintl3(n,n))
chebdata%aintl2 = matmul(chebdata%aintl,chebdata%aintl)
chebdata%aintl3 = matmul(chebdata%aintl,chebdata%aintl2)

!
!  Form the spectral differentiation matrix
!

allocate(chebdata%adiff(n,n))

do j=1,n
xx(1,j) = (-1)**(j) * (j-1)**2
xx(n,j) = (j-1)**2
end do

do i=2,n-1
x = chebdata%xs(i)
call chebs(x,n,pols)
do j=1,n
xx(i,j) = (j-1) * x * pols(j) - (j-1) * pols(j+1) 
xx(i,j) = xx(i,j) / (1-x**2)
end do
end do

chebdata%adiff = matmul(xx,chebdata%u)

end subroutine



subroutine chebs(x,n,pols)
implicit double precision (a-h,o-z)

integer          :: n 
double precision :: pols(0:n),x

!
!  Evaluate the Chebyshev polynomials of degree 0 through n at a specified point
!  using the standard 3-term recurrence relation.
!
!  Input parameters:
!
!    x - point at which to evaluate the polynomials
!    n - the order of the polynomials to evaluate
!
!  Output parameters:
!
!    pols - this user-supplied and allocated array of length n+1 will
!      contain the values of the polynomials of order 0 through n upon
!      return
!

if (x .eq. 1.0d0) then
do i=0,n
pols(i) = 1.0d0
end do
return
endif

if (x .eq. -1.0d0) then
pols(0) = 1.0d0
do i=1,n
pols(i) = -pols(i-1)
end do
return
endif

pols(0) = 1.0d0
if (n .eq. 0) return

pols(1) = x
if (n .eq. 1) return

xx1 = 1.0d0
xx2 = x

do i=1,n-1
xx        = 2*x*xx2 - xx1
pols(i+1) = xx
xx1       = xx2
xx2       = xx
end do

end subroutine



subroutine chebeval(a,b,n,xs,vals,x,val)
implicit double precision (a-h,o-z)

integer :: n
double precision ::  xs(n),vals(n),x,val

!
!  Use the barycentric formula to evaluate a function given its values at the 
!  n-point Chebyshev grid on an interval [a,b].
!
!  Input parameters:
!
!    (a,b) - the interval on which the function is given
!    n - the number of nodes in the Chebyshev grid
!    xs - an array specifying the n Chevyshev node on the interval [-1,1]
!    vals - the values of the function on the n Chebyshev nodes on the
!      interval [-1,1]
!    x - the point in the interval (a,b) at which the function is to be
!      evaluated
!     
!  Output parameters:
!
!   val - the approximate value of the function at the point x
!

eps0 = epsilon(0.0d0)

xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0

dd1 = 1.0d0

do i=1,n
dd=1.0d0
if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0

diff = xx-xs(i)

!
!  Handle the case in which the target node coincide with one of
!  of the Chebyshev nodes.
!

if(abs(diff) .le. eps0) then
val = vals(i)
return
endif

!
!  Otherwise, construct the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i)
sum2 = sum2+dd
dd   = - dd
end do

val = sum1/sum2


end subroutine


subroutine chebeval2(a,b,n,coefs,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)             :: n
double precision, intent(in)    :: a,b,x,coefs(n)
double precision, intent(out)   :: val

!
!  Use the Clenshaw algorithm in order to evaluate a Chebyshev expansion on the 
!  interval [a,b].
!
!  Input parameters:
!    (a,b) - the interval on which the expansion is given
!    n - the number of terms in the Chebyshev expansion
!    coefs - an array of length n specifying the expansion coefficients
!    x - the point at which to evaluate the expansion
!
!  Output parameters:
!    val - the value of the expansion at the point x
!

xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)

b2 = coefs(n)
b1 = 2*xx*b2+coefs(n-1)

do i=n-2,2,-1
b0  = coefs(i)+2*xx*b1-b2
b2 = b1
b1 = b0
end do

val = b1 * xx + (coefs(1)-b2)

end subroutine




subroutine chebpw_eval(nints,ab,k,xscheb,vals,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: xscheb(k),ab(2,nints),vals(k,nints)
double precision, intent(out) :: val

!
!  Evaluate a function represented via its values at the nodes of the k-point
!  Chebyshev grids on a collection of subintervals of [a,b].
!
!  Input parameters:
!
!    (nints,ab) - arrays specifying the collection of subintervals of [a,b]
!    k - the number of terms in the Chebyshev expansions
!    xscheb - the nodes of the k-point Clenshaw-Curtis quadrature on [-1,1]
!    vals - a (k,nints) array the jth column of which gives the values of the
!      function at the nodes of the k-point Chebyshev grid in the jth
!      subinterval
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the function at the point x
! 


!
!  Conduct several iterations of a binary search for the interval.
!

eps0 = epsilon(0.0d0)


!
!  Conduct a brute force check from here.
!

do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do
a = ab(1,int)
b = ab(2,int)


! call chebeval(a,b,k,xscheb,vals(1,int),x,val)
! return
xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0

dd1 = 1.0d0

do i=1,k
dd=1.0d0
if (i .eq. 1 .OR. i .eq. k) dd = 0.5d0

diff = xx-xscheb(i)

!
!  Handle the case in which the target node coincide with one of
!  of the Chebyshev nodes.
!

if(abs(diff) .le. eps0 ) then
val = vals(i,int)
return
endif

!
!  Otherwise, construct the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i,int)
sum2 = sum2+dd
dd   = - dd
end do

val = sum1/sum2


end subroutine





subroutine chebpw_eval_pair(nints,ab,k,xscheb,vals1,vals2,x,val1,val2)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: xscheb(k),ab(2,nints),vals1(k,nints),vals2(k,nints)
double precision, intent(out) :: val1,val2

!
!  Evaluate two functions represented via its values at the nodes of the k-point
!  Chebyshev grids on a collection of subintervals of [a,b].
!
!  Input parameters:
!
!    (nints,ab) - arrays specifying the collection of subintervals of [a,b]
!    k - the number of terms in the Chebyshev expansions
!    xscheb - the nodes of the k-point Clenshaw-Curtis quadrature on [-1,1]
!    vals - a (k,nints) array the jth column of which gives the values of the
!      function at the nodes of the k-point Chebyshev grid in the jth
!      subinterval
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the function at the point x
! 


!
!  Conduct several iterations of a binary search for the interval.
!

eps0 = epsilon(0.0d0)


do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do
a = ab(1,int)
b = ab(2,int)


xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0
sum3=0


dd1 = 1.0d0

do i=1,k
dd=1.0d0
if (i .eq. 1 .OR. i .eq. k) dd = 0.5d0

diff = xx-xscheb(i)

if(abs(diff) .le. eps0 ) then
val1 = vals1(i,int)
val2 = vals2(i,int)
return
endif


dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals1(i,int)
sum2 = sum2+dd
sum3 = sum3+dd*vals2(i,int)
dd   = - dd
end do

val1 = sum1/sum2
val2 = sum3/sum2

end subroutine


subroutine chebpw_eval2(nints,ab,k,coefs,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: ab(2,nints),coefs(k,nints)
double precision, intent(out) :: val

!
!  Evaluate a function represented via piecewise chebyshev expansions on
!  a collection of subintervals.
!
!  Input parameters:
!    (nints,ab) - the 
!    k - an integer specifying the order of the Chebyshev expansions; on
!    
!    coefs - a (k,nints) array whose jth column specified the coefficients
!     of the function's Chebyshev expansion on the jth subinterval
!    x - the point at which to evaluate
!
!
!  Output parameters:
!    val - the value of the function at the point x
! 


double precision :: pols(k)
do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do
a = ab(1,int)
b = ab(2,int)


!
!  Evaluate the Chebyshev expansion
!


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
call chebs(xx,k-1,pols)

val = 0
do i=1,k
val = val + coefs(i,int)*pols(i)
end do


return

! call chebeval2(a,b,k,coefs(1,int),x,val)


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
xx2 = 2*xx

b2 = coefs(k,int)
b1 = xx2*b2+coefs(k-1,int)

do i=k-2,2,-1
b0  = coefs(i,int)+xx2*b1-b2
b2 = b1
b1 = b0
end do

val = b1 * xx + (coefs(1,int)-b2)

return
end subroutine



subroutine chebpw_ceval2(nints,ab,k,coefs,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: ab(2,nints)
double complex, intent(in)    :: coefs(k,nints)
double complex, intent(out)   :: val

!
!  Evaluate a function represented via piecewise chebyshev expansions on
!  a collection of subintervals.
!
!  Input parameters:
!    (nints,ab) - the 
!    k - an integer specifying the order of the Chebyshev expansions; on
!    
!    coefs - a (k,nints) array whose jth column specified the coefficients
!     of the function's Chebyshev expansion on the jth subinterval
!    x - the point at which to evaluate
!
!
!  Output parameters:
!    val - the value of the function at the point x
! 


double precision :: pols(k)

do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do
a = ab(1,int)
b = ab(2,int)


!
!  Evaluate the Chebyshev expansion
!


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
call chebs(xx,k-1,pols)

val = 0
do i=1,k
val = val + coefs(i,int)*pols(i)
end do


return
end subroutine

subroutine chebpw_eval22(nints,ab,k,coefs1,coefs2,x,val1,val2)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: ab(2,nints),coefs1(k,nints),coefs2(k,nints)
double precision, intent(out) :: val1,val2

!
!  Evaluate a function represented via piecewise chebyshev expansions on
!  a collection of subintervals.
!
!  Input parameters:
!    (nints,ab) - the 
!    k - an integer specifying the order of the Chebyshev expansions; on
!    
!    coefs - a (k,nints) array whose jth column specified the coefficients
!     of the function's Chebyshev expansion on the jth subinterval
!    x - the point at which to evaluate
!
!
!  Output parameters:
!    val - the value of the function at the point x
! 


double precision :: pols(k)



do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do
a = ab(1,int)
b = ab(2,int)

!
!  Evaluate the Chebyshev expansion
!

xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
call chebs(xx,k-1,pols)

val1 = 0
val2 = 0
do i=1,k
val1 = val1 + coefs1(i,int)*pols(i)
val2 = val2 + coefs2(i,int)*pols(i)
end do

return
end subroutine


subroutine chebpw_eval23(nints,ab,k,coefs1,coefs2,coefs3,x,val1,val2,val3)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: ab(2,nints),coefs1(k,nints),coefs2(k,nints)
double precision, intent(in)  :: coefs3(k,nints)
double precision, intent(out) :: val1,val2,val3

!
!  Evaluate a function represented via piecewise chebyshev expansions on
!  a collection of subintervals.
!
!  Input parameters:
!    (nints,ab) - the 
!    k - an integer specifying the order of the Chebyshev expansions; on
!    
!    coefs - a (k,nints) array whose jth column specified the coefficients
!     of the function's Chebyshev expansion on the jth subinterval
!    x - the point at which to evaluate
!
!
!  Output parameters:
!    val - the value of the function at the point x
! 


double precision :: pols(k)

do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do
a = ab(1,int)
b = ab(2,int)


!
!  Evaluate the Chebyshev expansion
!


xx = (x - (b+a)/2.0d0) * 2.0d0/(b-a)
call chebs(xx,k-1,pols)

val1 = 0
val2 = 0
val3 = 0
do i=1,k
val1 = val1 + coefs1(i,int)*pols(i)
val2 = val2 + coefs2(i,int)*pols(i)
val3 = val3 + coefs3(i,int)*pols(i)

end do

return
end subroutine


subroutine chebpw_aint(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints), x(chebdata%k,nints)
double precision           :: ts(nts)
double precision           :: y(nts)

!
!  Apply the matrix interpolating a function whose values are given at the
!  nodes of a piecewise Chebyshev discretization scheme to a set of user-specified 
!  points.
!
!  IMPORTANT - the list of user-specified points must be sorted in ascending order.
!
!  Input parameters:
!    chebdata - chebyshev expansion and interpolation data
!    nints - the number of intervals in the discretization scheme
!    ab - a (2,nints) array giving the endpoints of the intervals in the discretization
!      scheme
!    nts - the number of user-specified nodes
!    ts - an array containing the user specified nodes SORTED IN ASCENDING ORDER
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!


double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k

do i=1,nts
t = ts(i)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int

xx    = (2*t - (b+a) ) /(b-a)
dsign = 1.0d0
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
val = x(l,int)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2
val     = sum(whts*x(:,int))/sum(whts)

1000 continue

y(i) = val
end do


end subroutine


subroutine chebpw_aintt(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double precision           :: x(nts),y(chebdata%k,nints)

!
!  Apply the transpose of the matrix interpolating a function whose values
!  are given at the nodes of a piecewise Chebyshev discretization scheme to
!  a set of user-specified points.
!
!  IMPORTANT - the list of user-specified points must be sorted in ascending order.
!
!  Input parameters:
!    chebdata - chebyshev expansion and interpolation data
!    nints - the number of intervals in the discretization scheme
!    ab - a (2,nints) array giving the endpoints of the intervals in the discretization
!      scheme
!    nts - the number of user-specified nodes
!    ts - an array containing the user specified nodes SORTED IN ASCENDING ORDER
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!

double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k
y    = 0

!
!  Traverse the list of ts
!

do i=1,nts
t = ts(i)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int

xx    = (2*t - (b+a) ) /(b-a)
dsign = 1.0d0
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
y(l,int) = y(l,int) + x(i)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do
whts(1) = whts(1)/2
whts(k) = whts(k)/2

y(:,int) = y(:,int) + x(i)*whts(:)/sum(whts)
1000 continue

end do

end subroutine


subroutine chebpw_aint_c(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double complex             :: y(nts), x(chebdata%k,nints)


!
!  Apply the matrix interpolating a function whose values are given at the
!  nodes of a piecewise Chebyshev discretization scheme to a set of user-specified 
!  points.
!
!  IMPORTANT - the list of user-specified points must be sorted in ascending order.
!
!  Input parameters:
!    chebdata - chebyshev expansion and interpolation data
!    nints - the number of intervals in the discretization scheme
!    ab - a (2,nints) array giving the endpoints of the intervals in the discretization
!      scheme
!    nts - the number of user-specified nodes
!    ts - an array containing the user specified nodes SORTED IN ASCENDING ORDER
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!


double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k

do i=1,nts
t = ts(i)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int

xx    = (2*t - (b+a) ) /(b-a)
dsign = 1.0d0
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
val = x(l,int)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(k) = whts(k)/2
val     = sum(whts*x(:,int))/sum(whts)

1000 continue

y(i) = val
end do


end subroutine


subroutine chebpw_aintt_c(chebdata,nints,ab,nts,ts,x,y)
implicit double precision (a-h,o-z)

type(chebexps_data)        :: chebdata
integer                    :: nints
double precision           :: ab(2,nints)
double precision           :: ts(nts)
double complex             :: x(nts),y(chebdata%k,nints)

!
!  Apply the transpose of the matrix interpolating a function whose values
!  are given at the nodes of a piecewise Chebyshev discretization scheme to
!  a set of user-specified points.
!
!  IMPORTANT - the list of user-specified points must be sorted in ascending order.
!
!  Input parameters:
!    chebdata - chebyshev expansion and interpolation data
!    nints - the number of intervals in the discretization scheme
!    ab - a (2,nints) array giving the endpoints of the intervals in the discretization
!      scheme
!    nts - the number of user-specified nodes
!    ts - an array containing the user specified nodes SORTED IN ASCENDING ORDER
!    x - the vector to which the transpose of the matrix should be applied
!
!  Output parameters:
!    y - the output vector
!

double precision           :: whts(chebdata%k)

eps0 = epsilon(0.0d0)
int0 = 1
k    = chebdata%k
y    = 0

!
!  Traverse the list of ts
!

do i=1,nts
t = ts(i)

do int=int0,nints-1
if (t .lt. ab(2,int)) exit
end do
a    = ab(1,int)
b    = ab(2,int)
int0 = int

xx    = (2*t - (b+a) ) /(b-a)
dsign = 1.0d0
do l=1,k
diff    = chebdata%xs(l) - xx
if (abs(diff) .lt. eps0) then
y(l,int) = y(l,int) + x(i)
goto 1000
endif
whts(l) = dsign/diff
dsign   = -dsign
end do
whts(1) = whts(1)/2
whts(k) = whts(k)/2

y(:,int) = y(:,int) + x(i)*whts(:)/sum(whts)
1000 continue

end do

end subroutine


subroutine chebpts(n,xs)
implicit double precision (a-h,o-z)

integer, intent(in)                          :: n
double precision, allocatable, intent(out)   :: xs(:)

!
!  Return the nodes of the n-point Chebyshev grid on the interval [-1,1].
!
!  Input parameters:
!    n  - the number of points in the grid
!
!  Output parameters:
!    xs - an array of length n containing the nodes
!

data pi /3.14159265358979323846264338327950288d0/

allocate(xs(n))

if (n .eq. 1) then
xs(1) = 0.0d0
return
endif

h = pi/(n-1)
do i=1,n
xs(n-i+1) = cos(h*(i-1))
end do

end subroutine



subroutine chebadap(ier,eps,a,b,fun,k,chebdata,nints,ab,userptr)
implicit double precision (a-h,o-z)

type(chebexps_data)                        :: chebdata
double precision, intent(in)               :: eps,a,b
integer, intent(in)                        :: k
integer, intent(out)                       :: nints
double precision, intent(out), allocatable :: ab(:,:)
procedure(chebadapfun)                     :: fun
type(c_ptr)                                :: userptr

!
!  This subroutine adaptively discretizes a user-specified function.  More 
!  specifically, it recursively subdivides the specified interval [a,b] until
!  it is represented via k-term polynomial expansions on each subinterval.
!
!  Input parameters:
!
!    eps - the precision with which the function should be represented
!    (a,b) - the interval on which the function is given
!    fun - an external procedure conforming to the interface chebadapfun
!      which returns the values of the function
!
!    k - the number of terms in the polynomial expansions
!    xs - the nodes of the k-point Clenshaw-Curtis quadrature returned by chebexps
!    u - the (k,k) matrix u which takes values to coefficients (as returned by
!      chebexps)
!   
!  Output parameters:
!
!    ier - an error return code;
!       ier =   0     indicates successful execution
!       ier =   4     means that the present maximum number of intervals
!                     was exceeded before an appropriate discretization scheme
!                     could be obtained
!       ier =   1024  means that the maximum number of recursion steps was
!                     exceeded 
!
!    nints - the number of subintervals in the discretization scheme
!    ab - an (2,nints) array each column of which specifies one of the
!      subintervals in the discretization scheme
!

double precision, allocatable :: ab0(:,:),ab1(:,:),vals(:),coefs(:),ts(:)

call elapsed(t1)

ier     = 0
maxints = 10000
ntail   = k * dtail
!timemax = 100

if (a .eq. b) then
nints = 0
allocate(ab(2,nints))
return
endif


allocate(ab0(2,maxints),ab1(2,maxints),vals(k),coefs(k),ts(k))

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b

nints1   = 0

do while (nints0 .gt. 0) 

call elapsed(t2)
!time = t2-t1
!print *,"chebadap: ",time,timemax

! if (time .gt. timemax) then
! call prina("in chebadap, maximum time exceeded")
! call prin2("eps = ",eps)
! call prin2("a   = ",a)
! call prin2("b   = ",b)
! call prini("nints1 = ",nints1)
! call prin2("ab1 = ",ab1(:,1:nints1))
! call prini("nints0 = ",nints0)
! call prin2("ab0 = ",ab0(:,1:nints0))

! stop
! endif

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
c0 = (a0+b0)/2

if (b0-a0 .eq. 0) then
ier = 1024
return
endif

nints0 = nints0-1

!
!  Evaluate the function at each of the Chebyshev nodes on [a0,b0]
!

ts = chebdata%xs*(b0-a0)/2 + (b0+a0)/2

do i=1,k
call fun(ts(i),vals(i),userptr)
end do

coefs = matmul(chebdata%u,vals)


!
!  Measure the relative error in the trailing coefficients
!


dd2 = maxval(abs(coefs))
dd1 = maxval(abs(coefs(k-ntail+1:k)))


ifsplit = 0
if (dd1 .gt. eps*dd2 ) ifsplit=1


if(ifsplit .eq. 1) then

if (nints0+2 .gt. maxints) then
ier = 4
return
endif

c0 = (a0+b0)/2

nints0 = nints0 + 1
ab0(1,nints0) = c0
ab0(2,nints0) = b0
nints0 = nints0 + 1
ab0(1,nints0) = a0
ab0(2,nints0) = c0

else

!print *,a0,b0,dd1/dd2
!call prin2("coefs = ",coefs)

if (nints1+1 .gt. maxints) then
ier = 4
return
endif

nints1 = nints1 + 1
ab1(1,nints1) = a0
ab1(2,nints1) = b0
endif

end do


nints = nints1
allocate(ab(2,nints))
ab = ab1(:,1:nints)


end subroutine




subroutine chebadap2(ier,eps,nintsin,abin,fun,k,chebdata,nints,ab,userptr)
implicit double precision (a-h,o-z)

type(chebexps_data)                        :: chebdata
double precision, intent(in)               :: eps
integer, intent(in)                        :: k
integer, intent(out)                       :: nints
double precision                           :: abin(2,nintsin)
double precision, intent(out), allocatable :: ab(:,:)
procedure(chebadapfun2)                    :: fun
type(c_ptr)                                :: userptr

!
!  This subroutine adaptively discretizes a user-specified function.  More 
!  specifically, it recursively subdivides the specified interval [a,b] until
!  it is represented via k-term polynomial expansions on each subinterval.
!
!  Input parameters:
!
!    eps - the precision with which the function should be represented
!    (a,b) - the interval on which the function is given
!    fun - an external procedure conforming to the interface chebadapfun
!      which returns the values of the function
!
!    k - the number of terms in the polynomial expansions
!    xs - the nodes of the k-point Clenshaw-Curtis quadrature returned by chebexps
!    u - the (k,k) matrix u which takes values to coefficients (as returned by
!      chebexps)
!   
!  Output parameters:
!
!    ier - an error return code;
!       ier =   0     indicates successful execution
!       ier =   4     means that the present maximum number of intervals
!                     was exceeded before an appropriate discretization scheme
!                     could be obtained
!       ier =   1024  means that the maximum number of recursion steps was
!                     exceeded 
!
!    nints - the number of subintervals in the discretization scheme
!    ab - an (2,nints) array each column of which specifies one of the
!      subintervals in the discretization scheme
!

integer                       :: OMP_GET_MAX_THREADS
double precision, allocatable :: ab0(:,:),ab1(:,:),vals(:),coefs(:),ts(:),vals2(:)

ier     = 0
maxints = 1000000
ntail   = k * dtail


! if (nthreads .eq. 0) then
! nthreads = OMP_GET_MAX_THREADS()
! endif

allocate(ab0(2,maxints),ab1(2,maxints),vals(k),coefs(k),ts(k),vals2(k))

nints0            = nintsin
ab0(:,nints0:1:-1)   = abin
nints1            = 0

do while (nints0 .gt. 0) 

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
c0 = (a0+b0)/2

if (b0-a0 .eq. 0) then
ier = 1024
return
endif

nints0 = nints0-1

!
!  Evaluate the function at each of the Chebyshev nodes on [a0,b0]
!

ts = chebdata%xs*(b0-a0)/2 + (b0+a0)/2

do i=1,k
call fun(ts(i),vals(i),vals2(i),userptr)
end do

ifsplit = 0

!
!  Measure the relative error in the trailing coefficients
!
coefs = matmul(chebdata%u,vals)
dd2 = maxval(abs(coefs))
dd1 = maxval(abs(coefs(k-ntail+1:k)))
if (dd1 .gt. eps*dd2) ifsplit = 1

coefs = matmul(chebdata%u,vals2)
dd2 = maxval(abs(coefs))
dd1 = maxval(abs(coefs(k-ntail+1:k)))
if (dd1 .gt. eps*dd2) ifsplit = 1


if (dd1 .gt. eps*dd2 ) then

if (nints0+2 .gt. maxints) then
ier = 4
return
endif

c0 = (a0+b0)/2


nints0 = nints0 + 1
ab0(1,nints0) = c0
ab0(2,nints0) = b0

nints0 = nints0 + 1
ab0(1,nints0) = a0
ab0(2,nints0) = c0

else

if (nints1+1 .gt. maxints) then
ier = 4
return
endif

nints1 = nints1 + 1
ab1(1,nints1) = a0
ab1(2,nints1) = b0
endif

end do

nints = nints1
allocate(ab(2,nints))
ab = ab1(:,1:nints)


end subroutine


subroutine chebpw_plot(ifshow,filename,nints,ab,k,xs,vals)
implicit double precision (a-h,o-z)

integer, intent(in)          :: nints,k,ifshow
double precision, intent(in) :: ab(2,nints),xs(k),vals(k,nints)
character(len=*)             :: filename

!
!  Produce a GNUPLOT file called "gnuplot.???" where ??? is a specified
!  integer which contains commnds for plotting a function specified
!  by its values at the nodes of a discretization scheme.
!  
!  Input parameters:
!    title = an asterick-terminated character string which specifies the
!      title for the plot
!    iplot = the index of the GNUPLOT ouput file
!    (nints,ab) - the subintervals of the discretization scheme used to
!      represent the input function
!    k - the length of the Chebyshev grid on each interval
!    xs - the nodes of the k-point Chebyshev grid on the interval [-1,1]   
!    vals - a (k,nints) array specifying the values of the input function
!
!  Output parameters:
!
!    N/A  
!
!

double precision, allocatable :: ts(:,:)

allocate(ts(k,nints))

do int = 1,nints
a0 = ab(1,int)
b0 = ab(2,int)
do i   = 1,k
ts(i,int) = xs(i)*(b0-a0)/2 + (b0+a0)/2
end do
end do

n = k *nints
call  plot_function(filename,"","",n,ts,vals)

end subroutine


subroutine chebders(x,n,pols,ders)
implicit double precision (a-h,o-z)

integer :: n 
double precision :: pols(n+2),ders(n+1)

!
!  Evaluate the Chebyshev polynomials order 0 through n and their derivatives
!  at a specified point using the standard 3-term recurrence relation.
!
!  WARNING: the array of polynomials values passed by the user must be of
!  length n+2!!!
!
!  Input parameters:
!
!    x - point at which to evaluate the polynomials
!    n - the order of the polynomials to evaluate
!
!  Output parameters:
!
!    pols - this user-supplied and allocated array of length n+2 will
!      contain the values of the polynomials of order 0 through n upon
!      return
!   ders - this user-supplied array of length n+1 will contain the values of
!      the derivative of orders 0 through n upon return
!

eps0 =epsilon(0.0d0)

if (abs(x-1) .lt. eps0) then
pols = 1
do i=0,n
ders(i+1) = i**2
end do
return
endif

if (abs(x+1) .lt. eps0) then
do i=0,n/2

pols(2*i+1) = 1
ders(2*i+1) = -(2*i)**2

pols(2*i+2) = -1
ders(2*i+2) = (2*i+1)**2

end do
return
endif


pols(1) = 1.0d0
ders(1) = 0.0d0

if (n .eq. 0) return

pols(1) = 1.0d0
pols(2) = x

ders(1) = 0
ders(2) = 1

if (n .eq. 1) return

xx1 = 1.0d0
xx2 = x

do i=1,n
xx        = 2*x*xx2 - xx1
pols(i+2) = xx
xx1       = xx2
xx2       = xx
end do

!
!  Compute the derivatives
!


do i=2,n
if (abs(x-1) .lt.  eps0*10) then
ders(i+1) = 1.0d0
else if(abs(x+1) .lt. eps0*10) then
ders(i+1) = (-1.0d0)**(i-1) * i**2
else
ders(i+1) = i*(x*pols(i+1)-pols(i+2))/(1-x**2)
endif
end do


end subroutine




subroutine chebpw_ceval(nints,ab,k,xscheb,vals,x,val)
implicit double precision (a-h,o-z)

integer, intent(in)           :: nints,k
double precision, intent(in)  :: xscheb(k),ab(2,nints)
double complex                :: vals(k,nints)
double complex, intent(out)   :: val

!
!  Evaluate a function represented via its values at the nodes of the k-point
!  Chebyshev grids on a collection of subintervals of [a,b].
!
!  Input parameters:
!
!    (nints,ab) - arrays specifying the collection of subintervals of [a,b]
!    k - the number of terms in the Chebyshev expansions
!    xscheb - the nodes of the k-point Clenshaw-Curtis quadrature on [-1,1]
!    vals - a (k,nints) array the jth column of which gives the values of the
!      function at the nodes of the k-point Chebyshev grid in the jth
!      subinterval
!    x - the point in [a,b] at which to evaluate the function
!
!  Output parameters:
!    val - the value of the function at the point x
! 

double complex :: sum1

!
!  Conduct several iterations of a binary search for the interval.
!

eps0 = epsilon(0.0d0)


intl   = 1
intr   = nints
do int = intl,intr-1
b = ab(2,int)
if (x .le. b) exit
end do


!
!  Call chebeval to evaluate the expansion and the save the index
!  of the interval containing x.
!

a = ab(1,int)
b = ab(2,int)


! call chebeval(a,b,k,xscheb,vals(1,int),x,val)
! return
xx   = (2*x - (b+a) ) /(b-a)

sum1=0
sum2=0

dd1 = 1.0d0

do i=1,k
dd=1.0d0
if (i .eq. 1 .OR. i .eq. k) dd = 0.5d0

diff = xx-xscheb(i)

!
!  Handle the case in which the target node coincide with one of
!  of the Chebyshev nodes.
!

if(abs(diff) .le. eps0 ) then
val = vals(i,int)
return
endif

!
!  Otherwise, construct the sums.
!

dd   = (dd1*dd)/diff
dd1  = - dd1
sum1 = sum1+dd*vals(i,int)
sum2 = sum2+dd
dd   = - dd
end do

val = sum1/sum2


end subroutine


subroutine chebpw_merge(nints0,ab0,nints,ab)
implicit double precision (a-h,o-z)
double precision              :: ab0(:,:)
double precision, allocatable :: ab(:,:)
!
!  Merge the list of intervals ab0 into the list of intervals ab
!
!  Input parameters:
!    nints0 - the number of intervals in the list of intervals to
!     merge into the main list
!    ab0 - the list of intervals in the list to merge
!
!    nints - the number of 
!
!  Output parameters:

double precision, allocatable :: xs(:)


allocate(xs(nints+nints0+2))
xs = 0

! call prini("nints0 = ",nints0)
! call prin2("ab0    = ",ab0)

! call prini("nints  = ",nints)
! call prin2("ab     = ",ab)


!
!  Add each partition point to a list
!
nn = 0

if (nints .gt. 0) then

do int=1,nints
nn      = nn + 1
xs(nn)  = ab(1,int)
end do

nn  = nn + 1
xs(nn) = ab(2,nints)

endif


if (nints0 .gt. 0) then
do i=1,nints0
nn      = nn + 1
xs(nn) = ab0(1,i)
end do

nn  = nn + 1
xs(nn) = ab0(2,nints0)

endif

!
!  Sort the list
!
call quicksort(nn,xs)
!call prin2("xs = ",xs(1:nn))

!
!  Count the length of the new list by removing duplicates
!

mm = 0
do i=1,nn-1
if (xs(i+1) .ne. xs(i)) mm = mm + 1
end do


!
!  Form the new list of intervals 
!

deallocate(ab)
nints = mm
allocate(ab(2,nints))
mm = 0

do i=1,nn-1

if (xs(i+1) .ne. xs(i)) then
mm = mm + 1
ab(1,mm) = xs(i)
ab(2,mm) = xs(i+1)
endif

end do

!call prin2("ab = ",ab)

end subroutine



subroutine chebdiff(k,xs,adiff)
implicit double precision (a-h,o-z)
integer, intent(in)                        :: k
double precision, intent(in)               :: xs(k)
double precision, allocatable, intent(out) :: adiff(:,:)
!
!  Return the (k,k) Chebyshev spectral differentiation matrix.
!
!  Input parameters:
!    k - the size of the desired spectral differentiation matirx
!    xs - the nodes of the k-point Chebyshev grid
!
!  Output parameters:
!   adiff - a (k,k) 
!

allocate(adiff(k,k))

adiff = 0

do i=0,k-1
do j=0,k-1

if (i .eq. j) then
if ( i .eq. 0) then
val = -(2.0d0 * (k-1.0d0)**2 + 1.0d0)/6.0d0
elseif (i .eq. k-1) then
val = (2.0d0 * (k-1.0d0)**2 + 1.0d0)/6.0d0
else
val =  -0.5d0 * xs(j+1)/(1-xs(j+1)**2)
endif

else

c1 = 1
c2 = 1
if (i .eq. 0 .OR. i .eq. k-1) c1 = 2
if (j .eq. 0 .OR. j .eq. k-1) c2 = 2

val = c1/c2 * (-1.0d0)**(i+j) / (xs(i+1) - xs(j+1))

endif

adiff(i+1,j+1) = val
end do
end do

end subroutine


! subroutine chebroots(a,b,k,coefs,nroots,roots)
! implicit double precision (a-h,o-z)

! double precision                           :: coefs(1:k)
! double precision, allocatable, intent(out) :: roots(:)
! !
! !  Find the real roots of a Chebyshev expansion on the interval [a,b]
! !  via the companion matrix method.
! !
! !  Input parameters:
! !    (a,b) - the interval on which the Chebyshev expansion is given
! !    k - the number of terms in the expansion
! !    coefs - an array of length k giving the coefficients
! !
! !  Output parameters:
! !    nroots - the number of real roots in the interval [a,b]
! !    roots - an array, allocated by this routine, containing the ordered
! !      list of the roots
! !


! double precision, allocatable :: amatr(:,:)
! double precision, allocatable :: wi(:),wr(:),w(:)

! !
! !  Find the first largest nonzero coefficient
! !

! eps0 = epsilon(0.0d0)*10
! dmax = maxval(abs(coefs))

! do n=k,1,-1
! if (coefs(n)/dmax .gt. eps0) exit
! end do


! if (n <= 1) then
! nroots = 0
! allocate(roots(nroots))
! return
! endif

! !
! !  Build the companion matrix
! !

! allocate(amatr(n-1,n-1))

! amatr        = 0
! amatr(1,2)   = 1.0d0

! do j = 2,n-2
! amatr(j,j-1) = 0.5d0
! amatr(j,j+1) = 0.5d0
! end do

! do j = 1,n-1
! amatr(n-1,j) = -coefs(j)/(2*coefs(n))
! end do

! amatr(n-1,n-2) = amatr(n-1,n-2) + 0.5d0


! !
! !  Compute all of its eigenvalues via LAPACK
! !

! lw = 10000
! allocate(wr(n-1),wi(n-1),w(lw))
! call dgeev('N','N',n-1,amatr,n-1,wr,wi,vl,n-1,vr,n-1,w,lw,info)

! nroots = 0
! do i=1,n-1
! if (abs(wi(i)) .le. eps0 .AND. abs(wr(i)) .le. 1.0d0) then
! nroots = nroots+1
! endif
! end do

! allocate(roots(nroots))

! nroots = 0
! do i=1,n-1
! if (abs(wi(i)) .le. eps0 .AND. abs(wr(i)) .le. 1.0d0) then
! nroots = nroots+1
! roots(nroots) = wr(i)
! endif
! end do

! call quicksort(nroots,roots)

! !
! !  Map the roots onto the original interval
! !

! do i=1,nroots
! roots(i) = (b-a)/2 * roots(i) + (b+a)/2
! end do

! end subroutine



! subroutine chebpw_roots(ier,k,nints,ab,coefs,nroots,roots)
! implicit double precision (a-h,o-z)
! double precision                           :: ab(2,nints)
! double precision                           :: coefs(k,nints)
! double precision, allocatable, intent(out) :: roots(:)
! !
! !  Find the real roots of a function represented via a piecewise Chebyshev 
! !  expansion.
! !
! !  Input parameters:
! !    k - the number of terms in each of the expansions
! !    (nints,ab) - the arrays describing the discretization intervals
! !    coefs - the (k,nints) array giving the coefficients on each interval
! !
! !  Output parameters:
! !    nroots - the number of roots found
! !    roots - the array of roots
! !

! double precision, allocatable              :: roots0(:)
! double precision, allocatable              :: roots00(:)

! ier = 0

! maxroots = 10000
! allocate(roots0(maxroots))

! nroots0 = 0

! do int=1,nints

! a = ab(1,int)
! b = ab(2,int)
! call chebroots(a,b,k,coefs(:,int),nroots00,roots00)

! if (nroots00 + nroots0 > maxroots) then
! ier = 4
! return
! endif

! roots0(nroots0+1:nroots0+nroots00) = roots00
! nroots0 = nroots0 + nroots00

! end do



! allocate(roots(nroots0))
! roots = roots0(1:nroots0)
! nroots = nroots0

! end subroutine

end module
