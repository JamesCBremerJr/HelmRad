!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains a solver for the Riccati equation
!
!    r'(t) + (r(t))^2 + q(t) = 0,                                                          (1)
!
!  where q(t) is a user-specified coefficient.  The logarithmic derivatives of
!  solutions of 
!
!   y''(t) + q(t) y(t) = 0                                                                 (2)
!
!  satisfy (1).  The solver operates via an adaptive Chebyshev method and solutions
!  are represented via their values on piecewise Chebyshev grids.  Code for
!  constructing a phase function for (2) from the solution of (1) is also provided. 
!  A function alpha(t) is a phase function for (2) provided
!
!       sin(alpha(t))                  cos(alpha(t))
!     ----------------     and       ----------------
!      sqrt(alpha'(t))                sqrt(alpha'(t))
!
!  form a basis in the space of its solutions.
!
!  The following subroutines should be regarded as public:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    riccati_ivp - solve an initial value problem for (1); the solution is represented
!      via a discretization scheme specified by the user
!
!    riccati_tvp - solve a terminal value problem for (1); the solution is represented
!      via a discretization scheme specified by the user
!
!    riccati_ivp_adap - solve an initial value problem for (1); a discretization scheme
!      for the solution is determined adaptively
!
!    riccati_tvp_adap - solve a terminal value problem for (1); a discretization scheme
!      for the solution is determined adaptively
!
!    riccati_int_forward - integrate the solution of (1) forward to obtain the
!      logarithm of a solution of (2)
!
!    riccati_int_back - integrate the solution of (1) backward to obtain the
!      logarithm of a solution of (2)
!
!    riccati_phase_back - integrate the solution of (1) backward to construct a
!      phase function representing the solutions of (2)
!
!    riccati_phase_forward integrate the solution of (1) forward to construct
!      a phase function representing the solutions of (2)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    rriccati_ivp_adap - solve an initival value problem for (1) adaptively in the case
!      in which the solution is real-valued
!
!    rriccati_tvp_adap - solve a terminal value problem for (1) adaptively in the case
!      of a real-valued solution
!
!    rriccati_int_forward - integrate the solution of (1) forward to obtain the
!      logarithm of a solution of (1) in the case in which the logarithm is
!      real-valued
!
!    rriccati_int_back - integrate the solution of (1) backward to obtain the
!      logarithm of a solution of (1) in the case of a real-valued solution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module riccati
use chebyshev
use iso_c_binding


interface

subroutine riccati_qfun(k,ts,qs,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
integer          :: k
double precision :: ts(k)
double complex   :: qs(k)
type(c_ptr)      :: userptr
!
!  This is the interface for the user-specified subroutine which provides
!  the values of the coefficient q in (1).
!
!  Input parameters:
!    k - the number of points at which to evaluate q
!    ts - an array containing the points at which to evaluate q
!    userptr - a user provided "void *" pointer for passing information
!      to this routine
!
!  Output parameters:
!    qs - the values of the coefficient q at each of the points in ts
!
end subroutine


subroutine rriccati_qfun(k,ts,qs,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
integer          :: k
double precision :: ts(k)
double precision :: qs(k)
type(c_ptr)      :: userptr
!
!  This is the interface for the user-specified subroutine which provides
!  the values of the coefficient q in (1).
!
!  Input parameters:
!    k - the number of points at which to evaluate q
!    ts - an array containing the points at which to evaluate q
!    userptr - a user provided "void *" pointer for passing information
!      to this routine
!
!  Output parameters:
!    qs - the values of the coefficient q at each of the points in ts
!
end subroutine

end interface

contains


subroutine riccati_ivp_adap(ier,eps,chebdata,a,b,qfun,ra,nints,ab,rs,rders,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision, allocatable, intent(out)    :: ab(:,:)
double complex, allocatable, intent(out)      :: rs(:,:),rders(:,:)
double complex                                :: ra
type(c_ptr)                                   :: userptr
procedure(riccati_qfun)                       :: qfun
!
!  Adaptively solve an initial value problem for (1).
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     ra - the initial value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 
double precision, allocatable       :: ab0(:,:),ts0(:),coefs00(:)
double complex, allocatable         :: rs0(:),rders0(:),coefs0(:),qs0(:),amatr0(:,:)
double complex, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double complex                      :: ra0

double precision, allocatable       :: about(:,:)
double complex, allocatable         :: rsout(:,:),rdersout(:,:)

ier      = 0
k        = chebdata%k
maxints  = 100000
ntail    = k-4
niters   = 4

allocate(rs0(k),rders0(k),coefs0(k),ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k),coefs00(k))

nints0   = 0
allocate(ab0(2,maxints))

nintsout = 0
allocate(about(2,maxints), rsout(k,maxints), rdersout(k,maxints))


!
!  Adaptively discretize the function q to form an initial set of intervals
!

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b

nintsout = 0

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in riccati_ivp_adap,  zero length interval enountered while discretizing q, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)
coefs0 = matmul(chebdata%u,qs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / maxval(coefs00)
dd      = maxval(coefs00(ntail:k))

if (dd .lt. sqrt(eps)) then
if(nintsout+1 .gt. maxints) then
ier = 4
return
endif

nintsout = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0
else

if (nints0 + 2 .gt. maxints) then
ier = 4
return
endif

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do

!
!  Copy out the intervals, reversing their order so that we process intervals
!  on the left first
!

nints0 = nintsout

do int=1,nints0
int0 = nints0-int+1
ab0(:,int)    = about(:,int0)
end do

nintsout = 0


!
!  Now adaptively solve the initial value problem for (1)
!

do while (nints0 > 0 )
a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in riccati_ivp_adap,  zero length interval enountered, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (nintsout .eq. 0) then
ra0 = ra
else
ra0 = rsout(k,nintsout)
endif

! use the trapezoidal rule to construct an initial guess for Newton iterations
call riccati_trap_ivp(ier,k,ts0,qs0,ra0,rs0,rders0)

ifsplit = 0

do i=1,k
if (isnan(real(rs0(i)))) then
ifsplit = 1
goto 1000
endif
end do


! integrate the derivative to obtain the values of rs --- this step is crucial
rs0 = ra0 + matmul(chebdata%aintl*(b0-a0)/2,rders0)


! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs0
fs0 = -rders0 - rs0**2 - qs0
ra0  = 0
call riccati_linear_ivp(k,a0,b0,ts0,chebdata%aintl,ps0,fs0,ra0,amatr0,delta,deltap)

rs0    = rs0 + delta
rders0 = rders0 + deltap

end do


do i=1,k
if (isnan(real(rs0(i)))) then
ifsplit = 1
goto 1000
endif
end do


ifsplit = 0

coefs0 = matmul(chebdata%u,rs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / maxval( coefs00 )
dd      = maxval(coefs00(ntail:k))
if (dd .gt. eps) ifsplit = 1


! coefs0 = matmul(chebdata%u,rders0)
! coefs00 = abs(coefs0)
! coefs00 = coefs00 / maxval( coefs00 )
! dd      = maxval(coefs00(ntail:k))
! if (dd .gt. eps) ifsplit = 1


1000 continue

if (ifsplit .eq. 0) then

if (nintsout+1 .gt. maxints) then
ier = 8
return
endif

nintsout             = nintsout+1
about(1,nintsout)    = a0
about(2,nintsout)    = b0
rsout(:,nintsout)    = rs0
rdersout(:,nintsout) = rders0
else


if (nints0+2 .gt. maxints) then
ier = 8
return
endif


nints0               = nints0+1
ab0(1,nints0)        = (a0+b0)/2
ab0(2,nints0)        = b0

nints0               = nints0+1
ab0(1,nints0)        = a0
ab0(2,nints0)        = (a0+b0)/2


endif

end do

!
!  Copy out the intervals 
!


nints = nintsout
allocate(ab(2,nintsout))
allocate(rs(k,nintsout))
allocate(rders(k,nintsout))


do int=1,nints
int0             = int
ab(:,int)        = about(:,int0)
rs(:,int)        = rsout(:,int0)
rders(:,int)     = rdersout(:,int0)
end do

! call prin2("in riccati_ivp_adap, ab = ",ab)

end subroutine



subroutine riccati_ivp(ier,eps,chebdata,nints,ab,qfun,ra,rs,rders,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision, intent(in)                  :: ab(2,nints)
double complex, allocatable, intent(out)      :: rs(:,:),rders(:,:)
double complex                                :: ra
type(c_ptr)                                   :: userptr
procedure(riccati_qfun)                       :: qfun
!
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!
!     nints - the number of intervals in the discretization
!     ab - 
!
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     ra - the initial value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 
double precision, allocatable       :: ab0(:,:),ts0(:)
double complex, allocatable         :: rs0(:),rders0(:),qs0(:)
double complex, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double complex                      :: ra0
double complex, allocatable         :: amatr0(:,:)

ier      = 0
k        = chebdata%k
maxints  = 100000
niters   = 4

allocate(rs(k,nints), rders(k,nints))

allocate(ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k))


do int=1,nints
a0 = ab(1,int)
b0 = ab(2,int)

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (int .eq. 1) then
ra0 = ra
else
ra0 = rs(k,int-1)
endif

! use the trapezoidal rule to construct an initial guess for Newton iterations
call riccati_trap_ivp(ier,k,ts0,qs0,ra0,rs(:,int),rders(:,int))

! integrate the derivative to obtain the values of rs --- this step is crucial
rs(:,int) = ra0 + matmul(chebdata%aintl*(b0-a0)/2,rders(:,int))


! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs(:,int)
fs0 = -rders(:,int) - rs(:,int)**2 - qs0
ra0  = 0
call riccati_linear_ivp(k,a0,b0,ts0,chebdata%aintl,ps0,fs0,ra0,amatr0,delta,deltap)

rs(:,int)    = rs(:,int) + delta
rders(:,int) = rders(:,int) + deltap

end do

end do

end subroutine



subroutine riccati_trap_ivp(ier,k,ts,qs,ra,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k)
double complex           :: rs(k),rders(k), qs(k), ra
!
!  Use the implicit trapezoidal method to solve the initial value problem
! 
!     r'(t) + (r(t))^2 + q(t) = 0      a < t < b
!     r(a)                    = ra
!
!  The solution is tabulated at a collection of nodes in the interval (a,b)
!  specified by the user the left-most of which must be the left endpoint
!  a.
!
!  Input parameters:
!    k - the number of nodes at which to approximate the solution
!    ts - an array containing the nodes
!    qs - the values of the function q at the nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    ier - an error return code;
!      ier =    0    indicates successful execution
!      ier = 1024    means that NaN was encountered and the procedure aborted
!
!    rs - the values of the solution at the specified nodes
!    rders - the values of the derivative of the solution at the specified
!      nodes
!    
double complex   :: f0,r,rhs,delta,q0,q1

ier    = 0
niters = 4

rs(1)    = ra
rders(1) = -qs(1) - ra**2

do j=2,k
t0     = ts(j-1)
t1     = ts(j)

q0     = qs(j-1)
q1     = qs(j)

r      = rs(j-1)
h      = t1-t0
f0     = -q0-r**2
rhs    = r + h/2 * f0 - h/2*q1

do iter=1,niters
delta  = (r+h/2*r**2-rhs) / (1 + h*r)
r      = r - delta

! if (r+1 .eq. r) then
! ier = 1024
! return
! endif

end do

rs(j)    = r
rders(j) = -q1-r**2
end do

end subroutine


subroutine riccati_linear_ivp(k,a,b,ts,chebint,qs,fs,ra,amatr,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k),chebint(k,k)
double complex           :: amatr(k,k),qs(k),fs(k),rs(k),rders(k),ra
!
!  Use a Chebyshev spectral method to solve the linear IVP
!
!    r'(t) + qs(t) r(t) = f(t)
!    r(a)               = ra
!
!  over the interval (a,b).
!
!  Input parameters:
!    k - the number of Chebyshev nodes
!    ts - the k Chebyshev nodes on the interval (a,b)
!    chebint - the "left" Chebyshev spectral integration matrix
!    qs - the values of the coefficient at the Chebyshev nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    rs - the values of the solution at the  Chebyshevnodes
!    rders - the values of the derivative of the solution at the 
!      Chebyshev nodes
!
!  Work arrays:
!    amatr - a (k,k) matrix 
!

!a = ts(1)
!b = ts(k)


!
!  We solve the integral equation 
!
!     \sigma(t) + q(t) \int_a^t \sigma(s) ds = f(t) - ra q(t)
!
!  obtained by letting
!
!     r(t) = ra + \int_a^t \sigma(s) ds.
! 


rders = 0
amatr = 0

do i=1,k
amatr(i,i) = 1.0d0
end do

do i=1,k
rders(i)   = fs(i) - ra * qs(i)
amatr(i,:) = amatr(i,:) + qs(i) * chebint(i,:)*(b-a)/2
end do


call riccati_cqrsolv(k,amatr,rders)
rs    = ra + (b-a)/2*matmul(chebint,rders)

end subroutine



subroutine riccati_tvp_adap(ier,eps,chebdata,a,b,qfun,rb,nints,ab,rs,rders,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints 
double precision, allocatable, intent(out)    :: ab(:,:)
double complex, allocatable, intent(out)      :: rs(:,:),rders(:,:)
double complex                                :: rb
type(c_ptr)                                   :: userptr
procedure(riccati_qfun)                       :: qfun
!
!  Adaptively solve a terminal value problem for (1).
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 
double precision, allocatable       :: ab0(:,:),ts0(:),coefs00(:)
double complex, allocatable         :: rs0(:),rders0(:),coefs0(:),qs0(:),amatr0(:,:)
double complex, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double complex                      :: rb0

double precision, allocatable       :: about(:,:)
double complex, allocatable         :: rsout(:,:),rdersout(:,:)

ier       = 0
k         = chebdata%k
ntail     = k-4
niters    = 4
maxints   = 10000


allocate(rs0(k),rders0(k),coefs0(k),ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k),coefs00(k))

nints0   = 0
allocate(ab0(2,maxints))

nintsout = 0
allocate(about(2,maxints), rsout(k,maxints), rdersout(k,maxints))


!
!  Adaptively discretize the function q to form an initial set of intervals
!

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b

nintsout = 0

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in riccati_tvp_adap,  zero length interval enountered while discretizing q, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

coefs0  = matmul(chebdata%u,qs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / (maxval(coefs00)+1)
dd      = maxval(coefs00(ntail:k))


if (dd .lt. sqrt(eps)) then
if(nintsout+1 .gt. maxints) then
ier = 4
return
endif


nintsout = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0
else


if (nints0 + 2 .gt. maxints) then
ier = 4
return
endif

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do

!
!  Copy out the intervals
!

nints0 = nintsout

do int=1,nints0
int0 = int
ab0(:,int)    = about(:,int0)
end do

nintsout = 0



!
!  Now adaptively solve the initial value problem for (1)
!

do while (nints0 > 0 )
a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in riccati_tvp_adap,  zero length interval enountered, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (nintsout .eq. 0) then
rb0 = rb
else
rb0 = rsout(1,nintsout)
endif




! use the trapezoidal rule to construct an initial guess for Newton iterations
call riccati_trap_tvp(ier,k,ts0,qs0,rb0,rs0,rders0)


ifsplit = 0

do i=1,k
if (isnan(real(rs0(i)))) then
ifsplit = 1
goto 1000
endif
end do


if (ier .eq. 1024) then
call prin2("in riccati_tvp_adap, NaN encounted, a0 = ",a0)
call prin2("in riccati_tvp_adap, NaN encounted, b0 = ",b0)
return
endif



! integrate the derivative to obtain the values of rs --- this step is crucial
rs0 = rb0 + (b0-a0)/2*matmul(chebdata%aintr,rders0)

! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs0
fs0 = -rders0 - rs0**2 - qs0
rb0  = 0
call riccati_linear_tvp(k,a0,b0,ts0,chebdata%aintr,ps0,fs0,rb0,amatr0,delta,deltap)

rs0    = rs0 + delta
rders0 = rders0 + deltap

end do


do i=1,k
if (isnan(real(rs0(i)))) then
ifsplit = 1
goto 1000
endif
end do


coefs0  = matmul(chebdata%u,rs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / (maxval(coefs00)+1)
dd      = maxval(coefs00(ntail:k))
if (dd .gt. eps) ifsplit = 1



1000 continue


if(ifsplit .eq. 0) then

if (nintsout+1 .gt. maxints) then
ier = 8
return
endif

nintsout             = nintsout+1
about(1,nintsout)    = a0
about(2,nintsout)    = b0
rsout(:,nintsout)    = rs0
rdersout(:,nintsout) = rders0

else


if (nints0+2 .gt. maxints) then
ier = 8
return
endif


nints0               = nints0+1
ab0(1,nints0)        = a0
ab0(2,nints0)        = (a0+b0)/2

nints0               = nints0+1
ab0(1,nints0)        = (a0+b0)/2
ab0(2,nints0)        = b0



endif

end do

!
!  Copy out the intervals in reverse order
!


nints = nintsout
allocate(ab(2,nintsout))
allocate(rs(k,nintsout))
allocate(rders(k,nintsout))


do int=1,nints
int0             = nints-int+1
ab(:,int)        = about(:,int0)
rs(:,int)        = rsout(:,int0)
rders(:,int)     = rdersout(:,int0)
end do

!call prin2("in riccati_tvp_adap, ab = ",ab)

end subroutine



subroutine riccati_tvp(ier,eps,chebdata,nints,ab,qfun,rb,rs,rders,userptr)
implicit double precision (a-h,o-z)

type(chebexps_data)                           :: chebdata
integer                                       :: nints 
double precision, intent(in)                  :: ab(2,nints)
double complex, allocatable, intent(out)      :: rs(:,:),rders(:,:)
double complex                                :: rb
type(c_ptr)                                   :: userptr
procedure(riccati_qfun)                       :: qfun

!
!  
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 

double precision, allocatable       :: ab0(:,:),ts0(:)
double complex, allocatable         :: qs0(:),amatr0(:,:)
double complex, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double complex                      :: rb0

ier       = 0
k         = chebdata%k
niters    = 4
maxints   = 10000

allocate(rs(k,nints),rders(k,nints))
allocate(ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k))

!
!  Now adaptively solve the initial value problem for (1)
!

do int=nints,1,-1

a0 = ab(1,int)
b0 = ab(2,int)

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (int .eq. nints) then
rb0 = rb
else
rb0 = rs(1,int+1)
endif


! use the trapezoidal rule to construct an initial guess for Newton iterations
call riccati_trap_tvp(ier,k,ts0,qs0,rb0,rs(:,int),rders(:,int))



! integrate the derivative to obtain the values of rs --- this step is crucial
rs(:,int) = rb0 + (b0-a0)/2*matmul(chebdata%aintr,rders(:,int))


! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs(:,int)
fs0 = -rders(:,int) - rs(:,int)**2 - qs0
rb0  = 0
call riccati_linear_tvp(k,a0,b0,ts0,chebdata%aintr,ps0,fs0,rb0,amatr0,delta,deltap)

rs(:,int)    = rs(:,int) + delta
rders(:,int) = rders(:,int) + deltap
end do
end do

end subroutine


subroutine riccati_trap_tvp(ier,k,ts,qs,rb,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k)
double complex           :: rs(k),rders(k), qs(k), rb
!
!  Use the implicit trapezoidal method to solve the terminal value problem
! 
!     r'(t) + (r(t))^2 + q(t) = 0      a < t < b
!     r(b)                    = rb
!
!  The solution is tabulated at a collection of nodes in the interval (a,b)
!  specified by the user the right-most of which must be the right  endpoint
!  b.
!
!  Input parameters:
!    k - the number of nodes at which to approximate the solution
!    ts - an array containing the nodes
!    qs - the values of the function q at the nodes
!    rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code;
!      ier =    0    indicates successful execution
!      ier = 1024    means that NaN was encountered and the procedure aborted
!
!    rs - the values of the solution at the specified nodes
!    rders - the values of the derivative of the solution at the specified
!      nodes
!    
double complex   :: f0,r,rhs,delta,q0,q1

ier    = 0
niters = 4
rs(k)    = rb
rders(k) = -qs(k) - rb**2

do j=k-1,1,-1

t0 = ts(j)
t1 = ts(j+1)
h  = t1-t0

q0 = qs(j)
q1 = qs(j+1)

r   = rs(j+1)
rhs = rs(j+1) + h/2 * rs(j+1)**2 + h/2 * q1 + h/2*q0

do iter=1,niters
delta = (rhs - r +h/2*r**2 ) / (1 - h*r)
r     = r + delta
end do

rs(j)    = r
rders(j) = -r**2 - q0
end do


end subroutine


subroutine riccati_linear_tvp(k,a,b,ts,chebint,qs,fs,rb,amatr,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k),chebint(k,k)
double complex           :: amatr(k,k),qs(k),fs(k),rs(k),rders(k),rb
!
!  Use a Chebyshev spectral method to solve the terminal value problem
!
!    r'(t) + qs(t) r(t) = f(t)
!    r(b)               = rb
!
!  over the interval (a,b).
!
!  Input parameters:
!    k - the number of Chebyshev nodes
!    ts - the k Chebyshev nodes on the interval (a,b)
!    chebint - the "right" Chebyshev spectral integration matrix
!    qs - the values of the coefficient at the Chebyshev nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    rs - the values of the solution at the  Chebyshevnodes
!    rders - the values of the derivative of the solution at the 
!      Chebyshev nodes
!
!  Work arrays:
!    amatr - a (k,k) matrix 
!



!
!  We solve the integral equation 
!
!     \sigma(t) + q(t) \int_b^t \sigma(s) ds = f(t) - rb q(t)
!
!  obtained by letting
!
!     r(t) = rb + \int_b^t \sigma(s) ds.
! 


rders = 0
amatr = 0

do i=1,k
amatr(i,i) = 1.0d0
end do

do i=1,k
rders(i)   = fs(i) - rb * qs(i)
amatr(i,:) = amatr(i,:) + qs(i) * chebint(i,:)*(b-a)/2
end do


call riccati_cqrsolv(k,amatr,rders)
rs    = rb + (b-a)/2*matmul(chebint,rders)

end subroutine



subroutine riccati_cqrsolv(n,a,rhs)
implicit complex *16 (a-h,o-z)
double complex         :: a(n,n),u(2,2),rhs(n)
double precision       :: rcond
! 
!  This subroutine uses a version of QR-decomposition to solve the user-supplied
!  system of linear algebraic equations Ax = y.   The matrix  a is destroyed 
!  in the process and the parameter supplying the values of the right-hand
!  side y is overwritten with the values of the solution x.
! 
!  Input parameters:
!    n - the dimensionality of the system being solved
!    a - the (n,n) complex-valued coefficient matrix WHICH WILL
!      BE DESTROYED BY THIS ROUUTINE
!    rhs - the right-hand side of the system to be solved, which
!      is overwritten by this subroutine
! 
!  Output parameters:
!    rhs - the solution of the linear system
!

double precision       ::  dmax,dmin,size22,dd,eps0
double complex         ::  aa(2)
integer, allocatable   ::  ipiv(:)

eps0 = epsilon(0.0d0)

! if (eps0 .gt. 1.0d-17) then
! allocate(ipiv(n))
! call zgesv(n,1,a,n,ipiv,rhs,n,info)
! if (info .ne. 0) then
! call prini("after zgesv, info = ",info)
! stop
! endif
! return
! endif

!
! Transpose the input matrix a
! 

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)*conjg(a(j,i))
size22=size22+a(i,j)*conjg(a(i,j))
end do
end do

! 
! Eliminate the upper right triangle
! 
do i=1,n-1
! 
do j=n,i+1,-1
! 
aa(1)=a(i,j-1)
aa(2)=a(i,j)
u21=-aa(2)
u22=aa(1)
dd=u22*conjg(u22)+u21*conjg(u21)
if(dd .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
dd=sqrt(dd)
u(2,2)=u22/dd
u(2,1)=u21/dd
u(1,1)=-conjg(u(2,2))
u(1,2)=conjg(u(2,1))
endif


do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j)
a(ii,j-1) = d1
a(ii,j)   = d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do
! 
! Estimate the condition number
! 
dmax=-1
dmin=abs(a(1,1))
do i=1,n
if(dmax .lt. abs(a(i,i)) ) dmax=abs(a(i,i))
if(dmin .gt. abs(a(i,i)) ) dmin=abs(a(i,i))
end do
if(dmin .eq. 0) then
rcond=-1
return
endif
rcond=dmax/dmin

!
!  Apply the inverse of the triangular matrix to rhs
!

!call cqrtrin(a,rhs,n)

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

end subroutine



subroutine riccati_phase_back(chebdata,nints,ab,rs,rders,aval,avals,apvals,appvals)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision, allocatable, intent(in)     :: ab(:,:)
double complex, allocatable, intent(in)       :: rs(:,:),rders(:,:)
double precision, allocatable, intent(out)    :: avals(:,:),apvals(:,:),appvals(:,:)
double precision                              :: aval

!
!  Integrate the solution of (1) backward to construct a phase function for (2).
!
!  Input parameters:
!    chebdata - the structure returned by chebexps
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!    aval - the value of the phase function at the point b
!
!  Output parameters:
!    avals - a (k,nints) array giving the values of the phase function at the
!     Chebyshev nodes
!    apvals - a (k,nints) array giving the values of the derivative of the
!     of the phase function at the Chebyshev nodes


k = chebdata%k
allocate(avals(k,nints),apvals(k,nints),appvals(k,nints))

aval0 = aval

do int = nints,1,-1
a = ab(1,int)
b = ab(2,int)

apvals(:,int)  = imag(rs(:,int))
appvals(:,int) = imag(rders(:,int))

avals(:,int)   = aval0 + (b-a)/2 * matmul(chebdata%aintr,apvals(:,int))
aval0          = avals(1,int)


end do


end subroutine



subroutine riccati_phase_forward(chebdata,nints,ab,rs,rders,aval,avals,apvals,appvals)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision, allocatable, intent(in)     :: ab(:,:)
double complex, allocatable, intent(in)       :: rs(:,:),rders(:,:)
double precision, allocatable, intent(out)    :: avals(:,:),apvals(:,:),appvals(:,:)
double precision                              :: aval

!
!  Integrate the solution of (1) forward to construct a phase function for (2).
!
!  Input parameters:
!    chebdata - the structure returned by chebexps
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!    aval - the value of the phase function at the point b
!
!  Output parameters:
!    avals - a (k,nints) array giving the values of the phase function at the
!     Chebyshev nodes
!    apvals - a (k,nints) array giving the values of the derivative of the
!     of the phase function at the Chebyshev nodes


k = chebdata%k
allocate(avals(k,nints),apvals(k,nints),appvals(k,nints))

aval0 = aval

do int = 1,nints
a = ab(1,int)
b = ab(2,int)

apvals(:,int)  = imag(rs(:,int))
appvals(:,int) = imag(rders(:,int))

avals(:,int)   = aval0 + (b-a)/2 * matmul(chebdata%aintl,apvals(:,int))
aval0          = avals(k,int)


end do

end subroutine



subroutine riccati_int_forward(chebdata,nints,ab,rs,rders,psival,psis)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision, allocatable, intent(in)     :: ab(:,:)
double complex, allocatable, intent(in)       :: rs(:,:),rders(:,:)
double complex, allocatable, intent(out)      :: psis(:,:)
double complex                                :: psival
!
!  Integrate the solution of (1) forward; that is, construct the function
!
!                          t
!    psi(t) = psival + \int r(u) du 
!                          a
!
!  Input parameters:
!    chebdata - the structure returned by chebexps
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!    psival - the value of psi at a
!
!  Output parameters:
!    psis  the values of the function psi at the Chebyshev nodes


double complex                  :: psival0


k = chebdata%k
allocate(psis(k,nints))

psival0 = psival

do int = 1,nints

a = ab(1,int)
b = ab(2,int)

psis(:,int) = psival0 + (b-a)/2 * matmul(chebdata%aintl,rs(:,int))
psival0     = psis(k,int)
end do

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine rriccati_ivp_adap(ier,eps,chebdata,a,b,qfun,ra,nints,ab,rs,rders,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision, allocatable, intent(out)    :: ab(:,:)
double precision, allocatable, intent(out)    :: rs(:,:),rders(:,:)
double precision                              :: ra
type(c_ptr)                                   :: userptr
procedure(rriccati_qfun)                      :: qfun
!
!  Adaptively solve an initial value problem for (1).
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     ra - the initial value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 
double precision, allocatable         :: ab0(:,:),ts0(:),coefs00(:)
double precision, allocatable         :: rs0(:),rders0(:),coefs0(:),qs0(:),amatr0(:,:)
double precision, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double precision                      :: ra0

double precision, allocatable         :: about(:,:)
double precision, allocatable         :: rsout(:,:),rdersout(:,:)

ier      = 0
k        = chebdata%k
maxints  = 100000
ntail    = k-4
niters   = 4
eps0     = epsilon(0.0d0)

allocate(rs0(k),rders0(k),coefs0(k),ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k),coefs00(k))

nints0   = 0
allocate(ab0(2,maxints))

nintsout = 0
allocate(about(2,maxints), rsout(k,maxints), rdersout(k,maxints))


!
!  Adaptively discretize the function q to form an initial set of intervals
!

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b

nintsout = 0

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1

if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in rriccati_ivp_adap,  zero length interval enountered while discretizing q, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)
coefs0 = matmul(chebdata%u,qs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / (maxval(coefs00)+1)
dd      = maxval(coefs00(ntail:k))


if (dd .lt. sqrt(eps)) then

if(nintsout+1 .gt. maxints) then
ier = 4
return
endif

nintsout = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0
else

if (nints0 + 2 .gt. maxints) then
ier = 4
return
endif

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do

!
!  Copy out the intervals, reversing their order so that we process intervals
!  on the left first
!

nints0 = nintsout

do int=1,nints0
int0 = nints0-int+1
ab0(:,int)    = about(:,int0)
end do

nintsout = 0


!
!  Now adaptively solve the initial value problem for (1)
!

do while (nints0 > 0 )
dd     = 1
a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1

if ( (b0 - a0) .lt. eps0*a0) then
ier = 16
call prin2("in rriccati_ivp_adap,  zero length interval enountered, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (nintsout .eq. 0) then
ra0 = ra
else
ra0 = rsout(k,nintsout)
endif

! use the trapezoidal rule to construct an initial guess for Newton iterations
call rriccati_trap_ivp(ier,k,ts0,qs0,ra0,rs0,rders0)


ifsplit = 0

do i=1,k
if (isnan(rs0(i))) then
ifsplit = 1
goto 1000
endif
end do


! integrate the derivative to obtain the values of rs --- this step is crucial
rs0 = ra0 + matmul(chebdata%aintl*(b0-a0)/2,rders0)


! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs0
fs0 = -rders0 - rs0**2 - qs0
ra0  = 0
call rriccati_linear_ivp(k,a0,b0,ts0,chebdata%aintl,ps0,fs0,ra0,amatr0,delta,deltap)

rs0    = rs0    + delta
rders0 = rders0 + deltap

end do


do i=1,k
if (isnan(rs0(i))) then
ifsplit = 1
goto 1000
endif
end do


ifsplit = 0

coefs0  = matmul(chebdata%u,rs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / (maxval( coefs00 )+1)
dd      = maxval(coefs00(ntail:k))
if (dd .gt. eps) ifsplit = 1


1000 continue

if (ifsplit .eq. 0) then

if (nintsout+1 .gt. maxints) then
ier = 8
return
endif

nintsout             = nintsout+1
about(1,nintsout)    = a0
about(2,nintsout)    = b0
rsout(:,nintsout)    = rs0
rdersout(:,nintsout) = rders0

else

if (nints0+2 .gt. maxints) then
ier = 8
return
endif


nints0               = nints0+1
ab0(1,nints0)        = (a0+b0)/2
ab0(2,nints0)        = b0

nints0               = nints0+1
ab0(1,nints0)        = a0
ab0(2,nints0)        = (a0+b0)/2


endif

end do

!
!  Copy out the intervals 
!


nints = nintsout
allocate(ab(2,nintsout))
allocate(rs(k,nintsout))
allocate(rders(k,nintsout))


do int=1,nints
int0             = int
ab(:,int)        = about(:,int0)
rs(:,int)        = rsout(:,int0)
rders(:,int)     = rdersout(:,int0)
end do

! call prin2("in riccati_ivp_adap, ab = ",ab)

end subroutine


subroutine rriccati_ivp(ier,eps,chebdata,nints,ab,qfun,ra,rs,rders,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints
double precision                              :: ab(2,nints)
double precision, allocatable, intent(out)    :: rs(:,:),rders(:,:)
double precision                              :: ra
type(c_ptr)                                   :: userptr
procedure(rriccati_qfun)                      :: qfun
!
!  Adaptively solve an initial value problem for (1).
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     ra - the initial value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 
double precision, allocatable         :: ts0(:)
double precision, allocatable         :: qs0(:),amatr0(:,:)
double precision, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double precision                      :: ra0

double precision, allocatable         :: about(:,:)
double precision, allocatable         :: rsout(:,:),rdersout(:,:)

ier      = 0
k        = chebdata%k
maxints  = 100000
niters   = 4

allocate(rs(k,nints),rders(k,nints))
allocate(ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k))


do int=1,nints
a0 = ab(1,int)
b0 = ab(2,int)

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (int .eq. 1) then
ra0 = ra
else
ra0 = rs(k,int-1)
endif


! use the trapezoidal rule to construct an initial guess for Newton iterations
call rriccati_trap_ivp(ier,k,ts0,qs0,ra0,rs(:,int),rders(:,int))

if (ier .ne. 0) then
call prini("after rriccati_trap_ivp, ier = ",ier)
stop
endif

! integrate the derivative to obtain the values of rs --- this step is crucial
rs(:,int) = ra0 + matmul(chebdata%aintl*(b0-a0)/2,rders(:,int))


! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs(:,int)
fs0 = -rders(:,int) - rs(:,int)**2 - qs0
ra0  = 0
call rriccati_linear_ivp(k,a0,b0,ts0,chebdata%aintl,ps0,fs0,ra0,amatr0,delta,deltap)

rs(:,int)    = rs(:,int) + delta
rders(:,int) = rders(:,int) + deltap

end do


end do

end subroutine




subroutine rriccati_trap_ivp(ier,k,ts,qs,ra,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k)
double precision          :: rs(k),rders(k), qs(k), ra
!
!  Use the implicit trapezoidal method to solve the initial value problem
! 
!     r'(t) + (r(t))^2 + q(t) = 0      a < t < b
!     r(a)                    = ra
!
!  The solution is tabulated at a collection of nodes in the interval (a,b)
!  specified by the user the left-most of which must be the left endpoint
!  a.
!
!  Input parameters:
!    k - the number of nodes at which to approximate the solution
!    ts - an array containing the nodes
!    qs - the values of the function q at the nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    ier - an error return code;
!      ier =    0    indicates successful execution
!      ier = 1024    means that NaN was encountered and the procedure aborted
!
!    rs - the values of the solution at the specified nodes
!    rders - the values of the derivative of the solution at the specified
!      nodes
!    
!double complex   :: f0,r,rhs,delta,q0,q1

ier    = 0
niters = 4

rs(1)    = ra
rders(1) = -qs(1) - ra**2

do j=2,k
t0     = ts(j-1)
t1     = ts(j)

q0     = qs(j-1)
q1     = qs(j)

r      = rs(j-1)
h      = t1-t0
f0     = -q0-r**2
rhs    = r + h/2 * f0 - h/2*q1

do iter=1,niters
delta  = (r+h/2*r**2-rhs) / (1 + h*r)
r      = r - delta

end do

rs(j)    = r
rders(j) = -q1-r**2
end do

end subroutine


subroutine rriccati_linear_ivp(k,a,b,ts,chebint,qs,fs,ra,amatr,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k),chebint(k,k)
double precision         :: amatr(k,k),qs(k),fs(k),rs(k),rders(k),ra
!
!  Use a Chebyshev spectral method to solve the linear IVP
!
!    r'(t) + qs(t) r(t) = f(t)
!    r(a)               = ra
!
!  over the interval (a,b).
!
!  Input parameters:
!    k - the number of Chebyshev nodes
!    ts - the k Chebyshev nodes on the interval (a,b)
!    chebint - the "left" Chebyshev spectral integration matrix
!    qs - the values of the coefficient at the Chebyshev nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    rs - the values of the solution at the  Chebyshevnodes
!    rders - the values of the derivative of the solution at the 
!      Chebyshev nodes
!
!  Work arrays:
!    amatr - a (k,k) matrix 
!

!a = ts(1)
!b = ts(k)


!
!  We solve the integral equation 
!
!     \sigma(t) + q(t) \int_a^t \sigma(s) ds = f(t) - ra q(t)
!
!  obtained by letting
!
!     r(t) = ra + \int_a^t \sigma(s) ds.
! 


rders = 0
amatr = 0

do i=1,k
amatr(i,i) = 1.0d0
end do

do i=1,k
rders(i)   = fs(i) - ra * qs(i)
amatr(i,:) = amatr(i,:) + qs(i) * chebint(i,:)*(b-a)/2
end do


! call riccati_cqrsolv(k,amatr,rders)
call qrsolv(amatr,k,rders)
rs    = ra + (b-a)/2*matmul(chebint,rders)

end subroutine



subroutine rriccati_int_forward(chebdata,nints,ab,rs,rders,psival,psis)
implicit double precision (a-h,o-z)
type(chebexps_data)                             :: chebdata
integer                                         :: nints
double precision, allocatable, intent(in)       :: ab(:,:)
double precision, allocatable, intent(in)       :: rs(:,:),rders(:,:)
double precision, allocatable, intent(out)      :: psis(:,:)
double precision                                :: psival
!
!  Integrate the solution of (1) forward; that is, construct the function
!
!                          t
!    psi(t) = psival + \int r(u) du 
!                          a
!
!  Input parameters:
!    chebdata - the structure returned by chebexps
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!    psival - the value of psi at a
!
!  Output parameters:
!    psis  the values of the function psi at the Chebyshev nodes


!double complex                  :: psival0


k = chebdata%k
allocate(psis(k,nints))

psival0 = psival

do int = 1,nints

a = ab(1,int)
b = ab(2,int)

psis(:,int) = psival0 + (b-a)/2 * matmul(chebdata%aintl,rs(:,int))
psival0     = psis(k,int)
end do

end subroutine




subroutine rriccati_tvp_adap(ier,eps,chebdata,a,b,qfun,rb,nints,ab,rs,rders,userptr)
implicit double precision (a-h,o-z)
type(chebexps_data)                           :: chebdata
integer                                       :: nints 
double precision, allocatable, intent(out)    :: ab(:,:)
double precision, allocatable, intent(out)    :: rs(:,:),rders(:,:)
double precision                              :: rb
type(c_ptr)                                   :: userptr
procedure(rriccati_qfun)                      :: qfun
!
!  Adaptively solve a terminal value problem for (1).
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 
double precision, allocatable       :: ab0(:,:),ts0(:),coefs00(:)
double precision, allocatable         :: rs0(:),rders0(:),coefs0(:),qs0(:),amatr0(:,:)
double precision, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double precision                      :: rb0

double precision, allocatable       :: about(:,:)
double precision, allocatable         :: rsout(:,:),rdersout(:,:)

ier       = 0
k         = chebdata%k
ntail     = k-4
niters    = 4
maxints   = 10000


allocate(rs0(k),rders0(k),coefs0(k),ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k),coefs00(k))

nints0   = 0
allocate(ab0(2,maxints))

nintsout = 0
allocate(about(2,maxints), rsout(k,maxints), rdersout(k,maxints))


!
!  Adaptively discretize the function q to form an initial set of intervals
!

nints0   = 1
ab0(1,1) = a
ab0(2,1) = b

nintsout = 0

do while (nints0 > 0)

a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in rriccati_tvp_adap,  zero length interval enountered while discretizing q, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

coefs0  = matmul(chebdata%u,qs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / (maxval(coefs00)+1)
dd      = maxval(coefs00(ntail:k))


if (dd .lt. sqrt(eps)) then
if(nintsout+1 .gt. maxints) then
ier = 4
return
endif


nintsout = nintsout+1
about(1,nintsout) = a0
about(2,nintsout) = b0
else


if (nints0 + 2 .gt. maxints) then
ier = 4
return
endif

nints0 = nints0+1
ab0(1,nints0) = (a0+b0)/2
ab0(2,nints0) = b0

nints0 = nints0+1
ab0(1,nints0) = a0
ab0(2,nints0) = (a0+b0)/2

endif

end do

!
!  Copy out the intervals
!

nints0 = nintsout

do int=1,nints0
int0 = int
ab0(:,int)    = about(:,int0)
end do

nintsout = 0



!
!  Now adaptively solve the initial value problem for (1)
!

do while (nints0 > 0 )
a0 = ab0(1,nints0)
b0 = ab0(2,nints0)
nints0 = nints0-1


if ( (b0 - a0) .eq. 0) then
ier = 16
call prin2("in rriccati_tvp_adap,  zero length interval enountered, a0 = ",a0)
return
endif

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (nintsout .eq. 0) then
rb0 = rb
else
rb0 = rsout(1,nintsout)
endif




! use the trapezoidal rule to construct an initial guess for Newton iterations
call rriccati_trap_tvp(ier,k,ts0,qs0,rb0,rs0,rders0)


ifsplit = 0

do i=1,k
if (isnan(real(rs0(i)))) then
ifsplit = 1
goto 1000
endif
end do


if (ier .eq. 1024) then
call prin2("in riccati_tvp_adap, NaN encounted, a0 = ",a0)
call prin2("in riccati_tvp_adap, NaN encounted, b0 = ",b0)
return
endif



! integrate the derivative to obtain the values of rs --- this step is crucial
rs0 = rb0 + (b0-a0)/2*matmul(chebdata%aintr,rders0)

! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs0
fs0 = -rders0 - rs0**2 - qs0
rb0  = 0
call rriccati_linear_tvp(k,a0,b0,ts0,chebdata%aintr,ps0,fs0,rb0,amatr0,delta,deltap)

rs0    = rs0 + delta
rders0 = rders0 + deltap

end do


do i=1,k
if (isnan(rs0(i))) then
ifsplit = 1
goto 1000
endif
end do


coefs0  = matmul(chebdata%u,rs0)
coefs00 = abs(coefs0)
coefs00 = coefs00 / (maxval(coefs00)+1)
dd      = maxval(coefs00(ntail:k))
if (dd .gt. eps) ifsplit = 1

1000 continue

if(ifsplit .eq. 0) then

if (nintsout+1 .gt. maxints) then
ier = 8
return
endif

nintsout             = nintsout+1
about(1,nintsout)    = a0
about(2,nintsout)    = b0
rsout(:,nintsout)    = rs0
rdersout(:,nintsout) = rders0

else


if (nints0+2 .gt. maxints) then
ier = 8
return
endif


nints0               = nints0+1
ab0(1,nints0)        = a0
ab0(2,nints0)        = (a0+b0)/2

nints0               = nints0+1
ab0(1,nints0)        = (a0+b0)/2
ab0(2,nints0)        = b0



endif

end do

!
!  Copy out the intervals in reverse order
!


nints = nintsout
allocate(ab(2,nintsout))
allocate(rs(k,nintsout))
allocate(rders(k,nintsout))


do int=1,nints
int0             = nints-int+1
ab(:,int)        = about(:,int0)
rs(:,int)        = rsout(:,int0)
rders(:,int)     = rdersout(:,int0)
end do

!call prin2("in riccati_tvp_adap, ab = ",ab)

end subroutine



subroutine rriccati_tvp(ier,eps,chebdata,nints,ab,qfun,rb,rs,rders,userptr)
implicit double precision (a-h,o-z)

type(chebexps_data)                           :: chebdata
integer                                       :: nints 
double precision, intent(in)                  :: ab(2,nints)
double precision, allocatable, intent(out)    :: rs(:,:),rders(:,:)
double precision                              :: rb
type(c_ptr)                                   :: userptr
procedure(rriccati_qfun)                      :: qfun

!
!  
!
!  Input parameters:
!     eps - precision for the adaptive discretization procedure
!     chebdata - the structure returned by chebexps
!     (a,b) - the interval on which (1) is given
!     qfun - user-specified function conforming to the riccati_qfun
!       interface which provides the values of q
!     rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code
!      ier = 0       indicates successful execution
!      ier = 4       means that the adaptive discretization of q failed
!      ier = 8       means that the maximum number of intervals was exceeded
!      ier = 16      means that an interval of length 0 was encountered; its location
!                    will be displayed
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!
! 

double precision, allocatable        :: ab0(:,:),ts0(:)
double precision, allocatable         :: qs0(:),amatr0(:,:)
double precision, allocatable         :: ps0(:),deltap(:),delta(:),fs0(:)
double precision                      :: rb0

ier       = 0
k         = chebdata%k
niters    = 4
maxints   = 10000

allocate(rs(k,nints),rders(k,nints))
allocate(ts0(k),qs0(k),delta(k),deltap(k))
allocate(ps0(k),fs0(k),amatr0(k,k))

!
!  Now adaptively solve the initial value problem for (1)
!

do int=nints,1,-1

a0 = ab(1,int)
b0 = ab(2,int)

ts0 = (b0-a0)/2 * chebdata%xs + (b0+a0)/2
call qfun(k,ts0,qs0,userptr)

if (int .eq. nints) then
rb0 = rb
else
rb0 = rs(1,int+1)
endif


! use the trapezoidal rule to construct an initial guess for Newton iterations
call rriccati_trap_tvp(ier,k,ts0,qs0,rb0,rs(:,int),rders(:,int))


! integrate the derivative to obtain the values of rs --- this step is crucial
rs(:,int) = rb0 + (b0-a0)/2*matmul(chebdata%aintr,rders(:,int))


! perform Newton iterations 
do iter=1,niters

ps0 = 2 *rs(:,int)
fs0 = -rders(:,int) - rs(:,int)**2 - qs0
rb0  = 0
call rriccati_linear_tvp(k,a0,b0,ts0,chebdata%aintr,ps0,fs0,rb0,amatr0,delta,deltap)

rs(:,int)    = rs(:,int) + delta
rders(:,int) = rders(:,int) + deltap
end do
end do

end subroutine


subroutine rriccati_trap_tvp(ier,k,ts,qs,rb,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k)
double precision         :: rs(k),rders(k), qs(k), rb
!
!  Use the implicit trapezoidal method to solve the terminal value problem
! 
!     r'(t) + (r(t))^2 + q(t) = 0      a < t < b
!     r(b)                    = rb
!
!  The solution is tabulated at a collection of nodes in the interval (a,b)
!  specified by the user the right-most of which must be the right  endpoint
!  b.
!
!  Input parameters:
!    k - the number of nodes at which to approximate the solution
!    ts - an array containing the nodes
!    qs - the values of the function q at the nodes
!    rb - the terminal value for the solution
!
!  Output parameters:
!    ier - an error return code;
!      ier =    0    indicates successful execution
!      ier = 1024    means that NaN was encountered and the procedure aborted
!
!    rs - the values of the solution at the specified nodes
!    rders - the values of the derivative of the solution at the specified
!      nodes
!    
double complex   :: f0,r,rhs,delta,q0,q1

ier    = 0
niters = 4

! the number of Newton iterations for each step

rs(k)    = rb
rders(k) = -qs(k) - rb**2

do j=k-1,1,-1

t0 = ts(j)
t1 = ts(j+1)
h  = t1-t0

q0 = qs(j)
q1 = qs(j+1)

r   = rs(j+1)
rhs = rs(j+1) + h/2 * rs(j+1)**2 + h/2 * q1 + h/2*q0

do iter=1,niters
delta = (rhs - r +h/2*r**2 ) / (1 - h*r)
r     = r + delta
end do

rs(j)    = r
rders(j) = -r**2 - q0
end do


end subroutine


subroutine rriccati_linear_tvp(k,a,b,ts,chebint,qs,fs,rb,amatr,rs,rders)
implicit double precision (a-h,o-z)
integer                  :: k
double precision         :: ts(k),chebint(k,k)
double precision         :: amatr(k,k),qs(k),fs(k),rs(k),rders(k),rb
!
!  Use a Chebyshev spectral method to solve the terminal value problem
!
!    r'(t) + qs(t) r(t) = f(t)
!    r(b)               = rb
!
!  over the interval (a,b).
!
!  Input parameters:
!    k - the number of Chebyshev nodes
!    ts - the k Chebyshev nodes on the interval (a,b)
!    chebint - the "right" Chebyshev spectral integration matrix
!    qs - the values of the coefficient at the Chebyshev nodes
!    ra - the initial value for the solution
!
!  Output parameters:
!    rs - the values of the solution at the  Chebyshevnodes
!    rders - the values of the derivative of the solution at the 
!      Chebyshev nodes
!
!  Work arrays:
!    amatr - a (k,k) matrix 
!



!
!  We solve the integral equation 
!
!     \sigma(t) + q(t) \int_b^t \sigma(s) ds = f(t) - rb q(t)
!
!  obtained by letting
!
!     r(t) = rb + \int_b^t \sigma(s) ds.
! 


rders = 0
amatr = 0

do i=1,k
amatr(i,i) = 1.0d0
end do

do i=1,k
rders(i)   = fs(i) - rb * qs(i)
amatr(i,:) = amatr(i,:) + qs(i) * chebint(i,:)*(b-a)/2
end do


call qrsolv(amatr,k,rders)
rs    = rb + (b-a)/2*matmul(chebint,rders)

end subroutine


subroutine rriccati_int_back(chebdata,nints,ab,rs,rders,psival,psis)
implicit double precision (a-h,o-z)
type(chebexps_data)                             :: chebdata
integer                                         :: nints
double precision, allocatable, intent(in)       :: ab(:,:)
double precision, allocatable, intent(in)       :: rs(:,:),rders(:,:)
double precision, allocatable, intent(out)      :: psis(:,:)
double precision                                :: psival
!
!  Integrate the solution of (1) backward; that is, construct the function
!
!                          t
!    psi(t) = psival + \int r(u) du 
!                          b
!
!  Input parameters:
!    chebdata - the structure returned by chebexps
!    nints - the number of subintervals used to discretize the solution
!    ab - a (2,nints) array specifying the endpoints of said subintervals
!    rs - a (k,nints) array giving the values of the solution at the Chebyshev
!       nodes on each subinterval
!    rders - a (k,nints) array giving the values of the derivatives of the
!       solution at the Chebyshev nodes on each subinterval
!    psival - the value of psi at a
!
!  Output parameters:
!    psis  the values of the function psi at the Chebyshev nodes


!double complex                  :: psival0


k = chebdata%k
allocate(psis(k,nints))

psival0 = psival

do int = nints,1,-1

a = ab(1,int)
b = ab(2,int)

psis(:,int) = psival0 + (b-a)/2 * matmul(chebdata%aintr,rs(:,int))
psival0     = psis(1,int)
end do

end subroutine


end module
