!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for finding a solution of the differential equation
!
!              (                             1/4 - dnu^2 )
!    u''(t) +  ( dlambda^2 ( 1 + q(t) ) +     ---------  )  u(t) = 0                     (1)
!              (                                t^2      )
!
!  on the interval (0,R) such that
!
!      u(t) = O( t^(dnu+1/2) ) as t -> 0.                                                (2)
!
!  Because just one boundary condition is being imposed, the solution is only
!  determined up to a multiplicative constant.
!
!  Since  (1) becomes the normal form Bessel's differential equation when q(t) = 0,
!  we refer to it as the perturbed Bessel equation.  In this case, the solution is  
!  a multiple of the function J_dnu(dlambda t) \sqrt{t} .
!
!  The coefficient q(t) is supplied by the user via an external subroutine and it is
!  assumed to be piecewise smooth and >= -1 in the interval (0,R).  The points 
!
!    x0 < x1 < x2 < ... < xn
!
!  at which q or its derivatives are discontinuous must be given as inputs to the 
!  solver.    Because this solver uses discretization grids which include the
!  endpoints of the interval under consideration, we use a kludge toindicate which
!  of th subinterval contains the points at which q is to be evaluated.  More
!  explicitly, the external subroutine provided by the user is supplied with
!  the index of the interval
! 
!    [0,x1], [x1,x2], [x2,x3], ..., [xn, R]
!
!  in which the points at which is being asked to evaluate q are located.  
!
!
!  The solution of (1) is represented (more-or-les
!  On regions in which the sign of the coefficient
!
!    Q(t) = dlambda^2 q(t)  + (1/4 - dnu^2) / t^2                                         (3)
!
!  is positive, the solution of (1) is represented indirectly, via 
!  That is, in the form
!
!        sin(alpha(t) )              cos( alpha(t) ) 
!   a1  ---------------------  +  a2 ------------------
!         sqrt(alpha'(t)             sqrt (alpha'(t) )
!
!  with alpha(t) a smooth, nonoscillatory function.  In the regions in which Q(t) is
!  negative, the solution of (1) is represented via the
!  form
!
!     b1 exp(psi(t))
!
!  with psi(t) a real-valued function and b1 either equal to either 1 or -1.  
!
!  The running time of this solver is essentially independent of the parameters
!  dnu and dlambda.
!
!  The following subroutines should be regarded as publicly callable:
!
!    pbessel_solve - given an external subroutine for evaluating q(t), construct a 
!     solution of (1) which satisfies the boundary condition (2)
!
!    pbessel_eval - evaluate a solution of (1) constructed by pbessel_solve (and
!     its derivative) at a specified point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pbessel

use utils
use chebyshev
use riccati
use iso_c_binding

interface

subroutine pbessel_qfun(ising,n,ts,qs,userptr)
use iso_c_binding
integer          :: ising
double precision :: ts(n),qs(n)
type(c_ptr)      :: userptr
!
!  This defines the interface for the external function specifying the coefficient q
!  in (1).  The the value of the function q(t) in (1) at each of the n points specified 
!  by the ts array should be returned in the qs array.  The integer parameter ising
!  gives the "singularity" interval in which the points ts are contained.
!
end subroutine
end interface


type        pbessel_riccati_qfun_data
double precision                         :: dnu
double precision                         :: dlambda
double precision                         :: R
double precision                         :: a
double precision                         :: b
double precision                         :: deta
integer                                  :: ising
procedure(pbessel_qfun), pointer, nopass :: qfun
type(c_ptr)                              :: userptr
end type    pbessel_riccati_qfun_data


! The pbessel_interval structure describes the solution of (1) over a single interval
! on which it is either oscillatory or nonoscillatory
type        pbessel_interval
integer                                  :: ifosc       ! indicates whether we are in the
                                                        ! oscillatory regime (ifosc = 1)
integer                                  :: ifleft      ! indicates whether this is the
                                                        ! leftmost interval (ifleft = 1) 
                                                        ! or not  (ifleft = 0)
integer                                  :: k           ! size of the Chebyshev expansions
                                                        ! used to represent functions

integer                                  :: ising

! redundant copies of the parameters, mainly for debugging purposes
double precision                         :: a
double precision                         :: b
double precision                         :: dnu
double precision                         :: dlambda

!  data for the oscillatory case -- the solution is represented via a phase function
integer                                  :: nints_a
double precision, allocatable            :: ab_a(:,:)
double precision, allocatable            :: acoefs(:,:)
double precision, allocatable            :: apcoefs(:,:)
double precision, allocatable            :: appcoefs(:,:)

! data for the nonoscillatory case -- the solution is represented via the logarithms
! of two solutions, one exponentially increasing and other exponentially decreasing

integer                                  :: nints_psi1
double precision, allocatable            :: ab_psi1(:,:)
double precision, allocatable            :: psi1coefs(:,:)
double precision, allocatable            :: psi1pcoefs(:,:)
double precision, allocatable            :: psi1ppcoefs(:,:)

integer                                  :: nints_psi2
double precision, allocatable            :: ab_psi2(:,:)
double precision, allocatable            :: psi2coefs(:,:)
double precision, allocatable            :: psi2pcoefs(:,:)
double precision, allocatable            :: psi2ppcoefs(:,:)


! the coefficients which specify which linear cominbation of the two solutions is
! the desired solution

real*16                                  :: a1, a2

end type    pbessel_interval

! The pbessel_solution structure describes a solution of (1)
type        pbessel_solution
double precision                         :: dnu
double precision                         :: dlambda
integer                                  :: nints
integer                                  :: ifnull 
double precision, allocatable            :: ab(:,:)
type(pbessel_interval), allocatable      :: sols(:)
end type    pbessel_solution


contains



subroutine pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata)
implicit double precision (a-h,o-z)
type(chebexps_data)                         :: chebdata
type(pbessel_solution), intent(out)         :: soldata
type(c_ptr)                                 :: userptr
procedure(pbessel_qfun)                     :: qfun
double precision, intent(in)                :: xsings(:)
!
!  Construct a solution of the ODE (1) which satisfies the boundary condition (2).
!  The user supplies the function q(t) in the coefficient of (1) via an external
!  subroutine.
!
!  IMPORTANT NOTE - if this routine encounters an error, it will halt program
!  execution with an appropriate message indicating what has gone wrong
!  rather than return an error code
!
!  Input parameters:
!    eps - the desired accuracy for the solution
!    chebdata - data structure generated by chebexps which, among other things,
!     specifies the order of the piecewise chebyshev expansion used to
!     represent the solution through its element "k"
!    xsings - an array which gives the list of points at which the input
!      function q(t) is singular
!    R - the radius of the interval on which (1) is to be solved
!    dlambda - the value of the parameter lambda in (1)
!    dnu - the value of the parameter nu in (1)
!    qfun - a user-supplied external subroutine conforming to the pbessel_qfun 
!      interface which supplies the values of q(t) and its derivative q'(t)
!    userptr - a "void *" pointer which is passed to the user-supplied external
!      subroutine qfun
!
!  Output parameters:
!    soldata - a data structure describing the obtained solution of (1)
!

double precision, allocatable    :: xsings0(:),roots0(:)
double precision, allocatable    :: ab(:,:)
double precision, allocatable    :: ts(:),qvals(:)
integer, allocatable             :: isings(:)

real*16                          :: uval1,uder1,vval1,vder1
real*16                          :: uval2,uder2,vval2,vder2,det
real*16                          :: a1,a2,c1,c2,b1,b2

!
!  The precomputed grids are only accurate to about 13 digits of precision ...
!  use the adaptive procedure when more accuracy is required
!

ier        = 0
eps0       = epsilon(0.0d0)
pi         = acos(-1.0d0)
k          = chebdata%k

allocate ( ts(1),qvals(1) )

nsings = size(xsings)
allocate ( xsings0(nsings+3) )
allocate ( ab(2,1000) )
allocate ( isings(1000) )

nints   = 0

! a = 1.0d-7

dd = -300
a   = exp( ( dd + dnu*log(2.0d0) - dnu*log(dlambda) +  log_gamma(1+dnu) ) / (dnu+0.5d0) ) 
a   = min(1.0d-2,a)
a   = max(1.0d-7,a)
b   = R

if (dnu .eq. 0)                        a = 1.0d-12
if (eps .lt. 1.0d-16 .AND. dnu .eq. 0) a = 1.0d-15

nn                  = nsings+1
xsings0(1)          = a
xsings0(2:nsings+1) = xsings
xsings0(nsings+2)   = b

do i=1,nn
a0 = xsings0(i)
b0 = xsings0(i+1)

ifleft = 0
if (i .eq. 1) ifleft = 1
call pbessel_roots(i,ifleft,a0,b0,dlambda,dnu,qfun,userptr,nroots0,roots0)


if (nroots0 .ge. 1) then
if ( abs(roots0(nroots0) - b0 ) .lt. 1.0d-6*b0) then
nroots0=nroots0-1
endif
endif

if (nroots0 .ge. 1) then
if (abs(roots0(1) - a0) .lt. 1.0d-6*a0) then
nroots0=nroots0-1
roots0(1:nroots0-1) = roots0(2:nroots0)
endif
endif

nints         = nints+1
ab(1,nints)   = a0
isings(nints) = i

do j=1,nroots0
ab(2,nints)   = roots0(j)
nints         = nints+1
ab(1,nints)   = roots0(j)
isings(nints) = i
end do
ab(2,nints) = b0

end do

!
!  If the last interval is too small and oscillatory, the windowing
!  algorithm can run into trouble.  Best to simply make it bigger if need be.
!
if (ab(2,nints) - ab(1,nints) .lt. 0.1d0) then
ab(2,nints) = ab(1,nints) + 0.1d0
endif

!
! Traverse the intervals, solving the problem on each interval
!

soldata%nints = nints
allocate(soldata%ab(2,nints))
allocate(soldata%sols(nints))
soldata%ab = ab(:,1:nints)

do int=1,nints
a0    = ab(1,int)
b0    = ab(2,int)
ising = isings(int)

ts(1) = (a0+b0)/2
call qfun(ising,1,ts,qvals,userptr)
qvals = dlambda**2*(qvals+1) + (0.25d0 - dnu**2)/ts**2

ifosc = 0
if (qvals(1) > 0) ifosc = 1


! indicate to the solvers whether we are on the left-most interval or not; this
! influences the choice of discretization grid in the nonoscillatory case

ifleft  = 0
if (int .eq. 1)     ifleft = 1

if (ifosc .eq. 0) then
call pbessel_nonoscillatory(ising,ifleft,eps,a0,b0,chebdata,dlambda,dnu,qfun,userptr, &
  soldata%sols(int))

else

call pbessel_oscillatory(ising,ifleft,eps,a0,b0,chebdata,dlambda,dnu,qfun,userptr, &
  soldata%sols(int))

endif
end do

!
!  Find the appropriate values for the left interval
!

do int=1,nints
a0     = ab(1,int)
b0     = ab(2,int)

ifleft = soldata%sols(int)%ifleft
ifosc  = soldata%sols(int)%ifosc


if (ifleft .eq. 1 .AND. ifosc .eq. 0) then
a1 = 1
a2 = 0
goto 1000
endif

if (ifleft .eq. 1) then
nn  = dnu
c1  = sqrt(pi/2)*bessel_jn(nn,dlambda*a0)*sqrt(a0)
c2  = sqrt(pi/2)*&
  (2*a0*dlambda * bessel_jn(nn-1,dlambda*a0) + (1-2*dnu)*bessel_jn(nn,dlambda*a0)) / (2*sqrt(a0))

else


b1 = soldata%sols(int-1)%a1
b2 = soldata%sols(int-1)%a2


call pbessel_eval_interval2(soldata%sols(int-1),a0,uval1,uder1,vval1,vder1)

c1 = b1*uval1 + b2*vval1
c2 = b1*uder1 + b2*vder1


endif

call pbessel_eval_interval2(soldata%sols(int),a0,uval2,uder2,vval2,vder2)



det = uval2*vder2 - uder2 * vval2


if (abs(det) .ge. 1.0d-100) then

a1  = (vder2  * c1 - vval2 * c2)/det
a2  = (-uder2 * c1 + uval2 * c2)/det

else

a1 = 1
a2 = 0


do ii=int-1,1,-1
soldata%sols(ii)%a1 = 0
soldata%sols(ii)%a2 = 0
end do

endif



1000 continue

soldata%sols(int)%a1 = a1
soldata%sols(int)%a2 = a2



end do

!
!  Normalize the solution at t = R if it is large there!
!

a1 = soldata%sols(nints)%a1 
a2 = soldata%sols(nints)%a2 

call pbessel_eval_interval2(soldata%sols(nints),R,uval1,uder1,vval1,vder1)
c1 = a1 * uval1 + a2*vval1
!c1 = abs(c1)


if (abs(c1) .gt. 1.0d0) then

do int=1,nints

soldata%sols(int)%a1 = soldata%sols(int)%a1 / c1
soldata%sols(int)%a2 = soldata%sols(int)%a2 / c1



end do

endif


end subroutine


subroutine pbessel_eval(soldata,t,val,der)
implicit double precision (a-h,o-z)
type(pbessel_solution), intent(in)         :: soldata
!
!  Evaluate the solution of (1) obtained by pbessel_sol and its derivative at a 
!  specified point.
!
!  Input parameters:
!    eps - the desired accuracy for the solution
!    chebdata - data structure generated by chebexps which, among other things,
!     specifies the order of the piecewise chebyshev expansion used to
!     represent the solution through its element "k"
!    R - the radius of the interval on which (1) is to be solved
!    dlambda - the value of the parameter lambda in (1)
!    dnu - the value of the parameter nu in (1)
!    qfun - a user-supplied external subroutine conforming to the pbessel_qfun 
!      interface which supplies the values of q(t) and its derivative q'(t)
!    userptr - a "void *" pointer which is passed to the user-supplied external
!      subroutine qfun
!
!  Output parameters:
!    ier - an error return code;
!        ier = 0    indicates successful execution
!        ier = 4    indicates the riccati_tvp procedure failed
!        ier = 8    indicates the riccati_ivp procedure failed
!        ier = 1024 means that Q(t) appears to have multiple turning points
!
!    pbdata - a data structure describing a solution of (1)
!


if (t .lt. soldata%ab(1,1)) then
val = 0
der = 0
return
endif


nints = soldata%nints

do int=1,nints-1
if (t .le. soldata%ab(2,int)) exit
end do

a  = soldata%ab(1,int)
b  = soldata%ab(2,int)

a1 = soldata%sols(int)%a1
a2 = soldata%sols(int)%a2

call  pbessel_eval_interval(soldata%sols(int),t,uval,uder,vval,vder)


val = a1*uval+a2*vval
der = a1*uder+a2*vder


end subroutine



subroutine pbessel_oscillatory(ising,ifleft,eps,a,b,chebdata,dlambda,dnu,qfun,userptr,intdata)
implicit double precision (a-h,o-z)
type(chebexps_data)                         :: chebdata
type(pbessel_interval), intent(out)         :: intdata
type(c_ptr)                                 :: userptr
procedure(pbessel_qfun)                     :: qfun
!
!  Construct a basis in the space of solutions of the differential equation (1) 
!  over an interval (a,b) on which the coefficient Q(t) in (3) is positive.
!  The solutions are represented via a nonoscillatory phase function.
!
!  Input parameters:
!    ifleft - indiciate whether or not this is the leftmost interval on
!      which (1) is being solved
!    eps - the precision for the solution
!    (a,b) - the interval on which the solution is to be constructed
!    chebdata - data structure generated by chebexps which, among other things,
!     specifies the order of the piecewise chebyshev expansion used to
!     represent the solution through its element "k"
!    dlambda - the value of the parameter dlambda in (1)
!    dnu - the value of the parameter dnu in (1)
!    qfun - a user-supplied external function conforming to the interface
!      pbessel_qfun which supplies the values of q(t) in (1)
!    userptr - a  "void *" pointer which is passed by this routine to qfun
!
!  Output parameters:
!    intdata - a data structure describing the solution
!

type(pbessel_riccati_qfun_data), pointer    :: rdata
type(c_ptr)                                 :: rdataptr
double precision, allocatable               :: ts(:),qvals(:),ab(:,:),ab0(:,:)
double precision, allocatable               :: vals(:),coefs(:)

! data for the phase function
double complex, allocatable                 :: rs(:,:),rders(:,:)
double precision, allocatable               :: avals(:,:),apvals(:,:),appvals(:,:)

! data for the nonoscillatory regime
double complex, allocatable                 :: rs2(:,:),rders2(:,:),psis(:,:)
double complex                              :: ima, rb, ra

ier        = 0
eps0       = epsilon(0.0d0)
ima        = (0.0d0,1.0d0)
pi         = acos(-1.0d0)
k          = chebdata%k

!
!  Use the windowing algorithm to construct a phase function
!

deta = dlambda

allocate(rdata)
intdata%k      = k
rdata%dnu      = dnu
rdata%dlambda  = dlambda
rdata%deta     = deta
rdata%ising    = ising
rdata%a        = a 
rdata%b        = b
rdata%qfun    => qfun
rdata%userptr  = userptr
rdataptr       = c_loc(rdata)


rb = ima*deta
call riccati_tvp_adap(ier,eps,chebdata,a,b,pbessel_riccati_qfun2,rb,nints,ab, &
 rs,rders,rdataptr)

if  (ier .ne. 0) then
call prini("in pbessel_oscillatory, riccati_tvp_adap failed with ier = ",ier)
call prind("eps = ",eps)
call prin2("dlambda = ",dlambda)
call prin2("dnu = ",dnu)
call prin2("a = ",a)
call prin2("b = ",b)
stop
endif


ra = rs(1,1)
call riccati_ivp_adap(ier,eps,chebdata,a,b,pbessel_riccati_qfun3,ra,nints,ab, &
   rs,rders,rdataptr)
if  (ier .ne. 0) then
call prini("in pbessel_oscillatory, riccati_ivp_adap failed with ier = ",ier)
call prind("eps = ",eps)
call prin2("dlambda = ",dlambda)
call prin2("dnu = ",dnu)
call prin2("a = ",a)
call prin2("b = ",b)
stop
endif

aval    =   pi/4
call riccati_phase_forward(chebdata,nints,ab,rs,rders,aval,avals,apvals,appvals)



!
!  Copy the data into the output structure
!

intdata%ifosc   = 1
intdata%ifleft  = ifleft
intdata%a       = a
intdata%b       = b
intdata%dlambda = dlambda
intdata%dnu     = dnu


allocate(intdata%ab_a(2,nints))
allocate(intdata%acoefs(k,nints),intdata%apcoefs(k,nints),intdata%appcoefs(k,nints))

intdata%nints_a = nints
intdata%ab_a    = ab

do int=1,nints
intdata%acoefs (:,int)    = matmul(chebdata%u,avals(:,int))
intdata%apcoefs (:,int)   = matmul(chebdata%u,apvals(:,int))
intdata%appcoefs (:,int)  = matmul(chebdata%u,appvals(:,int))
end do

end subroutine



subroutine pbessel_nonoscillatory(ising,ifleft,eps,a,b,chebdata,dlambda,dnu,qfun,userptr,intdata)
implicit double precision (a-h,o-z)
type(chebexps_data)                         :: chebdata
type(pbessel_interval), intent(out)         :: intdata
type(c_ptr)                                 :: userptr
procedure(pbessel_qfun)                     :: qfun

!
!  Build a basis in the space of solutions of (1) on an interval (a,b) in
!  which the coefficient Q(t) in (3) is negative.
!

type(pbessel_riccati_qfun_data), pointer    :: rdata
type(c_ptr)                                 :: rdataptr
double precision, allocatable               :: ts(:),qvals(:),vals(:),coefs(:)

double precision, allocatable               :: ab(:,:),psi(:,:),psip(:,:),psipp(:,:)
double precision, allocatable               :: ab2(:,:),psi2(:,:),psi2p(:,:),psi2pp(:,:)

double complex                              :: ima

ier        = 0
eps0       = epsilon(0.0d0)
ima        = (0.0d0,1.0d0)
pi         = acos(-1.0d0)
k          = chebdata%k

allocate(rdata)
rdata%dnu      =  dnu
rdata%dlambda  =  dlambda
rdata%ising    = ising
rdata%qfun     => qfun
rdata%userptr  =  userptr
rdataptr       =  c_loc(rdata)

!
!  Solve the Riccati equation forward to find the logarithm of a function
!  y(t) which behaves like an increasing exponential
!

if (ifleft .eq. 1) then 
ra = (dnu+0.5d0)/a
else
ra = 0
endif

call rriccati_ivp_adap(ier,eps,chebdata,a,b,pbessel_riccati_qfun1,ra,&
  nints,ab,psip,psipp,rdataptr)

if (ier .ne. 0) then
call prini("in pbessel_nonoscillatory, rrriccati_ivp_adap ier = ",ier)
call prin2("dlambda = ",dlambda)
call prin2("dnu = ",dnu)
call prin2("a = ",a)
call prin2("b = ",b)
stop
endif

rb = 0
call rriccati_int_back(chebdata,nints,ab,psip,psipp,rb,psi)

rb = 0
call rriccati_tvp_adap(ier,eps,chebdata,a,b,pbessel_riccati_qfun1,rb,&
  nints2,ab2,psi2p,psi2pp,rdataptr)

if (ier .ne. 0) then
call prini("in pbessel_nonoscillatory, rrriccati_tvp_adap ier = ",ier)
call prin2("dlambda = ",dlambda)
call prin2("dnu = ",dnu)
call prin2("a = ",a)
call prin2("b = ",b)
stop
endif

ra = 0
call rriccati_int_forward(chebdata,nints2,ab2,psi2p,psi2pp,ra,psi2)

!
!  Copy the data into the output structure
!

intdata%ifosc   = 0
intdata%ifleft  = ifleft
intdata%dnu     = dnu
intdata%dlambda = dlambda
intdata%a       = a
intdata%b       = b
intdata%k       = k

intdata%nints_psi1 = nints
intdata%nints_psi2 = nints2

allocate(intdata%ab_psi1(2,nints))
allocate(intdata%ab_psi2(2,nints2))

intdata%ab_psi1 = ab
intdata%ab_psi2 = ab2

allocate(intdata%psi1coefs(k,nints))
allocate(intdata%psi1pcoefs(k,nints))
allocate(intdata%psi1ppcoefs(k,nints))

allocate(intdata%psi2coefs(k,nints2))
allocate(intdata%psi2pcoefs(k,nints2))
allocate(intdata%psi2ppcoefs(k,nints2))


allocate(vals(k))

do int=1,nints
vals                          = psi(:,int)
intdata%psi1coefs (:,int)     = matmul(chebdata%u,vals)

vals                          = psip(:,int)
intdata%psi1pcoefs (:,int)    = matmul(chebdata%u,vals)

vals                          = psipp(:,int)
intdata%psi1ppcoefs (:,int)   = matmul(chebdata%u,vals)
end do


do int=1,nints2
vals                          = psi2(:,int)
intdata%psi2coefs (:,int)     = matmul(chebdata%u,vals)

vals                          = psi2p(:,int)
intdata%psi2pcoefs (:,int)    = matmul(chebdata%u,vals)

vals                          = psi2pp(:,int)
intdata%psi2ppcoefs (:,int)   = matmul(chebdata%u,vals)

end do

end subroutine




subroutine pbessel_eval_interval(intdata,t,uval,uder,vval,vder)
implicit double precision (a-h,o-z)
type(pbessel_interval), intent(in)         :: intdata
!
!
!


ifosc = intdata%ifosc
a1    = intdata%a1
a2    = intdata%a2
k     = intdata%k


if (ifosc .eq. 1) then
call chebpw_eval23(intdata%nints_a,intdata%ab_a,k,intdata%acoefs,intdata%apcoefs,&
  intdata%appcoefs,t,aval,apval,appval)

uval = sin(aval)/sqrt(apval)
uder = cos(aval)*sqrt(apval) - sin(aval)*appval/(2*apval**1.5d0)

vval = cos(aval)/sqrt(apval)
vder =-sin(aval)*sqrt(apval) - cos(aval)*appval/(2*apval**1.5d0)


else

call chebpw_eval22(intdata%nints_psi1,intdata%ab_psi1,k,intdata%psi1coefs,intdata%psi1pcoefs, &
  t,psi1,psi1p)

call chebpw_eval22(intdata%nints_psi2,intdata%ab_psi2,k,intdata%psi2coefs,intdata%psi2pcoefs, &
  t,psi2,psi2p)



uval = exp(psi1)
uder = exp(psi1)*psi1p

vval = exp(psi2)
vder = exp(psi2)*psi2p 


endif

end subroutine





subroutine pbessel_eval_interval2(intdata,t,uval,uder,vval,vder)
implicit double precision (a-h,o-z)
type(pbessel_interval), intent(in)         :: intdata
real*16                                    :: uval,uder
real*16                                    :: vval,vder

!
!
!

real*16                                    :: x1,x2,a1,a2

ifosc = intdata%ifosc
a1    = intdata%a1
a2    = intdata%a2
a     = intdata%a
b     = intdata%b
k     = intdata%k


if (ifosc .eq. 1) then
call chebpw_eval23(intdata%nints_a,intdata%ab_a,k,intdata%acoefs,intdata%apcoefs,&
  intdata%appcoefs,t,aval,apval,appval)


uval = sin(aval)/sqrt(apval)
uder = cos(aval)*sqrt(apval) - sin(aval)*appval/(2*apval**1.5d0)

vval = cos(aval)/sqrt(apval)
vder =-sin(aval)*sqrt(apval) - cos(aval)*appval/(2*apval**1.5d0)


else

call chebpw_eval22(intdata%nints_psi1,intdata%ab_psi1,k,intdata%psi1coefs,intdata%psi1pcoefs, &
  t,psi1,psi1p)

call chebpw_eval22(intdata%nints_psi2,intdata%ab_psi2,k,intdata%psi2coefs,intdata%psi2pcoefs, &
  t,psi2,psi2p)



x1   = psi1
x2   = psi1p

uval = exp(x1)
uder = exp(x1)*x2

x1   = psi2
x2   = psi2p

vval = exp(x1)
vder = exp(x1)*x2

!print *,psi1,psi1p,psi2,psi2p

endif

end subroutine


subroutine pbessel_riccati_qfun1(k,ts,qs,userptr)
implicit double precision (a-h,o-z)
integer                                  :: k
double precision                         :: ts(k)
double precision                         :: qs(k)
type(c_ptr)                              :: userptr
type(pbessel_riccati_qfun_data), pointer :: rdata
double precision                         :: qs0(k)

call c_f_pointer(userptr,rdata)
dnu     = rdata%dnu
dlambda = rdata%dlambda
ising   = rdata%ising

call rdata%qfun(ising,k,ts,qs0,rdata%userptr)
qs = dlambda**2*(qs0+1) + (0.25d0 - dnu**2)/ ts**2

end subroutine


subroutine pbessel_riccati_qfun2(k,ts,qs,userptr)
implicit double precision (a-h,o-z)
integer                                  :: k
double precision                         :: ts(k)
double complex                           :: qs(k)
type(c_ptr)                              :: userptr
type(pbessel_riccati_qfun_data), pointer :: rdata

double precision                         :: qs1(k),qs2(k)
double precision                         :: phis(k)

call c_f_pointer(userptr,rdata)

dnu      = rdata%dnu
dlambda  = rdata%dlambda
deta     = rdata%deta    
a        = rdata%a
b        = rdata%b
ising    = rdata%ising

c        = (a+b)/2
dd       = -18.0d0/(b-a)

call rdata%qfun(ising,k,ts,qs1,rdata%userptr)
qs1      = dlambda**2*(qs1+1) + (0.25d0 - dnu**2)/ ts**2
qs2      = deta**2

phis     = (1-erf(dd*(ts-c)))/2.0d0
qs       = (1-phis)*qs1 + (phis)*qs2


end subroutine


subroutine pbessel_riccati_qfun3(k,ts,qs,userptr)
implicit double precision (a-h,o-z)
integer                                  :: k
double precision                         :: ts(k)
double complex                           :: qs(k)
type(c_ptr)                              :: userptr
type(pbessel_riccati_qfun_data), pointer :: rdata
double precision                         :: qs1(k)

call c_f_pointer(userptr,rdata)

dnu     = rdata%dnu
dlambda = rdata%dlambda
ising   = rdata%ising

call rdata%qfun(ising,k,ts,qs1,rdata%userptr)
qs = dlambda**2*(qs1+1) + (0.25d0 - dnu**2)/ ts**2

end subroutine



subroutine pbessel_roots(ising,ifleft,a,b,dlambda,dnu,qfun,userptr,nroots,roots)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: roots(:)
type(c_ptr)                                 :: userptr
procedure(pbessel_qfun)                     :: qfun

!
!  Use a very simple procedure to attempt to find all of the roots
!  of the coefficient Q(t) = dlambda^2 ( 1+ q(t) ) + (1/4-dnu^2)/t^2
!  on a specified interval.
!
!  This procedure is relatively fast but could potentially fail;
!  a foolproof version could be constructed with the chebpw_roots
!  code.
!

double precision, allocatable  :: ts(:),vals(:)
double precision, allocatable  :: ab(:,:)

nbisect  = 200
maxroots = 1000
nn       = 4*maxroots

allocate(vals(nn),ts(nn),ab(2,maxroots))
ab = 0
ts = 0


if (ifleft .eq. 1) then

m0 = 20
bb = a + (b-a)*.1d0
do i=1,m0
dd    = 2.0d0**(-m0+i)
ts(i) = a + (bb-a)*dd
end do


do i=m0+1,nn
dd    = (i-m0-1.0d0) / (nn-m+00.0d0)
ts(i) = bb + (b-bb)*dd
end do

ts(nn) = b-1.0d-12


else

do i=1,nn
dd    = (i+0.0d0)/(nn+1.0d0)
ts(i) = a + (b-a)*dd
end do

ts(nn) = b-1.0d-12

endif


call qfun(ising,nn,ts,vals,userptr)
vals = dlambda**2*(vals+1) + (0.250 - dnu**2)/ts**2

m = 0
do i=1,nn-1

ifint = 0

if (vals(i)*vals(i+1) < 0) then
m = m +1

if (m .gt. maxroots) then
call prina("pbessel_find_roots failed: too many roots")
stop
endif

ab(1,m) = ts(i)
ab(2,m) = ts(i+1)
endif
end do


nroots = m
allocate(roots(nroots))


!
!  Use bisection to refine the brackets 
!


do iter=1,nbisect
do int=1,m

a0 = ab(1,int)
b0 = ab(2,int)
c0 = (a0+b0)/2

ts(1) = a0
ts(2) = c0

call qfun(ising,2,ts,vals,userptr)
vals(1:2) = dlambda**2*(vals(1:2)+1) + (0.250 - dnu**2)/ts(1:2)**2

val1 = vals(1)
val2 = vals(2)

if (val1*val2 < 0) then
ab(2,int) = c0
else
ab(1,int) = c0
endif

end do
end do

!
!  Copy the roots out
!
do int=1,m
a0         = ab(1,int)
b0         = ab(2,int)
t          = (a0+b0)/2
roots(int) = t
end do


end subroutine



end module
