module test_pbessel_functions

use utils
use chebyshev
use iso_c_binding

integer          :: iwhich

contains

subroutine qfun(ising,n,rs,qs,userptr)
implicit double precision (a-h,o-z)
double precision          :: rs(n), qs(n), qders(n)
type(c_ptr)               :: userptr
!qs = rs**2-1 

qs = 0
if (ising .eq. 1) qs = 1.0d0
if (ising .eq. 3) qs = 2.0d0

end subroutine



end module


program test_pbessel
use utils
use chebyshev
use pbessel
use iso_c_binding
use test_pbessel_functions
implicit double precision (a-h,o-z)

type(chebexps_data)                        :: chebdata
type(c_ptr)                                :: userptr
type(pbessel_solution)                     :: soldata,soldata2
double precision, allocatable              :: xsings(:)
double precision, allocatable              :: xs(:),ys(:)

k = 30
call chebexps(k,chebdata)

eps0    =    epsilon(0.0d0)
eps     =    1.0d-12
pi      =    acos(-1.0d0)
R       =    4

dlambda = 2.0d0**17
dnu     = 454545

allocate(xsings(3))
xsings(1) = 1
xsings(2) = 2
xsings(3) = 3

call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata)

stop


!
!  Perform an artifical test using the fact that 
!
!    { J_(n/2) ( lambda/2 * t^2) sqrt(t),  Y_(n/2) ( lambda/2 * t^2) sqrt(t) }
!
!  is a basis in the space of solutions of the differential equation
!
!    y''(t) + (dlambda^2 t^2 + (1/4-n^2)/t^2) y(t) = 0
!
!

allocate(xsings(0))

eps    = 1.0d-12
R      = 1
pi     = acos(-1.0d0)

write (*,"(A10,3X,A10,3X,A10,3XA10)")       "wavenumber","time","ratio","max error"
write (*,*)                                 ""

write (13,"(A10,3X,A10,3X,A10,3XA10)")       "wavenumber","time","ratio","max error"
write (13,*)                                 ""

dtime  = 0

jj1    = 6
jj2    = 20

do jj=jj1,jj2
 
dtime2  = dtime
dlambda = 2**jj
m1      = 0
m2      = dlambda
dmax    = 0

call elapsed(t1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,dnu,soldata,t,val,der,val0,cc,derr,ier,t0,t1,eps1)
!$OMP DO
do n = m1,m2,2
dnu = n

call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata)

t     = 1.0d0
call pbessel_eval(soldata,t,val,der)
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)
cc    = val0 / val

t     = 0.9d0
call pbessel_eval(soldata,t,val,der)
val   = cc*val
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)
derr  = abs(val-val0)

!$OMP ATOMIC
dmax = max(derr,dmax)

end do

!$OMP END DO
!$OMP END PARALLEL

call elapsed(t2)
dtime = t2-t1

if (jj .eq. jj1) dtime00 = dtime

lambda = dlambda
if (dtime2 .ne. 0) then
dratio= dtime/dtime2
else
dratio = -1
endif

est = dtime00 * 2**(jj-jj1)

write (*,"(I10,3X,D10.3,3X,D10.3,3X,D10.3,3X,D10.3)")   lambda,dtime,dratio,dmax
write (13,"(I10,3X,D10.3,3X,D10.3,3X,D10.3,3X,D10.3)")   lambda,dtime,dratio,dmax



dtime2 = dtime

end do

end program
