!
!  Perform various experiments measuring the performance of the solver for the
!  perturbed Bessel equation using the fact that
!
!    { J_(n/2) ( lambda/2 * t^2) sqrt(t),  Y_(n/2) ( lambda/2 * t^2) sqrt(t) }
!
!  is a basis in the space of solutions of the differential equation
!
!    y''(t) + (dlambda^2 t^2 + (1/4-n^2)/t^2) y(t) = 0.
!
!
!  In this test, the wavenumber is fixed and the cost of solving the 
!

module test_pbessel_functions

use utils
use chebyshev
use iso_c_binding

double precision :: dnu, dlambda

contains

subroutine qfun(ising,n,rs,qs,userptr)
implicit double precision (a-h,o-z)
double precision          :: rs(n), qs(n)
type(c_ptr)               :: userptr
integer                   :: ising

qs = rs**2-1

end subroutine


! subroutine fit_constant(nn,xs,ys,c)
! implicit double precision (a-h,o-z)

! double precision              ::  xs(nn),ys(nn),dd(1,1)
! double precision, allocatable :: amatr(:,:),bmatr(:,:),zs(:,:)

! allocate(amatr(nn,1),bmatr(1,nn),zs(nn,1))

! do i=1,nn
! amatr(i,1)  = log(xs(i))
! bmatr(1,i)  = amatr(i,1)
! zs(i,1)     = ys(i)
! end do

! dd  = matmul(bmatr,amatr)
! dd  = matmul(bmatr,zs)/dd(1,1)

! c = dd(1,1)


! end subroutine


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
type(pbessel_solution), allocatable        :: soldata(:)
type(pbessel_solution)                     :: soldata0
double precision, allocatable              :: xsings(:)
double precision, allocatable              :: xs(:),ys(:)
double precision, allocatable              :: derrs(:),dtimes(:),dnus(:),xs0(:),ys0(:),dlog(:),dlog2(:)
double precision, allocatable              :: dlambdas(:)

double precision, allocatable              :: save0(:)

eps = 1.0d-12
k   = 30

call chebexps(k,chebdata)
allocate(xsings(0))


pi     = acos(-1.0d0)

!
!  Measure the performance of the solver for fixed wavenumber and increasing n
!

dlambda = 2**17
delta   = 2**8

nn      = dlambda/delta
R       = 2

allocate(dtimes(1:nn))
allocate(derrs(1:nn))
allocate(dnus(1:nn))
allocate(dlog(1:nn))
allocate(soldata(1:nn))

!

ncount = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,lambda,dnu,t1,t2,dtime,derr,idx,t,val,val0,cc, &
!$OMP  dmax,a0,b0,i,mm)
!$OMP DO
do idx=1,nn

n      = delta*idx
lambda = dlambda
dnu    = n

call elapsed(t1)
call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata(idx))
call elapsed(t2)
dtime  = t2-t1

! normalize the solution
t     = 1
call pbessel_eval(soldata(idx),t,val,der)
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)
cc    = val0 / val
dmax  = 0

 ! measure the error at 10 points on the interval
dmax  = 0.0d0
mm    = 100
a0    = 1.0d0
b0    = 2.0d0

do i=1,mm
t     = a0 + (b0-a0)/(mm-1.0d0) * (i-1.0d0)
call pbessel_eval(soldata(idx),t,val,der)
val   = cc*val
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)

derr   = abs(val-val0)


dmax  = max(dmax,derr)

end do

write (*,"(I10,3X,I13,3X,D10.3,3X,D10.3,A)",advance='no')   lambda,n,dtime,dmax,13

dtimes(idx) = dtime
derrs(idx)  = dmax
dnus(idx)   = dnu
!dlog(idx)   = log(dnu)/log(dnus(1)) * dtimes(1)

end do
!$OMP END DO
!$OMP END PARALLEL



nn0 = 2
do idx = 1,nn
dnu           = dnus(idx)
dlog(idx)     = log(dnu) * dtimes(nn0)/log(dnus(nn0))
!dlog(idx)     = log(dnu)*c
end do


!
!  Produce the timing graph
!


iw = 101
open(iw,FILE='pbessel_graph1.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"Runtime",nn,dnus,dtimes)
call pyplot_add_function(istatus,10,"k log(k)",nn,dnus,dlog)

call pyplot_ylabel(istatus,"Solve time (seconds)")
call pyplot_xlabel(istatus,"Order n")
call pyplot_xlogscale2(istatus)


! write(iw,"(A)") "ax.text(0.90, 0.90, 'log(n)',"
! write(iw,"(A)") "verticalalignment='bottom', horizontalalignment='right',"
! write(iw,"(A)") "transform=ax.transAxes,"
! write(iw,"(A)") "color='black', fontsize=12)"

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph1.pdf")
call system("python pbessel_graph1.py")

!
!  Produce the error graph
!

iw = 101
open(iw,FILE='pbessel_graph2.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"",nn,dnus,derrs)
! call pyplot_add_function(istatus,10,"",nn,dnus,dlog)

call pyplot_ylabel(istatus,"Maximum absolute error")
call pyplot_xlabel(istatus,"Order n")
call pyplot_xlogscale2(istatus)
call pyplot_ylogscale(istatus)

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph2.pdf")
call system("python pbessel_graph2.py")


deallocate(dtimes,derrs,dnus,dlog,soldata)

!
!  Measure the performance of the solver for n=0 and increasing dlambda
!


dlambda = 2**17
delta   = 2**8

nn      = dlambda/delta
R       = 2


allocate(soldata(1:nn))
allocate(dtimes(1:nn))
allocate(derrs(1:nn))
allocate(dlog(1:nn))
allocate(dlambdas(1:nn))

do idx=1,nn

dlambda = delta*idx
n       = 0
dnu     = n 

call elapsed(t1)
call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata(idx))
call elapsed(t2)
dtime  = t2-t1

! normalize the solution

t     = 1
call pbessel_eval(soldata(idx),t,val,der)
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)
cc    = val0 / val



!  ! measure the error at 100  points on the interval
dmax  = 0.0d0
mm    = 100
a0    = 1.0d0
b0    = 2.0d0

do i=1,mm
t     = a0 + (b0-a0)/(mm-1.0d0) * (i-1.0d0)
call pbessel_eval(soldata(idx),t,val,der)
val   = cc*val
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)

derr   = abs(val-val0)
dmax  = max(dmax,derr)
end do

dtimes(idx)   = dtime
derrs(idx)    = dmax
dlambdas(idx) = dlambda

lambda=dlambda
write (*,"(I10,3X,I13,3X,D10.3,3X,D10.3,A)",advance='no')   lambda,n,dtime,dmax,13

end do

nn0 = 1
do idx = 1,nn
dlambda       = dlambdas(idx)
dlog(idx)     = log(dlambda)/log(dlambdas(nn0)) * dtimes(nn0)
end do


!
!  Produce the timing graph
!


iw = 101
open(iw,FILE='pbessel_graph3.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"Runtime",nn,dlambdas,dtimes)
call pyplot_add_function(istatus,10,"k log(k)",nn,dlambdas,dlog)

call pyplot_ylabel(istatus,"Solve time (seconds)")
call pyplot_xlabel(istatus,"Wavenumber k")
call pyplot_xlogscale2(istatus)

! write(iw,"(A)") "ax.text(0.90, 0.90, 'log(k)',"
! write(iw,"(A)") "verticalalignment='bottom', horizontalalignment='right',"
! write(iw,"(A)") "transform=ax.transAxes,"
! write(iw,"(A)") "color='black', fontsize=12)"

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph3.pdf")
call system("python pbessel_graph3.py")

!
!  Produce the error graph
!

iw = 101
open(iw,FILE='pbessel_graph4.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"",nn,dlambdas,derrs)

call pyplot_ylabel(istatus,"Maximum absolute error")
call pyplot_xlabel(istatus,"Wavenumber k")
call pyplot_xlogscale2(istatus)
call pyplot_ylogscale(istatus)

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph4.pdf")
call system("python pbessel_graph4.py")


deallocate(dtimes,derrs,dlambdas,dlog,soldata)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Measure the performance of the solver for n = dlambda/2 incerasing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


dlambda = 2**17
delta   = 2**8

nn      = dlambda/delta
R       = 2


allocate(soldata(1:nn))
allocate(dtimes(1:nn))
allocate(derrs(1:nn))
allocate(dlog(1:nn))
allocate(dlambdas(1:nn))

do idx=1,nn

dlambda = delta*idx
n       = dlambda/2
dnu     = n 

call elapsed(t1)
call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata(idx))
call elapsed(t2)
dtime  = t2-t1

! normalize the solution

t     = 1
call pbessel_eval(soldata(idx),t,val,der)
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)
cc    = val0 / val

 ! measure the error at 100 points on the interval
dmax  = 0.0d0
mm    = 100
a0    = 1.0d0
b0    = 2.0d0

do i=1,mm
t     = a0 + (b0-a0)/(mm-1.0d0) * (i-1.0d0)
call pbessel_eval(soldata(idx),t,val,der)
val   = cc*val
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)

derr   = abs(val-val0)
dmax  = max(dmax,derr)
end do

dtimes(idx)   = dtime
derrs(idx)    = dmax
dlambdas(idx) = dlambda

lambda = dlambda
write (*,"(I10,3X,I13,3X,D10.3,3X,D10.3,A)",advance='no')   lambda,n,dtime,dmax,13

end do


nn0 = 1
do idx = 1,nn
dlambda       = dlambdas(idx)
dlog(idx)     = log(dlambda)/log(dlambdas(nn0)) * dtimes(nn0)
end do


!
!  Produce the timing graph
!


iw = 101
open(iw,FILE='pbessel_graph5.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"Runtime",nn,dlambdas,dtimes)
call pyplot_add_function(istatus,10,"k log(k)",nn,dlambdas,dlog)

call pyplot_ylabel(istatus,"Solve time (seconds)")
call pyplot_xlabel(istatus,"n=k/2")
call pyplot_xlogscale2(istatus)


! write(iw,"(A)") "ax.text(0.90, 0.90, 'log(n)',"
! write(iw,"(A)") "verticalalignment='bottom', horizontalalignment='right',"
! write(iw,"(A)") "transform=ax.transAxes,"
! write(iw,"(A)") "color='black', fontsize=12)"

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph5.pdf")
call system("python pbessel_graph5.py")

!
!  Produce the error graph
!

iw = 101
open(iw,FILE='pbessel_graph6.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"",nn,dlambdas,derrs)

call pyplot_ylabel(istatus,"Maximum absolute error")
call pyplot_xlabel(istatus,"n=k/2")
call pyplot_xlogscale2(istatus)
call pyplot_ylogscale(istatus)

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph6.pdf")
call system("python pbessel_graph6.py")

deallocate(dtimes,derrs,dlambdas,dlog,soldata)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Measure the performance of the solver for n = dlambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


dlambda = 2**17
delta   = 2**8

nn      = dlambda/delta
R       = 2


allocate(soldata(1:nn))
allocate(dtimes(1:nn))
allocate(derrs(1:nn))
allocate(dlog(1:nn))
allocate(dlambdas(1:nn))

do idx=1,nn

dlambda = delta*idx
n       = dlambda
dnu     = n 

call elapsed(t1)
call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata(idx))
call elapsed(t2)
dtime  = t2-t1

! normalize the solution

t     = 1
call pbessel_eval(soldata(idx),t,val,der)
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)
cc    = val0 / val

 ! measure the error at 100 points on the interval
dmax  = 0.0d0
mm    = 100
a0    = 1.0d0
b0    = 2.0d0

do i=1,mm
t     = a0 + (b0-a0)/(mm-1.0d0) * (i-1.0d0)
call pbessel_eval(soldata(idx),t,val,der)
val   = cc*val
val0  = bessel_jn(n/2,dlambda/2*t**2)
val0  = val0 * sqrt(t)

derr   = abs(val-val0)
dmax  = max(dmax,derr)
end do

dtimes(idx)   = dtime
derrs(idx)    = dmax
dlambdas(idx) = dlambda

lambda=dlambda
write (*,"(I10,3X,I13,3X,D10.3,3X,D10.3,A)",advance='no')   lambda,n,dtime,dmax,13

end do

nn0 = 1
do idx = 1,nn
dlambda       = dlambdas(idx)
dlog(idx)     = log(dlambda)/log(dlambdas(nn0)) * dtimes(nn0)
end do


!
!  Produce the timing graph
!


iw = 101
open(iw,FILE='pbessel_graph7.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"Runtime",nn,dlambdas,dtimes)
call pyplot_add_function(istatus,10,"k log(k)",nn,dlambdas,dlog)

call pyplot_ylabel(istatus,"Solve time (seconds)")
call pyplot_xlabel(istatus,"n=k")
call pyplot_xlogscale2(istatus)


! write(iw,"(A)") "ax.text(0.90, 0.90, 'log(n)',"
! write(iw,"(A)") "verticalalignment='bottom', horizontalalignment='right',"
! write(iw,"(A)") "transform=ax.transAxes,"
! write(iw,"(A)") "color='black', fontsize=12)"

idelta=  delta
lambda= dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph7.pdf")
call system("python pbessel_graph7.py")

!
!  Produce the error graph
!

iw = 101
open(iw,FILE='pbessel_graph8.py')
call pyplot_begin(iw,istatus)

istyle = 1
call pyplot_add_function(istatus,1,"",nn,dlambdas,derrs)

call pyplot_ylabel(istatus,"Maximum absolute error")
call pyplot_xlabel(istatus,"k=n")
call pyplot_xlogscale2(istatus)
call pyplot_ylogscale(istatus)

idelta =  delta
lambda = dlambda
write(iw,"(A,I12,A,I12,A)") "plt.xlim(",idelta,",",lambda,")"

call pyplot_end(istatus,"pbessel_graph8.pdf")
call system("python pbessel_graph8.py")


deallocate(dtimes,derrs,dlambdas,dlog,soldata)




end program
