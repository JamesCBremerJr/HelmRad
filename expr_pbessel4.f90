module test_pbessel2_functions

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

if  (ising .eq. 1) qs = 0
if  (ising .eq. 2) qs = 3


end subroutine


end module


program test_pbessel2
use utils
use chebyshev
use pbessel
use iso_c_binding
use test_pbessel2_functions
implicit double precision (a-h,o-z)

type(chebexps_data)                        :: chebdata
type(c_ptr)                                :: userptr
type(pbessel_solution), allocatable        :: soldata(:)
type(pbessel_solution)                     :: soldata0
double precision, allocatable              :: xsings(:)
double precision, allocatable              :: xs(:),ys(:)
double precision, allocatable              :: dlambdas(:),dlog(:),derrs(:),dtimes(:),xs0(:),ys0(:)

double precision, allocatable              :: save0(:)

eps = 1.0d-12
k   = 30

call chebexps(k,chebdata)
allocate(xsings(1))
xsings(1) = 1


pi     = acos(-1.0d0)

write (*,"(A10,3X,A10,3X,A10,3XA10)")       "wavenumber","time","ratio","max error"
write (*,*)                                 ""

write (13,"(A10,3X,A10,3X,A10,3XA10)")       "wavenumber","time","ratio","max error"
write (13,*)                                 ""




dtime  = 0
jj1    = 8
jj2    = 17

allocate(dtimes(jj1:jj2),derrs(jj1:jj2))
allocate(dlog(jj1:jj2),dlambdas(jj1:jj2))

do jj=jj1,jj2

eps     = 1.0d-12 
R       = 2.0d0
dtime2  = dtime
dlambda = 2**jj
m1      = 0
m2      = dlambda
dmax    = 0
ifadap  = 0

allocate(soldata(m1:m2))

!
!  Solve the perturbed Bessel equation for n=0,...,dlambda
!

call elapsed(t1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,dnu)
!$OMP DO
do n = m1,m2
dnu = n
call pbessel_solve(eps,chebdata,xsings,R,dlambda,dnu,qfun,userptr,soldata(n))
end do
!$OMP END DO
!$OMP END PARALLEL
call elapsed(t2)

dtime = t2-t1

!
!
!  Test the accuracy of the solutions
!


dmax = 0
do n=1,m2
dnu   = n 
dmax0 = 0

! normalize the solution
t     = 1.0d0
call pbessel_eval(soldata(n),t,val,der)
val0 = (dlambda*Pi*sqrt(t)*(bessel_jn(n,0.2d1*dlambda*t)*(0.2d1*bessel_jn(n,dlambda)*bessel_yn(-1  &
  +n,2*dlambda)-0.1d1*bessel_jn(-1+n,dlambda)*bessel_yn(n,2*dlambda))+  &
  (-0.2d1*bessel_jn(-1+n,2*dlambda)*bessel_jn(n,dlambda)+bessel_jn(-1+  &
  n,dlambda)*bessel_jn(n,0.2d1*dlambda))*bessel_yn(n,0.2d1*dlambda*t)))/0.2d1

cc = val0/val


!
! measure the error at 100 points
dmax0 = 0
nn    = 100
a0    = 1.0d0
b0    = 2.0d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,t,val,der,val0,derr)
!$OMP DO
do i=1,nn

t     = a0 + (b0-a0)/(nn-1.0d0) * (i-1.0d0)

call pbessel_eval(soldata(n),t,val,der)
val   = cc*val

val0 = (dlambda*Pi*sqrt(t)*(bessel_jn(n,0.2d1*dlambda*t)*(0.2d1*bessel_jn(n,dlambda)*bessel_yn(-1  &
  +n,2*dlambda)-0.1d1*bessel_jn(-1+n,dlambda)*bessel_yn(n,2*dlambda))+  &
  (-0.2d1*bessel_jn(-1+n,2*dlambda)*bessel_jn(n,dlambda)+bessel_jn(-1+  &
  n,dlambda)*bessel_jn(n,0.2d1*dlambda))*bessel_yn(n,0.2d1*dlambda*t)))/0.2d1


derr   = abs(val-val0) 
dmax0  = max(derr,dmax0)

end do
!$OMP END DO
!$OMP END PARALLEL


dmax         = max(dmax0,dmax)

end do



dlambdas(jj) = dlambda
dtimes(jj)   = dtime
dlog(jj)     = log(dlambda)*dlambda * dtimes(jj1) / ( log(dlambdas(jj1)) * dlambdas(jj1))
derrs(jj)    = dmax

lambda = dlambda
if (dtime2 .ne. 0) then
dratio= dtime/dtime2
else
dratio = -1
endif

write (*,"(I10,3X,D10.3,3X,D10.3,3X,D10.3,3X,D10.3)")   lambda,dtime,dratio,dmax
write (13,"(I10,3X,D10.3,3X,D10.3,3X,D10.3,3X,D10.3)")   lambda,dtime,dratio,dmax


dtime2 = dtime

deallocate(soldata)

end do


!
!  Produce LaTeX source code for a table
!

iw = 25
open(iw,FILE='pbessel_table3.tex')


write(iw,'(A)') "\begin{tabular}{cccc}"
write(iw,'(A)') "\toprule" 
write (iw,*) "\addlinespace[.25em]"
write (iw,*) "\addlinespace[.125em]"

write(iw,'(A)', advance='no') "$k$                                &"
write(iw,'(A)', advance='no') "Total time                         &"
write(iw,'(A)', advance='no') "Ratio of                          &" 
write(iw,'(A)'              ) "Maximum absolute                          \\"


write(iw,'(A)', advance='no') "                                   &"
write(iw,'(A)', advance='no') "(in seconds)                       &"
write(iw,'(A)', advance='no') "times                              &"
write(iw,'(A)'              ) "error                               \\"


write(iw,'(A)')               "\midrule"
write (iw,*) "\addlinespace[.25em]"


do jj=jj1,jj2


if (jj .ge. 10) then
write(iw,"('$2^{',I2,'}$ ')",advance='no') jj
else
write(iw,"('$2^{',I1,'}$')",advance='no') jj
endif

call write_table_next(iw)


call write_table_double(iw,dtimes(jj))
call write_table_next(iw)

if (jj .eq. jj1) then
write (iw,*) "-"
else
call write_table_double(iw,dtimes(jj)/dtimes(jj-1))
endif

call write_table_next(iw)

call write_table_double(iw,derrs(jj))

call write_table_nextline(iw)

write (iw,*) "\addlinespace[.125em]"

end do


write (iw,'(A)') "\bottomrule"
write (iw,'(A)') "\end{tabular}"
close(iw)

!
!  Produce the timing graph
!

nn = jj2-jj1+1



iw = 101
open(iw,FILE='pbessel_graph13.py')

call pyplot_begin(iw,istatus)


call pyplot_add_function(istatus,1,"Runtime",nn,dlambdas,dtimes)
call pyplot_add_function(istatus,10,"k log(k)",nn,dlambdas,dlog)

! write(iw,"(A)") "ax.text(0.85, 0.90, 'k log(k)',"
! write(iw,"(A)") "verticalalignment='bottom', horizontalalignment='right',"
! write(iw,"(A)") "transform=ax.transAxes,"
! write(iw,"(A)") "color='black', fontsize=12)"

call pyplot_ylabel(istatus,"Solve time (seconds)")
call pyplot_xlabel(istatus,"Wavenumber k")
call pyplot_xlogscale2(istatus)

call pyplot_end(istatus,"pbessel_graph13.pdf")
call system("python pbessel_graph13.py")

!
!  Plot the function q(t)
!

nn = 1000
allocate(xs0(nn),ys0(nn))

a0 = 0
b0 = R

do i=1,nn
xs0(i) = a0 + (i-1.0d0) /(nn-1.0d0) * (b0-a0)
end do


call qfun(1,nn/2,xs0,ys0,userptr)
call qfun(2,nn/2,xs0(nn/2+1:nn),ys0(nn/2+1:nn),userptr)


iw = 101
open(iw,FILE='pbessel_graph14.py')
call pyplot_begin(iw,istatus)

call pyplot_add_function(istatus,1,"",nn,xs0,ys0)
call pyplot_ylabel(istatus,"q(t)")
call pyplot_xlabel(istatus,"t")

write(iw,"(A,I12,A,I12,A)") "plt.ylim(-1,4)"

call pyplot_end(istatus,"pbessel_graph14.pdf")
call system("python pbessel_graph14.py")


end program
