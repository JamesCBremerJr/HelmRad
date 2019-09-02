program test_helmop
use utils
use chebyshev
use helmop

implicit double precision (a-h,o-z)

double precision, allocatable        :: ts(:), rs(:)
double complex, allocatable          :: vals(:), ders(:),ders0(:)
double complex, allocatable          :: x(:,:),y(:,:),y0(:,:)
double complex                       :: ima

pi = acos(-1.0d0)

ima = (0.0d0,1.0d0)
n   = 120
m   = 240

R1  = 1.0d-1
R2  = 2

call helmop_init(R1,R2,n,m,ts,rs)
call prind("R1 = ",R1)
call prind("R2 = ",R2)
call prini("n = ",n)
call prini("m = ",m)
nn = (2*n+1)*m
call prini("total # discretization nodes = ",nn)
call prina("")

!
!  Test the routine for applying (d/dt)^2
!

allocate(vals(2*n+1),ders(2*n+1),ders0(2*n+1))
vals  = exp(sin(ts))
ders0 = exp(sin(ts))*cos(ts)**2-exp(sin(ts))*sin(ts)


vals  = 1
ders0 = 0

call elapsed(t1)
call helmop_apply_ddt(R1,R2,n,ts,vals,ders)
call elapsed(t2)

derr = maxval(abs(ders-ders0))

call prin2("(d/dt)^2 error = ",derr)
call prin2("apply_ddt time = ",t2-t1)
call prina("")

deallocate(vals,ders,ders0)


!
!  Test the routine for applying d/dr
!

allocate(vals(m+1),ders(m+1),ders0(m+1))

vals  = rs**2 + ima * rs
vals  = ima*rs
ders0 = ima

vals  = rs**2  + ima * rs**3
ders0 = 2*rs   + ima * 3*rs**2

call elapsed(t1)
call helmop_apply_dr(R1,R2,m,rs,vals,ders)
call elapsed(t2)

derr = maxval(abs(ders-ders0))

call prin2("d/dr error = ",derr)
call prin2("apply_dr time = ",t2-t1)
call prina("")


deallocate(vals,ders,ders0)

!
!  Test the application of the Laplacian
!

allocate(x(2*n+1,m+1),y(2*n+1,m+1),y0(2*n+1,m+1))

dlambda = 3
k       = 5

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,t,r)
!$OMP DO
do i=1,2*n+1
do j=1,m+1
t = ts(i)
r = rs(j)

x(i,j)   = bessel_jn(k,dlambda*r)*exp(ima*k*t)
y0(i,j)  = -dlambda**2*bessel_jn(k,dlambda*r)*exp(ima*k*t)

end do
end do
!$OMP END DO
!$OMP END PARALLEL

call elapsed(t1)
call helmop_delta(R1,R2,n,m,ts,rs,x,y)
call elapsed(t2)


derr = maxval( abs(y - y0) )

call prin2("delta error = ",derr)
call prin2("delta time  = ",t2-t1)

deallocate(x,y0,y)



end program

