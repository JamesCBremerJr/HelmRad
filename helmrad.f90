!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains a solver for the boundary value problem
!
!   {
!   {   \Delta u_s(r,t) + dlambda^2 (1 + q(r) ) u_s(r,t) 
!   {                                                              for all 0 < r < R
!   {                              =  - lambda^2 q(r) u_i(t)               0 < t < 2pi     (1)
!   {
!   {                          ( d u_s(r,t)                      )                         
!   {   lim  \sqrt{r}   sup    ( ----------  - i lambda u_s(r,t) ) = 0,
!   {    r             0<t<2pi (     dr                          )
!   {
!
!  where \Delta denotes the Laplacian, q(r) is positive and piecewise smooth with 
!  compact support in (0,R), dlambda is real-valued and u_i is an incident wave 
!  which satisfies the constant coefficient Helmholtz equation
!
!   \Delta u_i(r,t)  + dlambda^2 u_i(r,t) = 0    for all 0 < r < \infty and 0 < t < 2 \pi.
!
!  In addition to proving an external subroutine for evaluating q, the user must
!  supply a list of the points at which q is not smooth.  This  list of points
!  partitions (0,R) into subintervals and whenever the user-supplied external subroutine 
!  defining q is called, the index of the interval containing all of these arguments
!  is provided.
!
!  The run time of this solver is O( dlambda log(dlambda) ). It operates via separation of 
!  variables; that is, the total solution u, which is the sum of  u_s + u_i, is represented 
!  in the form
!
!       m                            
!     \sum  a_n u_n(r ) \exp(i n t)                                                       (2)
!      n=-m                        
!
!  where, for each n, u_n(r) is a  solution of the boundary value problem
!
!    {
!    {  z^2 * u_n''(r) + z u_n(r) + ( dlambda^21 (1 + q(r)) r^2 - n^2 ) u_n(r) = 0
!    {                                                                                    (3)
!    {  u_n(t) = O (t^dnu) as t --> 0                                                     
!    {
!
!  and m is chosen sufficiently large to represent the solution in (2).  It suffices
!  to choose m sufficiently large so that the restriction of the incident wave to the
!  unit circle can be represented via a Fourier series of the form
!
!      m
!    \sum c_n exp(int).
!     n=-m
!
!  The routine helmrad_init performs the initialization phase for the solver
!  in which the solutions of the boundary value problem (3) for n = 0,...,m are
!  are computed.  The raison d'etre for this code is that this procedure is carried
!  out in time which grows at worst logarithmically with k.
!  
!  We use Paul N. Swarztrauber's fftpack library to perform the fast Fourier 
!  transform.
!
!  The following subroutines are publicly callable:
!
!    helmrad_init - initialize the solver for a given scattering potential and
!      wavenumber
!
!    helmrad_solve - given an external subrouting specifying the incident wave
!      u_i, construct a solution of the boundary value problem (1)
!
!    helmrad_eval - evaluate the total field and the scattered field at a specified 
!       point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module helmrad

use utils
use chebyshev
use pbessel
use iso_c_binding

! The helmrad_data structure stores the data necessary to solve the boundary
! value problems (*) and (**) for a given potential q and wavenumber lambda
type       helmrad_data
integer                                    :: m
integer                                    :: k
double precision                           :: dlambda
double precision                           :: R
type(pbessel_solution), allocatable        :: sols(:)

! data for the fft
 double precision, allocatable              :: wsave(:), ts(:)

end type   helmrad_data


! The helmrad_scat structure describes a solution of the BVP (1)
type       helmrad_scat

integer                                   :: m
double precision                          :: dlambda

! Fourier expansions of the incident wave and its normal derivative
double complex, allocatable               :: in_coefs(:),inder_coefs(:)

! coefficients in the H expansion of the scattered wave
double complex, allocatable               :: scat_coefs(:)

! coefficients in the expansion of the total field 
double complex, allocatable               :: tot_coefs(:)

end type   helmrad_scat

interface


subroutine helmrad_potfun(ising,n,rs,vals,userptr)
use iso_c_binding
double precision :: rs(n)
double precision :: vals(n)
!
!  This is the interface for the external routines supplying the values
!  of the radially symmetric scattering potential.
!
!  Input parameters:
!    ising - an integer parameter indicating which of the intervals formed
!       by the singular points of q the input arguments are contained in 
!    n - the number of points at which to evaluate the incident wave and
!      its derivative
!    rs - an array containing the list of points r at which to evaluate
!      the scattering potential
!   userptr - a "void *" pointer which provided by the user to the helmrad_init
!      routine and passed on to this routine
!
!  Output parameters:
!    vals - the values of the scattering potential q(r) at the points
!      specified in the rs array
!
type(c_ptr)      :: userptr
end subroutine


subroutine helmrad_wavefun(dlambda,n,rs,ts,vals,ders,userptr)
use iso_c_binding
double precision :: dlambda
double precision :: rs(n)
double precision :: ts(n)
double complex   :: vals(n)        
double complex   :: ders(n)
type(c_ptr)      :: userptr

!
!  This is the interface for the external routines which supply the values
!  of the incident wave u_i(r,t) on the unit circle and the normal derivative
!  of u_i.
!
!  Input parameters:
!    dlambda - the wavenumber for the problem
!    n - the number of points at which to evaluate the incident wave and
!      its derivative
!    ts - an array containg the points on the unit circle (expressed as
!      angles) at which to evaluate u_i and its derivative with respect
!      to r
!   userptr - a "void *" pointer which provided by the user to the 
!      helmrad_scat_solve routine and passed along to this subroutine
!
!  Output parameters:
!    vals - the values of the incident wave u_i(1,t) at the points t
!    ders - the values of the the d/dr u_i(1,t) 
!
end subroutine

end interface

contains


subroutine helmrad_init(ifout,eps,m,xs0,R,dlambda,qfun,userptr,helmdata)
implicit double precision (a-h,o-z)
integer                         :: m
procedure(helmrad_potfun)       :: qfun
type(c_ptr)                     :: userptr
type(helmrad_data), intent(out) :: helmdata
double precision, allocatable   :: xs0(:)
!
!  Perform the initialization stage by solving the sequence of boundary value 
!  problems of the form (3) with n = 0,1,...,m.  The scattering potential
!  q(t) is specified via an external subroutine.
!
!  Input parameters:
!    eps - the precision for the solutions of (3)
!    m - a parameter controlling the number of terms in the representation (2)
!      of the solution of (1)
!    xs0 - a list of the singular points of q
!    dlambda - the wavenumber of the problem
!    R - the radius of a disk containing the support of the scattering potential q(r)
!    qfun - a user-supplied external subroutine conforming to the pbessel_qfun 
!      interface which supplies the values of the scattering potential q(t)
!    userptr - a "void *" pointer which is passed to the user-supplied external
!      subroutine qfun
!
!  Output parameters:
!    helmdata - a data structure produced by the solver which stores the
!      information necessary to solve (1) 
!

data pi / 3.14159265358979323846264338327950288d0 /

type(chebexps_data)             :: chebdata

k       = 30
call chebexps(k,chebdata)


helmdata%m       = m
helmdata%k       = k
helmdata%dlambda = dlambda
helmdata%R       = R

allocate(helmdata%sols(0:m))
allocate(helmdata%ts(2*m+1))

dd = 2*pi/(2*m+1.0d0)
do i=1,2*m+1
helmdata%ts(i) = dd* (i-1)
end do

ncount = 0
call elapsed(t1)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,dnu,ier,tt1,tt2)
!$OMP DO SCHEDULE(DYNAMIC)
do n=0,m

dnu = n 
call  pbessel_solve(eps,chebdata,xs0,R,dlambda,dnu,qfun,userptr,helmdata%sols(n))

if (ifout .eq. 1) then
!$OMP CRITICAL
ncount = ncount + 1
call elapsed(t2)
write (*,"(A,I9,' / ',I9,A,F13.4,A)",advance='no'), &
  'helmrad_init: ',ncount,m+1,',  elapsed time = ',t2-t1,' seconds'
write (*,"(A)",advance='no'),13
!$OMP END CRITICAL
endif

end do
!$OMP END DO
!$OMP END PARALLEL

if (ifout .eq. 1) then
write (*,"(A)") ""

write (13,"(A,I9,' / ',I9,A,F13.4,A)"), &
  'helmrad_init: ',ncount,m+1,',  elapsed time = ',t2-t1,' seconds'
endif



end subroutine


subroutine helmrad_solve(helmdata,wavefun,userptr,scatdata)
implicit double precision (a-h,o-z)
type(helmrad_data)                         :: helmdata
type(helmrad_scat), intent(out)            :: scatdata
procedure(helmrad_wavefun)                 :: wavefun
type(c_ptr)                                :: userptr
!
!  Compute the coefficients in the expansion (2) of the solution of (1)
!  given the data generated by helmrad_init and a user-supplied subroutine
!  which specifies the values of the incident wave and its derivative.
!
!  Input parameters:
!    helmdata - the data structure returned by helmrad_init
!    wavefun - an external subroutine conforming to the helmrad_wavefun
!       interface which supplies the values of the incident wave and its
!       normal derivative
!    userptr - a "void *" pointer which is passed to wavefun
! 
!  Output parameters:
!    scatdata - a data structure which stores information about the
!     solution
!

double precision, allocatable :: ts(:), rs(:)
double complex, allocatable   :: ders(:), wsave(:), vals(:), hanks(:)
double complex                :: ima
double complex                :: amatr(2,2)
double complex                :: za,zb,zc,zd,det,dd1,dd2,hval,hval1,hval2,hder,val0,z


data pi  / 3.14159265358979323846264338327950288d0 /
data ima / (0.0d0,1.0d0) /


m                = helmdata%m
dlambda          = helmdata%dlambda
R                = helmdata%R
scatdata%m       = m


allocate(scatdata%in_coefs(-m:m))
allocate(scatdata%inder_coefs(-m:m))
allocate(scatdata%tot_coefs(-m:m))
allocate(scatdata%scat_coefs(-m:m))

! compute the Fourier coefficient of the restriction of the incident wave to the unit circle
! and of its normal derivative

allocate(rs(2*m+1))
allocate(ders(2*m+1))
allocate(vals(2*m+1))


scatdata%in_coefs = 0
scatdata%inder_coefs = 0
scatdata%tot_coefs = 0
scatdata%scat_coefs = 0

rs = R
call wavefun(dlambda,2*m+1,rs,helmdata%ts,vals,ders,userptr)

allocate(helmdata%wsave(10*m+1000))
call zffti(2*m+1,helmdata%wsave)
call zfftf(2*m+1,vals,helmdata%wsave)
call zfftf(2*m+1,ders,helmdata%wsave)

vals    = vals*1.0d0/(2*m+1.0d0)
ders    = ders*1.0d0/(2*m+1.0d0)

do i=1,m+1
scatdata%in_coefs(i-1)    = vals(i)
scatdata%inder_coefs(i-1) = ders(i)
end do

do i=1,m
scatdata%in_coefs(-i)    = vals(2*m+2-i)
scatdata%inder_coefs(-i) = ders(2*m+2-i)
end do

! compute the coefficients in the expansions of the scattering wave and the
! total field

val0     = 0
val_scat = 0

z        = dlambda*R
ifexpon  = 1
allocate(hanks(0:m+1))
call hanks103(z,hanks,m+1,ifexpon)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,nn,ifosc,psi,psip,s1,aval,apval,appval, &
!$OMP   uval,uder,amatr,det,zc,zd,za,zb,hval,hval1,hval2,hder,uval0,uder0)
 !$OMP DO SCHEDULE(dynamic)
do n=0,m

nn = abs(n)

call pbessel_eval(helmdata%sols(nn),R,uval,uder)
uval = uval/sqrt(R)
uder = uder / sqrt(R) - uval/(2*R)

hval  = hanks(nn)
hval1 = hanks(nn+1)

if (nn .gt. 0) then
hval2 = hanks(nn-1)
else
hval2 = -hanks(1)
!hval2 = bessel_jn(nn-1,dlambda*R) + ima * bessel_yn(nn-1,dlambda*R)
endif

! hval  = bessel_jn(nn,dlambda*R)   + ima * bessel_yn(nn,dlambda*R)
! hval1 = bessel_jn(nn+1,dlambda*R) + ima * bessel_yn(nn+1,dlambda*R)
! hval2 = bessel_jn(nn-1,dlambda*R) + ima * bessel_yn(nn-1,dlambda*R)

hder  = dlambda/2 * (hval2-hval1)

!
!  Compute a_n and b_n for nonnegative n
!

amatr(1,1) = uval
amatr(1,2) = -hval
amatr(2,1) = uder
amatr(2,2) = -hder

det = amatr(1,1)*amatr(2,2) - amatr(1,2) * amatr(2,1)


zc = scatdata%in_coefs(n)
zd = scatdata%inder_coefs(n)

za = 1.0d0/det  * (amatr(2,2)*zc  - amatr(1,2)*zd)
zb = 1.0d0/det  * (-amatr(2,1)*zc + amatr(1,1)*zd)

scatdata%tot_coefs(n)  = za
scatdata%scat_coefs(n) = zb

if (n .eq. 0) continue

!
!  Compute a_n and b_n for negative n
!

amatr(1,1) = uval
amatr(1,2) = -hval
amatr(2,1) = uder
amatr(2,2) = -hder

det = amatr(1,1)*amatr(2,2) - amatr(1,2) * amatr(2,1)

zc = scatdata%in_coefs(-n)
zd = scatdata%inder_coefs(-n)

za = 1.0d0/det  * (amatr(2,2)*zc  - amatr(1,2)*zd)
zb = 1.0d0/det  * (-amatr(2,1)*zc + amatr(1,1)*zd)

scatdata%tot_coefs(-n)  = za
scatdata%scat_coefs(-n) = zb

end do
!$OMP END DO
!$OMP END PARALLEL


end subroutine



subroutine helmrad_eval(helmdata,scatdata,wavefun,userptr,r,t,val_tot,val_scat)
implicit double precision (a-h,o-z)
type(helmrad_data)                         :: helmdata
type(helmrad_scat)                         :: scatdata
procedure(helmrad_wavefun)                 :: wavefun
type(c_ptr)                                :: userptr
double complex                             :: val_tot, val_scat
!
!  Evaluate the scattered field and the total field at a point in the plane specified 
!  via its polar coordinates.
!
!  Input parameters:
!    helmdata - the data structure returned by helmrad_init routine
!    scatdata - the data structure returned by the helmrad_scat_solve routine
! 
!  Output parameters:
!    val_tot  - the value of the total wave at the point x  r exp(it)
!    val_scat - the value of the scattered wave at the point x = r exp(it)
!

double precision            :: ts(1),rs(1)
double complex              :: vals(1),ders(1),hval,ima,z
double complex, allocatable :: hanks(:)


data pi  / 3.14159265358979323846264338327950288d0 /
data ima / (0.0d0,1.0d0) /

ima     = (0.0d0,1.0d0)
m       = helmdata%m
RR      = helmdata%R
dlambda = helmdata%dlambda

!
!  Handle the case of r < R
!

if (r .le. RR) then

val_tot = 0

do n=0,m
nn = abs(n)
call pbessel_eval(helmdata%sols(nn),r,uval,uder)
uval = uval/sqrt(r)

val_tot  = val_tot + uval * exp(ima*n*t)*scatdata%tot_coefs(n)

if (n .gt. 0) then
val_tot  = val_tot + uval * exp(-ima*n*t)*scatdata%tot_coefs(-n)
endif

end do

rs(1) = r
ts(1) = t
call wavefun(dlambda,1,rs,ts,vals,ders,userptr)
val_scat = val_tot - vals(1)

return
endif



!
!  Handle r > RR
!


val_scat = 0
z        = dlambda*r
ifexpon  = 1

allocate(hanks(0:m))
call hanks103(z,hanks,m,ifexpon)

do n=0,m
nn   = abs(n)

val_scat  = val_scat + hanks(n) * exp(ima*n*t)*scatdata%scat_coefs(n)

if (n .gt. 0) then
val_scat  = val_scat + hanks(n) * exp(-ima*n*t)*scatdata%scat_coefs(-n)
endif

end do


rs(1) = r
ts(1) = t
call wavefun(dlambda,1,rs,ts,vals,ders,userptr)
val_tot = val_scat + vals(1)

end subroutine



subroutine helmrad_inwave(eps,dlambda,R,wavefun,userptr,m)
implicit double precision (a-h,o-z)
procedure(helmrad_wavefun)                 :: wavefun
type(c_ptr)                                :: userptr
!
!  Adaptively determine a value of m such that the sum (4) represents
!  the incoming wave to high accuracy.
!
!  Input parameters:
!    eps - the precision for the representation
!    dlambda - the wavenumber for the problem
!    R - the radius of a disk containing the support of the scattering potential q(r)
!    wavefun - an external subroutine conforming to the helmrad_wavefun
!       interface which supplies the values of the incident wave and its
!       normal derivative
!    userptr - a "void *" pointer which is passed to wavefun
! 
!  Output parameters:
!    m - an appropriate value for the integer m in the  representation (2)
!

double precision, allocatable :: ts(:), rs(:), coefs(:)
double complex, allocatable   :: wsave(:), vals(:), ders(:)

data pi  / 3.14159265358979323846264338327950288d0 /
data ima / (0.0d0,1.0d0) /

ntail            = 8
m                = dlambda*R
ifdone           = 0

call elapsed(t1)

do while (ifdone .eq. 0)

allocate(wsave(15*m+1000))
allocate(ts(2*m+1))
allocate(rs(2*m+1))
allocate(vals(2*m+1))
allocate(ders(2*m+1))
allocate(coefs(-m:m))

rs = R
dd = 2*pi/(2*m+1.0d0)
do i=1,2*m+1
ts(i) = dd* (i-1)
end do


call zffti(2*m+1,wsave)
call wavefun(dlambda,2*m+1,rs,ts,vals,ders,userptr)
call zfftf(2*m+1,vals,wsave)
call zfftf(2*m+1,ders,wsave)

vals    = vals*1.0d0/(2*m+1.0d0)
ders    = ders*1.0d0/(2*m+1.0d0)


do i=1,m+1
coefs(i-1) = abs(vals(i))
end do

do i=1,m
coefs(-i)    = abs(vals(2*m+2-i))
end do


coefs = coefs / maxval(coefs)
dmax1 = maxval(coefs(-m:-m+ntail))
dmax2 = maxval(coefs(m-ntail:m))
dmax  = max(dmax1,dmax2)


ifdone = 1
if (dmax .gt. eps) then
ifdone = 0
m = m*1.1d0
endif

deallocate(wsave,ts,rs,vals,ders,coefs)

end do


end subroutine


end module
