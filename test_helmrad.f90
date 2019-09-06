module test_helmrad_functions
use utils
use chebyshev
use pbessel
use helmrad
integer :: iwhich

contains


subroutine potfun(ising,n,rs,qs,userptr)
implicit double precision (a-h,o-z)
double precision          :: rs(n), qs(n), qders(n)
type(c_ptr)               :: userptr


! Gaussian bump
if (iwhich .eq. 1) then
qs           = exp(-5*rs**2)
endif

! volcano
if (iwhich .eq. 2) then 
qs = 14*rs**2 * exp(-5*rs**2)
endif

! discontinuous potential
if (iwhich .eq. 3) then 

qs = 0
if (ising .eq. 1) qs = 1.0d0
if (ising .eq. 3) qs = 2.0d0

endif

end subroutine


subroutine wavefun(dlambda,n,rs,ts,vals,ders,userptr)
implicit double precision (a-h,o-z)
double precision          :: ts(n),rs(n)
double complex            :: vals(n), ders(n)
type(c_ptr)               :: userptr
double complex            :: ima, z, h0, h1, z0
data                      ima / (0.0d0,1.0d0) /
data                      pi  / 3.14159265358979323846264338327950288d0 /

complex*32                :: zz1,zz2

! plane wave at angle t0 
t0    = pi /  4
vals  = exp(ima * rs * dlambda * cos(ts-t0) )
ders  = exp(ima * rs * dlambda * cos(ts-t0) ) * ima * dlambda * cos(ts-t0)


! ! circular wave centered at (x0,y0)
! x0 = 0
! y0 = 6

! do i=1,n
! r       = rs(i)
! t       = ts(i)
! x       = r*cos(t)
! y       = r*sin(t)
! z0      = x0 + ima*y0
! z       = x  + ima*y

! z       = dlambda*abs(z-z0)
! ifexpon = 1
! call hank103(z,h0,h1,ifexpon)

! vals(i) = h0
! ders(i) = dlambda*h1*(-r+x0*cos(t) + y0*sin(t))/sqrt((x0-x)**2+(y0-y)**2)

! end do

end subroutine


subroutine image_incident(ipart,filename,R,dlambda,wfun,userptr)
implicit double precision (a-h,o-z)
character(len=*)             :: filename
procedure(helmrad_wavefun)   :: wfun
type(c_ptr)                  :: userptr
!
!  Produce an image of the incident field.
!
!  Input parameters:
!    ipart - an integer parameter indicating 
!      ipart = 0  means plot the absolute value of the incoming wave
!      ipart = 1  means plot the real part of the incoming wave
!      ipart = 2  means plot the imaginary part of the incoming wave
!    filename - the name of the PDF which is to be generated
!    R - the radius for the plot
!    dlambda - the wavenumber for the incident field
!    wavefun - the external subroutine conforming to the helmrad_wavefun
!      interface which supplies the 
!    userptr - a "void *" pointer passed on to the external subroutine
!

double complex, allocatable   :: uin(:,:),vals(:),ders(:)
double precision, allocatable :: amatr(:,:) ,xs(:),ys(:),rs(:),ts(:)

nn = dlambda*10
nn = max(nn,512)
nn = min(nn,4096)


allocate(uin(nn,nn),amatr(nn,nn),xs(nn),ys(nn))
allocate(rs(1),ts(1),vals(1),ders(1))

do i=1,nn
x      = -R+ (i-1.0d0)/(nn-1.0d0) * 2*R
xs(i) = x
end do

do i=1,nn
y      = -R+ (i-1.0d0)/(nn-1.0d0) * 2*R
ys(i)  = y
end do


do i=1,nn
do j=1,nn
x     = xs(i)
y     = ys(j)

rs(1) = sqrt(x**2+y**2)
ts(1) = atan2(y,x)

call  wfun(dlambda,1,rs,ts,vals,ders,userptr)

uin(nn-j+1,nn-i+1) = vals(1)

end do
end do

if (ipart .eq. 0) amatr = abs(uin)
if (ipart .eq. 1) amatr = real(uin)
if (ipart .eq. 2) amatr = imag(uin)

call show_matrix(0,"","",-R,R,-R,R,filename,amatr)

end subroutine


subroutine image_output(helmdata,scatdata,ipart,filename1,filename2,R,dlambda,wfun,userptr)
implicit double precision (a-h,o-z)
character(len=*)             :: filename1,filename2
procedure(helmrad_wavefun)   :: wfun
type(c_ptr)                  :: userptr
type(helmrad_data)           :: helmdata
type(helmrad_scat)           :: scatdata

!
!  Produce images of the total field and the scattered field.
!
!  Input parameters:
!    ipart - an integer parameter indicating 
!      ipart = 0  means plot the absolute value of the incoming wave
!      ipart = 1  means plot the real part of the incoming wave
!      ipart = 2  means plot the imaginary part of the incoming wave
!    filename1 - the name of the PDF which will contain the image of the total field
!    filename2 - the name of the PDF which will contain the image of the scattered field
!    R - the radius for the plot
!    dlambda - the wavenumber for the incident field
!    wavefun - the external subroutine conforming to the helmrad_wavefun
!      interface which supplies the 
!    userptr - a "void *" pointer passed on to the external subroutine
!

double complex, allocatable   :: uscat(:,:),utot(:,:),vals(:),ders(:)
double precision, allocatable :: amatr(:,:) ,xs(:),ys(:),rs(:),ts(:)
double complex                :: val_tot,val_scat

nn = dlambda*10
nn = max(nn,512)
nn = min(nn,4096)

allocate(utot(nn,nn),uscat(nn,nn),amatr(nn,nn),xs(nn),ys(nn))
allocate(rs(1),ts(1),vals(1),ders(1))

do i=1,nn
x      = -R+ (i-1.0d0)/(nn-1.0d0) * 2*R
xs(i) = x
end do

do i=1,nn
y      = -R+ (i-1.0d0)/(nn-1.0d0) * 2*R
ys(i)  = y
end do


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,x,y,rs,ts,val_tot,val_scat)
!$OMP DO
do i=1,nn
do j=1,nn
x     = xs(i)
y     = ys(j)

rs(1) = sqrt(x**2+y**2)
ts(1) = atan2(y,x)

call helmrad_eval(helmdata,scatdata,wavefun,userptr,rs(1),ts(1),val_tot,val_scat)

if (isNaN(real(val_tot))) then
print *,rs(1),ts(1),x,y,val_tot,val_scat
stop
endif

utot(nn-j+1,nn-i+1)  = val_tot
uscat(nn-j+1,nn-i+1) = val_scat

end do
end do
!$OMP END DO
!$OMP END PARALLEL

if (ipart .eq. 0) amatr = abs(utot)
if (ipart .eq. 1) amatr = real(utot)
if (ipart .eq. 2) amatr = imag(utot)
call show_matrix(0,"x","y",-R,R,-R,R,filename1,amatr)

if (ipart .eq. 0) amatr = abs(uscat)
if (ipart .eq. 1) amatr = real(uscat)
if (ipart .eq. 2) amatr = imag(uscat)
call show_matrix(0,"x","y",-R,R,-R,R,filename2,amatr)

end subroutine


subroutine plot_q(filename,R,xsings,pfun,userptr)
implicit double precision (a-h,o-z)
character(len=*)              :: filename
double precision              :: xsings(:)
type(c_ptr)                   :: userptr
procedure(helmrad_potfun)     :: pfun

!
!  Plot the function q
!

double precision, allocatable  :: xs(:,:),ys(:,:),ab(:,:)
double precision               :: rs(1),qs(1)

nn     = 100
nsings = size(xsings)
nints  = nsings+1

allocate(ab(2,nints),xs(nn,nints),ys(nn,nints))

int     = 1
ab(1,int) = 1.0d-3
do i=1,nsings
ab(2,int) = xsings(i)
int       = int+1
ab(1,int) = ab(2,int-1)
end do
ab(2,int) = R

do int=1,nints
a = ab(1,int)
b = ab(2,int)

do i=1,nn
t     = a + (i-1.0d0)/(nn-1.0d0) * (b-a)
rs(1) = t
call potfun(int,1,rs,qs,userptr)
val   = qs(1)

ys(i,int) = val
xs(i,int) = t

end do
end do

call plot_function(filename,"t","q(t)",nn*nints,xs,ys)

end subroutine



end module

program test_helmrad
use utils
use chebyshev
use pbessel
use helmrad
use iso_c_binding
use test_helmrad_functions

implicit double precision (a-h,o-z)

type(helmrad_data)            :: helmdata
type(helmrad_scat)            :: scatdata

type(c_ptr)                   :: userptr
double precision, allocatable :: xsings(:)
double complex                :: val_tot, val_scat

ima     = (0.0d0,1.0d0)
pi      = acos(-1.0d0)

dlambda = 2**11
R       = 4                ! radius of  the scatterer
R2      = 8                ! dimensions of the region to plot

m       = pi/2*R*dlambda
eps     = 1.0d-12
ifout   = 1
ipart   = 1

iwhich = 3


! if (dlambda .lt. 32) m = m + 20



if (iwhich .eq. 3) then

allocate(xsings(3))
xsings(1) = 1.0d0
xsings(2) = 2.0d0
xsings(3) = 3.0d0
else

allocate(xsings(0))
endif

! call plot_q("q.pdf",R,xsings,potfun,userptr)
! call image_incident(1,"incoming.pdf",R2,dlambda,wavefun,userptr)

call helmrad_init(ifout,eps,m,xsings,R,dlambda,potfun,userptr,helmdata)
call helmrad_solve(helmdata,wavefun,userptr,scatdata)

rr = R/2
tt = pi/4
call helmrad_eval(helmdata,scatdata,wavefun,userptr,rr,tt,val_tot,val_scat)

print *,val_tot
print *,val_scat

stop


call prinz("scatdata%in_coefs  = ",scatdata%in_coefs)
call prinz("scatdata%tot_coefs  = ",scatdata%tot_coefs)
call prinz("scatdata%scat_coefs = ",scatdata%scat_coefs)



end program
