module expr_discont_functions

use utils
use iso_c_binding
use helmrad

contains

subroutine potfun(ising,n,rs,qs,userptr)
implicit double precision (a-h,o-z)
double precision          :: rs(n), qs(n), qders(n)
type(c_ptr)               :: userptr


qs = 0
if (ising .eq. 1) qs = 1.0d0
if (ising .eq. 3) qs = 2.0d0

end subroutine


subroutine wavefun(dlambda,n,rs,ts,vals,ders,userptr)
implicit double precision (a-h,o-z)
double precision          :: ts(n),rs(n)
double complex            :: vals(n), ders(n)
type(c_ptr)               :: userptr
double complex            :: ima
data                      ima / (0.0d0,1.0d0) /
data                      pi  / 3.14159265358979323846264338327950288d0 /

! plane wave at angle t0 
t0    = pi /  4
vals  = exp(ima * rs * dlambda * cos(ts-t0) )
ders  = exp(ima * rs * dlambda * cos(ts-t0) ) * ima * dlambda * cos(ts-t0)



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

call show_matrix(0,"x","y",-R,R,-R,R,filename,amatr)

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


end module


program expr_discont
use utils
use chebyshev
use helmrad
use helmop
use expr_discont_functions

implicit double precision (a-h,o-z)

type(helmrad_data)              :: helmdata
type(helmrad_scat)              :: scatdata
type(c_ptr)                     :: userptr

double precision, allocatable   :: coefs(:), rs(:), ts(:),qs(:)
double precision, allocatable   :: dinits(:),dsolves(:),derrs(:),dlambdas(:),dlogs1(:),dlogs2(:)
integer, allocatable            :: ms(:)

double complex                  :: val,ima,val0,val_tot,val_scat
double complex, allocatable     :: x(:,:), y(:,:)

double precision, allocatable   :: xsings(:),xs0(:),ys0(:),xsave(:,:),amatr(:,:),zs0(:)

eps0 = epsilon(0.0d0)
pi   = acos(-1.0d0)
ima  = (0.0d0,1.0d0)

R   = 4

jj1 = 4
jj2 = 8
jj3 = 17

allocate(xsings(3))
xsings(1) = 1
xsings(2) = 2
xsings(3) = 3

!
!  Extended precision precompuation
!

if (eps0 .lt. 1.0d-16) then

iw = 20
open(iw,FILE='discont.dat')

do jj=jj1,jj2

dlambda = 2.0d0**jj
eps     = eps0*100
m       = R*dlambda*pi
ifout   = 1

! produce a solution

call prin2("dlambda = ",dlambda)

call elapsed(t1)
call helmrad_init(ifout,eps,m,xsings,R,dlambda,potfun,userptr,helmdata)
call helmrad_solve(helmdata,wavefun,userptr,scatdata)
call elapsed(t2)

tsolve = t2-t1

! verify the accuracy of the obtained solution
call elapsed(t1)

dmax = 0.0d0
do ii=1,1

lambda  = dlambda
n        = m

if (jj .eq. 4) l = 40
if (jj .eq. 5) l = 60
if (jj .eq. 6) l = 80
if (jj .eq. 7) l = 140
if (jj .eq. 8) l = 200
if (ii .eq. 0) l = l *2

R1   =  0.0d0+1.0d0*ii
R2   =  1.000+1.0d0*ii

if (ii .eq. 0) R1 = 1.0d-3

allocate(x(2*n+1,l+1),y(2*n+1,l+1),qs(l+1))
call helmop_init(R1,R2,n,l,ts,rs)


ncount = 0
call elapsed(tt1)

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,i,rr,tt,val_tot,val_scat)
!!$OMP DO SCHEDULE(DYNAMIC)
do i=1,2*n+1
do j=1,l+1
tt = ts(i)
rr = rs(j)
call helmrad_eval(helmdata,scatdata,wavefun,userptr,rr,tt,val_tot,val_scat)
x(i,j) = val_tot
end do


!!$OMP CRITICAL
call elapsed(tt2)
ncount=ncount+1
write (*, '("     ",I5.5,"/",I5.5," ",D14.5,A)',advance='no')  ncount,2*n+1,tt2-tt1,13
!!$OMP END CRITICAL

end do
!!$OMP END DO
!!$OMP END PARALLEL

call elapsed(tt2)

ising=4
if (ii .eq. 0) ising = 1
if (ii .eq. 1) ising = 2
if (ii .eq. 2) ising = 3

call potfun(ising,l+1,rs,qs,userptr)
call helmop_helmholtz(dlambda,R1,R2,n,l,ts,rs,qs,x,y)

dmax0 = maxval(abs(y))
dmax  = max(dmax,dmax0)

write (*,*)  R1,R2,dmax0,tt2-tt1
write (13,*) R1,R2,dmax0,tt2-tt1

deallocate(x,y,qs)


end do

call elapsed(t2)

tverify = t2-t1


! sample the solution and output the samples to the disk
call elapsed(t1)
nn = 50
allocate(x(nn,nn),xs0(nn),ys0(nn))
ncount = 0


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,i,xx,yy,rr,tt,val_tot,val_scat)
!$OMP DO SCHEDULE(DYNAMIC)
do i=1,nn
xx     = -R + 2*R * (i-1.0d0) / (nn-1.0d0)
xs0(i) = xx

do j=1,nn
yy     = -R + 2*R * (j-1.0d0) / (nn-1.0d0)
ys0(j) = yy

rr  = sqrt(xx**2+yy**2)
tt  = atan2(yy,xx)
call helmrad_eval(helmdata,scatdata,wavefun,userptr,rr,tt,val_tot,val_scat)

x(i,j) = val_tot

end do

!$OMP CRITICAL
ncount=ncount+1
write (*,'("     ",I5.5,"/",I5.5,A)',advance='no')  ncount,nn,13
!$OMP END CRITICAL

end do
!$OMP END DO
!$OMP END PARALLEL


write (iw,"(D24.16)")   dlambda
write (iw,"(I9.9)")     nn
write (iw,"(D24.16)")   xs0
write (iw,"(D24.16)")   ys0
write (iw,"(D24.16)")   x


deallocate(x,xs0,ys0)

call elapsed(t2)

toutput = t2-t1

lambda =  dlambda
write (*,"(I8.8,' ',D10.3,' ',D10.3, ' ',D10.3, ' ',D10.3, ' ',D10.3)") &
  lambda,tsolve,tverify,toutput,dmax

write (13,"(I8.8,' ',D10.3,' ',D10.3, ' ',D10.3, ' ',D10.3, ' ',D10.3)") &
  lambda,tsolve,tverify,toutput,dmax

write (*,*) ""
write (13,*) ""

end do

close(iw)

stop

endif

!
!  Produce images of the incident field, the total field and the scatterd field
!


call prina("producing images at 16 wavelengths")

dlambda = 16
eps     = 1.0d-12
R       = 4
ifout   = 1
m       = R*dlambda*pi/2
R2      = 2*R
ipart   = 1              ! display real parts

call image_incident(ipart,"discont_incident.pdf",2*R,dlambda,wavefun,userptr)
call helmrad_init(ifout,eps,m,xsings,R,dlambda,potfun,userptr,helmdata)
call helmrad_solve(helmdata,wavefun,userptr,scatdata)
call image_output(helmdata,scatdata,ipart,"discont_total.pdf","discont_scattered.pdf",&
  2*R,dlambda,wavefun,userptr)


!
!  Plot the function q(t)
!

nn = 1000
mm = nn/4
allocate(xs0(nn),ys0(nn))

a0 = 0
b0 = R

do i=1,nn
xs0(i) = a0 + (i-1.0d0) /(nn-1.0d0) * (b0-a0)
end do


call potfun(1,mm,xs0(0*mm+1:1*mm),ys0(0*mm+1:1*mm),userptr)
call potfun(2,mm,xs0(1*mm+1:2*mm),ys0(1*mm+1:2*mm),userptr)
call potfun(3,mm,xs0(2*mm+1:3*mm),ys0(2*mm+1:3*mm),userptr)
call potfun(4,mm,xs0(3*mm+1:4*mm),ys0(3*mm+1:4*mm),userptr)

iw = 101
open(iw,FILE='discont_q.py')
call pyplot_begin(iw,istatus)

call pyplot_add_function(istatus,1,"",nn,xs0,ys0)
call pyplot_ylabel(istatus,"q(t)")
call pyplot_xlabel(istatus,"t")

write(iw,"(A,I12,A,I12,A)") "plt.ylim(-1,4)"

call pyplot_end(istatus,"discont_q.pdf")
call system("python discont_q.py")
deallocate(xs0,ys0)

!
!  Verification/timing in double precision
!

write (*,*) ""
write (13,*) ""

allocate(dinits(jj1:jj3),dsolves(jj1:jj3),ms(jj1:jj3),derrs(jj1:jj3))
allocate(dlambdas(jj1:jj3),dlogs1(jj1:jj3),dlogs2(jj1:jj3))

iw = 20
open(iw,FILE='discont.dat',STATUS='old')

do jj=jj1,jj3


dlambda = 2.0d0**jj
eps     = 1.0d-12
R       = 4
ifout   = 0
m       = R*dlambda*pi/2

call elapsed(t1)
call helmrad_init(ifout,eps,m,xsings,R,dlambda,potfun,userptr,helmdata)
call elapsed(t2)
tinit  = t2-t1

call elapsed(t1)
call helmrad_solve(helmdata,wavefun,userptr,scatdata)
call elapsed(t2)
tsolve = t2-t1



if (jj .le. jj2) then

read (iw,"(D24.16)")    dlambda0
read (iw,"(I9.9)")      nn

allocate(xs0(nn),ys0(nn),x(nn,nn),y(nn,nn))


read (iw,"(D24.16)")   xs0
read (iw,"(D24.16)")   ys0
read (iw,"(D24.16)")   x


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,xx,yy,rr,tt,val_tot,val_scat)
!$OMP DO
do i=1,nn
do j=1,nn

xx     = xs0(i)
yy     = ys0(j)

rr  = sqrt(xx**2+yy**2)
tt  = atan2(yy,xx)
call helmrad_eval(helmdata,scatdata,wavefun,userptr,rr,tt,val_tot,val_scat)

y(i,j) = val_tot


end do
end do
!$OMP END DO
!$OMP END PARALLEL



derr = maxval(abs(x-y))
deallocate(x,y,xs0,ys0)

else
derr = -1
endif


dlambda0     = 2.0d0**jj1

ms(jj)       = m
dinits(jj)   = tinit
dsolves(jj)  = tsolve
derrs(jj)    = derr
dlambdas(jj) = dlambda

dlogs1(jj)   = dinits(jj1)  *  log(dlambda) * dlambda / (dlambda0*log(dlambda0))
dlogs2(jj)   = dsolves(jj1) *  log(dlambda) * dlambda / (dlambda0*log(dlambda0))

if (jj .eq. jj1) then
ratio1=-1
ratio2=-1
else
ratio1 = dinits(jj)/dinits(jj-1)
ratio2 = dsolves(jj)/dsolves(jj-1)
endif


lambda = dlambda
write (*,"(I8.8,' ',I8.8,' ',D14.5,' ',D14.5,' ',D14.5,' ',D14.5,' ',D14.5)")  &
  lambda,m,tinit,ratio1,tsolve,ratio2,derr
write (13,"(I8.8,' ',I8.8,' ',D14.5,' ',D14.5,' ',D14.5,' ',D14.5,' ',D14.5)")  &
  lambda,m,tinit,ratio1,tsolve,ratio2,derr





end do

close(iw)




!
!  Produce a LaTeX table for the timing
!

iw = 25
open(iw,FILE='discont_table.tex')


write(iw,'(A)') "\begin{tabular}{crccc}"
write(iw,'(A)') "\toprule" 
write (iw,*) "\addlinespace[.25em]"

write(iw,'(A)', advance='no') "$k$                                &"
write(iw,'(A)', advance='no') "$m$                                &"
write(iw,'(A)', advance='no') "Maximum absolute                   &"
write(iw,'(A)', advance='no') "Precomp time                       &"

write(iw,'(A)'              ) "Solve time                        \\"
 
write(iw,'(A)', advance='no') "                                   &"
write(iw,'(A)', advance='no') "                                   &"
write(iw,'(A)', advance='no') "error                              &"
write(iw,'(A)', advance='no') "(in seconds)                       &"
write(iw,'(A)'              ) "(in seconds)                       \\"

write(iw,'(A)')               "\midrule"
write (iw,*) "\addlinespace[.25em]"


do jj=jj1,jj3


if (jj .ge. 10) then
write(iw,"('$2^{',I2,'}$ ')",advance='no') jj
else
write(iw,"('$2^{',I1,'}$')",advance='no') jj
endif

call write_table_next(iw)

write(iw,"(I8)",advance='no') ms(jj)
call write_table_next(iw)

if (derrs(jj) == -1) then
write (iw,"(A)",advance='no') "-"
else
call write_table_double(iw,derrs(jj))
endif

call write_table_next(iw)
call write_table_double(iw,dinits(jj))


call write_table_next(iw)
call write_table_double(iw,dsolves(jj))

call write_table_nextline(iw)

write (iw,*) "\addlinespace[.125em]"

end do


write (iw,'(A)') "\bottomrule"
write (iw,'(A)') "\end{tabular}"
close(iw)



!
! Generate the two timing graphs
!

nn = jj3-jj1+1
iw = 101

open(iw,FILE='discont_graph1.py')

call pyplot_begin(iw,istatus)


call pyplot_add_function(istatus,1,"",nn,dlambdas,dinits)
call pyplot_add_function(istatus,9,"",nn,dlambdas,dlogs1)


write(iw,"(A)") "ax.text(0.90, 0.90, 'k log(k)',"
write(iw,"(A)") "verticalalignment='bottom', horizontalalignment='right',"
write(iw,"(A)") "transform=ax.transAxes,"
write(iw,"(A)") "color='black', fontsize=12)"

call pyplot_xlabel(istatus,"Wavenumber k")
call pyplot_ylabel(istatus,"Precomputation time (seconds)")
call pyplot_xlogscale2(istatus)

call pyplot_end(istatus,"discont_graph1.pdf")
call system("python discont_graph1.py")





open(iw,FILE='discont_graph2.py')

call pyplot_begin(iw,istatus)


call pyplot_add_function(istatus,1,"",nn,dlambdas,dsolves)

call pyplot_xlabel(istatus,"Wavenumber k")
call pyplot_ylabel(istatus,"Solve time (seconds)")
call pyplot_xlogscale2(istatus)

call pyplot_end(istatus,"discont_graph2.pdf")
call system("python discont_graph2.py")





end program
