module test_riccati_functions
use utils
use chebyshev
use riccati

double precision                 :: dnu, dlambda

contains


subroutine qfun(k,ts,qs,userptr)
implicit double precision (a-h,o-z)
integer          :: k
double precision :: ts(k)
double complex   :: qs(k)
type(c_ptr)      :: userptr

double precision             :: qs1(k),qs2(k),phis(k)

phis     = (1-Erf(12*(ts-1.5d0)))/2.0d0

qs1      = dlambda**2*ts**2 + (0.25d0 - dnu**2)/ts**2
qs2      = dlambda**2

qs       = phis*qs1 + (1-phis)*qs2


end subroutine


subroutine qfun2(k,ts,qs,userptr)
implicit double precision (a-h,o-z)
integer          :: k
double precision :: ts(k)
double precision :: qs(k)
type(c_ptr)      :: userptr

double precision             :: qs1(k),qs2(k),phis(k)

phis     = (1-Erf(12*(ts-1.5d0)))/2.0d0

qs1      = dlambda**2*ts**2 + (0.25d0 - dnu**2)/ts**2
qs2      = dlambda**2

qs       = phis*qs1 + (1-phis)*qs2


end subroutine


subroutine turning_point(nroots,roots)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: roots(:)

double precision, allocatable  :: ts(:),vals(:),ders(:)
double precision, allocatable  :: ab(:,:)

a        = 1.0d-30
b        = 1.0d0

nbisect  = 10
nnewt    = 8

maxroots = 30
nn       = maxroots*4

allocate(ders(nn),vals(nn),ts(nn),ab(2,maxroots))
ab = 0

ts = 0
bb = 0.1d0

do i=1,nn/2

dd    = 2.0d0**(-nn/2+i)
ts(i) = a + (bb-a)*dd
end do


do i=nn/2+1,nn
dd = 1.0d0/(nn/2)*(i-nn/2)
ts(i) = bb + (b-bb)*dd
end do


! ts(1) = a
! do i=2,nn-1
! !call random_number(ts(i))
! !ts(i) = 1.0d0/(nn+0.0d0)*i
! end do

! ts(nn) = b

! call quicksort(nn,ts)

! call qfun(nn,ts,vals,ders,userptr)
! vals = vals + (0.250 - dnu**2)/ts**2
! ders = ders - 2 * (0.25d0 - dnu**2)/ts**3

vals = dlambda**2 + (0.250 - dnu**2)/ts**2
ders =  - 2 * (0.25d0 - dnu**2)/ts**3

m = 0
do i=1,nn-1

ifint = 0

if (vals(i)*vals(i+1) < 0) then
m = m +1
if (m .gt. maxroots) then
call prina("find_roots failed")
stop
endif

ab(1,m) = ts(i)
ab(2,m) = ts(i+1)
endif
end do


nroots = m
allocate(roots(nroots))


!
!  Use bisection to refine the brackets a bit
!


do iter=1,nbisect
do int=1,m

a0 = ab(1,int)
b0 = ab(2,int)
c0 = (a0+b0)/2

ts(1) = a0
ts(2) = c0



vals(1:2) = dlambda**2 + (0.250 - dnu**2)/ts(1:2)**2
ders(1:2) = - 2 * (0.25d0 - dnu**2)/ts(1:2)**3


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
!  Use Newton to refine the roots further
!
do int=1,m
a0    = ab(1,int)
b0    = ab(2,int)
t     = (a0+b0)/2

do iter=1,nnewt
ts(1)   = t 
vals(1) = dlambda**2 + (0.250 - dnu**2)/ts(1)**2
ders(1) = - 2 * (0.25d0 - dnu**2)/ts(1)**3


val = vals(1)
der = ders(1)

t = t - val/der
end do

roots(int) = t
end do


end subroutine



end module


program test_riccati
use utils
use chebyshev
use riccati
use iso_c_binding
use test_riccati_functions
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata

double precision, allocatable :: ab2(:,:),ab3(:,:),roots(:)
double complex, allocatable   :: psis(:,:),rs2(:,:),rders2(:,:)
double complex, allocatable   :: rs3(:,:),rders3(:,:)

double precision, allocatable              :: rs4(:,:),rders4(:,:),ab4(:,:)

double complex                :: ra, psival

double precision, allocatable :: ab(:,:),avals(:,:),apvals(:,:),appvals(:,:)
double complex, allocatable   :: rs(:,:),rders(:,:)
double complex                :: rb

type(c_ptr)                   :: userptr
double complex                :: ima,val,rval

ima     = (0.0d0,1.0d0)
k       = 30
pi      = acos(-1.0d0)
call chebexps(k,chebdata)


dlambda =     100.0d0
dnu     =      10.0d0


eps     =   1.0d-13

a       =   1.0d-15
b       =   2.0d0

call turning_point(nroots,roots)
call prini("nroots = ",nroots)
call prin2("roots = ",roots)

tp = roots(1)
call prind("tp = ",tp)


! if (tp .gt. 1.0d0) stop



ra     = 1.0d0
psival = 0.0d0
call riccati_ivp_adap(ier,eps,chebdata,1.0d-20,1.0d-15,qfun,ra,nints2,ab2,rs2,rders2,userptr)


stop

! ra0    = 1.0d0
! call rriccati_ivp_adap(ier,eps,chebdata,1.0d-20,1.0d-15,qfun2,ra0,nints4,ab4,rs4,rders4,userptr)

! do int=1,nints4
! do i=1,k
! print *,(real(rs2(i,int))-rs4(i,int))/(rs4(i,int))
! end do
! end do

! stop



! stop

! !
! !  Solve the Riccati equation backward
! !

! aval    =    0
! rb      =   ima*dlambda
! call riccati_tvp(ier,eps,chebdata,tp,b,qfun,rb,nints,ab,rs,rders,userptr)

! if (ier .ne. 0) then
! call prini("after riccati_tvp, ier = ",ier)
! stop
! endif

! call prin2("after riccati_tvp, ab = ",ab)

! call riccati_phase_forward(chebdata,nints,ab,rs,rders,aval,avals,apvals,appvals)


! !call chebpw_plot(ifshow,"plot_r.pdf",nints,ab,k,chebdata%xs,imag(rs))

! !
! !

! ra     = 0.0d0
! psival = 0.0d0
! call riccati_ivp(ier,eps,chebdata,1.0d-20,1.0d-15,qfun,ra,nints2,ab2,rs2,rders2,userptr)

! ! call prinz("rs2 = ",rs2)
!  ra     = rs2(k,nints2)
! ! print *,ra

! ! stop

! call riccati_ivp(ier,eps,chebdata,1.0d-15,tp,qfun,ra,nints2,ab2,rs2,rders2,userptr)


! psival = 0.0d0
! call riccati_int_forward(chebdata,nints2,ab2,rs2,rders2,psival,psis)
! psis = psis - psis(k,nints2)


! call prin2("after riccati_ivp, ab2 = ",ab2)
! !
! !  Find a1 and a2 such that
! !

! call chebpw_eval(nints,ab,k,chebdata%xs,avals,tp,aval)
! call chebpw_eval(nints,ab,k,chebdata%xs,apvals,tp,apval)
! call chebpw_eval(nints,ab,k,chebdata%xs,appvals,tp,appval)

! call chebpw_ceval(nints2,ab2,k,chebdata%xs,psis,tp,psival)
! call chebpw_ceval(nints2,ab2,k,chebdata%xs,rs2,tp,rval)

! y  = exp(real(psival))
! yp = y * real(rval)


! c1 = y*sqrt(apval)
! c2 = yp/sqrt(apval)+y*appval/(2*apval**(1.5d0))

! a1  = sqrt(c1**2+c2**2)
! a2  = atan2(c1,c2)

! if(a2 .gt. pi) then
! a1 = -a1
! a2 = a2-pi
! endif

! if(a2 .le. 0) then
! a1 = -a1
! a2 = a2 + pi
! endif

! !
! !  Adjust the constants and psi
! !

! t = 1.0d0
! call chebpw_eval(nints,ab,k,chebdata%xs,avals,t,aval)
! call chebpw_eval(nints,ab,k,chebdata%xs,apvals,t,apval)
! val = a1 * sin(aval+a2)/sqrt(apval)

! psis = psis - log(val)
! a1   = a1 / val
! c1   = c1 / val
! c2   = c2 / val

! call prind("a1 = ",a1)
! call prind("a2 = ",a2)

! ! val0 = -0.1008950293211436d0

! ! errrel = (val1-val0)/val0

! ! call prin2("relative error in u = ",errrel)


! t = 0.5d0
! call chebpw_eval(nints,ab,k,chebdata%xs,avals,t,aval)
! call chebpw_eval(nints,ab,k,chebdata%xs,apvals,t,apval)
! val    = a1 * sin(aval+a2)/sqrt(apval)
! val1   = -2.09373093703944368461408074822105562d0
! errrel = (val-val1)/val1
! call prin2("errrel in u = ",errrel)

! ! t = 0.1d0
! ! call chebpw_eval(nints,ab,k,chebdata%xs,avals,t,aval)
! ! call chebpw_eval(nints,ab,k,chebdata%xs,apvals,t,apval)
! ! val = a1 * sin(aval+a2)/sqrt(apval)
! ! print *,val
! ! ! val1 = -1.198798799060278d0
! ! ! errrel = (val-val1)/val1
! ! ! call prin2("errrel in u = ",errrel)


! t = 1.0d-3
! call chebpw_ceval(nints2,ab2,k,chebdata%xs,psis,t,psival)

! val = -2957.39685541535565586915260465030780d0 +   3.14159265358979323846264338327950288d0*ima
! print *,val,psival
! ! ! val0   = exp(psival)
! ! ! val1   = -3.51833290344740738624954508122321862d-15
! ! ! errrel =(val0-val1)/val1 

! ! call prin2("errrel in u = ",errrel)


end program
