program testit
use ieee_arithmetic
implicit double precision (a-h,o-z)

double precision :: x
real*10          :: x10
real*16          :: x16


print *,log10(HUGE(x))
print *,log10(HUGE(x10))
print *,log10(HUGE(x16))


end program