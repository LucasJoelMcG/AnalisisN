PROGRAM simpson

Implicit none
Integer m,i,n
Real,allocatable::x(:),y(:)
Real h,area
real, PARAMETER :: pi=3.14159265359

!cantidad de parabolas. Tiene que ser PAR
n=8
!cant subintervalos
m=2*n

Allocate(x(0:m),y(0:m))

!Ingrese limites del intervalo
x(0)=-1.0; x(m)=1.0
!----------------------------------------------------------------------
h=(x(m)-x(0))/m
do i=1,m-1
   x(i)=x(i-1)+h
End do
Do i=0,m
   y(i)=f(x(i))
End do
area=y(0)+4*y(1)+y(m)
Do i=2,m-2,2
   area=area+4*y(i+1)+2*y(i)
End do
area=(h/3.0)*area
Write(*,*) 'El area total es ',area

CONTAINS
Function f(xi)
Real xi,f

f=(1/(sqrt(2.0*pi)))*exp((-xi**2.0)/2.0)
End function

END PROGRAM
