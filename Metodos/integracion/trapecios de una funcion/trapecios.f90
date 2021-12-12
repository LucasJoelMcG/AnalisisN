PROGRAM trapecios

Implicit none
Integer m,i
Real,allocatable:: x(:),y(:)
Real area,h,xo
real, PARAMETER :: pi=3.14159265359

!Cantidad de trapecios
m=16

Allocate(x(0:m),y(0:m))

!Ingrese extremos intervalo
x(0)=-1.0;  x(m)=1.0

h=(x(m)-x(0))/m
xo=x(0)
Do i=0,m
   y(i)=f(xo)
   xo=xo+h
End do
area=(h*(y(0)+2*sum(y(1:m-1))+y(m)))/2.0
Write(*,*) 'El area es ',area

CONTAINS

Function f(xi)
Real xi,f

f=(1/(sqrt(2.0*pi)))*exp((-xi**2.0)/2.0)
End function
END PROGRAM
