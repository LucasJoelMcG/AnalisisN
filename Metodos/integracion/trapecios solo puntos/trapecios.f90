PROGRAM trapecios

Implicit none
Integer m
Real,allocatable::x(:),y(:)
Real h,area

!cantidad de puntos
m=9
m=m-1
Allocate(x(0:m),y(0:m))
!Ingrese puntos
x(0)=0.0;   y(0)=0.0
x(1)=0.1;   y(1)=2.1220
x(2)=0.2;   y(2)=3.0244
x(3)=0.3;   y(3)=3.2568
x(4)=0.4;   y(4)=3.1399
x(5)=0.5;   y(5)=2.8570
x(6)=0.6;   y(6)=2.514
x(7)=0.7;   y(7)=2.1369
x(8)=0.8;   y(8)=1.8358

h=x(1)-x(0)
area=h*(y(0)+2*sum(y(1:m-1)) +y(m))/2.0
Write(*,*) 'El area total es ',area
CONTAINS

END PROGRAM
