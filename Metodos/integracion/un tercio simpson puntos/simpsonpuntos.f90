PROGRAM simpson

Implicit none
Integer m,i,n
Real,allocatable::x(:),y(:)
Real h,area

!cantidad de puntos. Tiene que ser IMPAR
n=7
!cant subintervalos
m=n-1

Allocate(x(0:m),y(0:m))

!Ingrese puntos
x(0)=1.0; y(0)=-54.4021
x(1)=1.4; y(1)=38.6545
x(2)=1.8; y(2)=-20.5278
x(3)=2.2; y(3)=-20.3739
x(4)=2.6; y(4)=-9.58135
x(5)=3.0; y(5)=-2.11742
x(6)=3.4; y(6)=-1.57301
!----------------------------------------------------------------------
h=x(1)-x(0)
area=y(0)+4*y(1)+y(m)
Do i=2,m-2,2
   area=area+4*y(i+1)+2*y(i)
End do
area=(h/3.0)*area
Write(*,*) 'El area total es ',area

CONTAINS

END PROGRAM

