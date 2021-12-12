PROGRAM romberg
!=====================================================================================================
!                         >--IMPORTANTE--<
! La cantidad de puntos puede que no cumpla el 2**N +1. Eso quiere decir
! que SOBRAN puntos!!! En ese caso, calcular con romberg con la cant de puntos que cumpla la condicion
! y luego SUMAR la aproximacion (por cualquier metodo) de los puntos sobrantes
!=====================================================================================================
Implicit none
Integer i,j,m,g,k
Real,allocatable::x(:),y(:),mat(:,:)
Real h,xo,aux,hfijo,suma

!Ingrese cantidad de puntos
g=5
!cant trapecios totales
g=g-1
!Ingrese dimension de la tabla. Despejar m de: 2**m=g (m empieza en 0). OJO! pueden sobrar puntos
m=2
Allocate(mat(0:m,m+1),x(0:g),y(0:g))

!Ingrese puntos
x(0)=0.0;   y(0)=0.9162
x(1)=0.25;  y(1)=0.8109
x(2)=0.5;   y(2)=0.6931
x(3)=0.75;  y(3)=0.5596
x(4)=1.0;   y(4)=0.4055
!x(5)=2.25;  y(5)=-19.0482
!x(6)=2.5;   y(6)=-12.1088
!x(7)=2.75;  y(7)=-6.27875
!x(8)=3.0;   y(8)=-2.11742

mat=0; hfijo=x(1)-x(0)
!armo primera columna
mat(0,1)=((x(g)-x(0))*(y(g)+y(0)))/2.0
Do i=1,m
   mat(i,1)=mat(i-1,1)
   h=(x(g)-x(0))/(2**i)
   suma=0
   Do j=1,2**i -1,2
      aux=x(0)+h*j
      k=0;xo=x(0)
      Do while(xo<=aux)
         xo=xo+hfijo
         k=k+1
      End do
      suma=suma+y(k-1) 
   End do
   mat(i,1)=(mat(i-1,1)+((x(g)-x(0))/(2**(i-1)))*suma)/2.0
End do
Do j=2,m+1
   Do i=0,m+1-j
      mat(i,j)=(4**(j-1)*mat(i+1,j-1)-mat(i,j-1))/(4**(j-1)-1)
   End do
End do

Write(*,*) 'La tabla de Romberg es '
Write(*,*) ' '
Do i=0,m
   Do j=1,m+1
      Write(*,'(F10.5)',advance='no') mat(i,j)
   End do
   Write(*,*) ' '
End do

CONTAINS

END PROGRAM

