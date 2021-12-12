PROGRAM romberg

Implicit none
Integer n,i,j,m,g
Real,allocatable::x(:),y(:),mat(:,:)
Real h,xo

!Ingrese grado del polinomio interpolante (es la dimension de la tabla que se forma)
m=5
g=m-1
g=int(2**g)
Allocate(mat(m,m),x(g+1),y(g+1))

!Ingrese extremos del intervalo
x(1)=1.0; x(g+1)=2.0

mat=0
!armo primera columna
Do i=1,m
   xo=x(1)
   y=0
   n=int(2.0**(i-1))!nro trapecios
   h=(x(g+1)-x(1))/n
   Do j=1,n+1
      y(j)=f(xo)
      xo=xo+h
   End do
   mat(i,1)=(h*(y(1)+2*sum(y(2:n))+y(n+1)))/2.0
End do
Do j=2,m
   Do i=1,m+1-j
      mat(i,j)=(4**(j-1)*mat(i+1,j-1)-mat(i,j-1))/(4**(j-1)-1)
   End do
End do

Write(*,*) 'El ultimo valor es ',mat(1,m)

CONTAINS
Function f(x)
Real x,f

f=log10(1+x)/x
End function
END PROGRAM
