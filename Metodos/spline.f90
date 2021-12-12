PROGRAM spline_cubico

!CAMBIAR EN EL SCRIPT  SI HAY QUE CALCULAR INTERSECCION CON OTRA CURVA


IMplicit none
Integer ni,ii
Real,allocatable::xi(:),yi(:),hi(:),bi(:),cji(:),ui(:),di(:),li(:),dji(:),bji(:)

!Cantidad de puntos
ni=7

Allocate(xi(ni),yi(ni),hi(ni),bi(ni),cji(ni),ui(ni),di(ni),li(ni),dji(ni),bji(ni))
!Ingrese puntos
xi(1)= 0.;   yi(1)= 19.33;
xi(2)= 9.35;   yi(2)= 35.64;
xi(3)= 18.35;   yi(3)= 30.88;
xi(4)= 30.75;   yi(4)= 21.66;
xi(5)= 41.75;  yi(5)= 35.6;
xi(6)= 58.92;   yi(6)= 8.31;
xi(7)= 79.47;  yi(7)= 30.12;
!xi(8)= 52.58;  yi(8)= 38.37;
!xi(9)= 63.45;  yi(9)= 20.11;
!xi(10)= 75.14;  yi(10)= 34.32;
!xi(11)= 81.47;  yi(11)= 52.94;
!xi(12)= 97.33;  yi(12)= 40.67;


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!calcula los hi
hi=0.0
Do ii=1,ni-1
   hi(ii)=xi(ii+1)-xi(ii)
End do

!Armo sist tridiagonal
di(1)=1.0
di(ni)=1.0
Do ii=2,ni-1
   di(ii)=2*(hi(ii-1)+hi(ii))  
End do
li(1)=0.0;li(ni)=0.0
Do ii=2,ni-1
   li(ii)=hi(ii-1)
End do
ui(ni)=0.0; ui(1)=0.0
Do ii=2,ni-1
   ui(ii)=hi(ii)
End do

cji=0
call arma_b(ni,hi,bi,yi)
call thomas(ui,di,li,ni,cji,bi)
Do ii=1,ni-1
   dji(ii)=(cji(ii+1)-cji(ii))/(3.0*hi(ii))
   bji(ii)=(yi(ii+1)-yi(ii))/hi(ii) - hi(ii)*(2*cji(ii)+cji(ii+1))/3.0
End do
Write(*,*) 'i      cj         dj          bj'
Write(*,*) ' '
Do ii=1,ni-1
   Write(*,'(I2, A2, F9.5, A2, F9.5, A2, F9.5)') ii,'  ',cji(ii),'  ',dji(ii),'  ',bji(ii)
End do
Write(*,*) ' '
call guarda_datos(xi,yi,ni)
call grafica(ni,cji,bji,dji,xi)
call aproxima(ni,bji,cji,dji,xi,yi)

CONTAINS
Subroutine arma_b(n,h,b,y)
Integer n,i
Real h(n),b(n),y(n)

b(1)=0.0; b(n)=0.0
Do i=2,n-1
   b(i)=3*(y(i+1)-y(i))/h(i) - 3*(y(i)-y(i-1))/h(i-1)
End do
End subroutine

Subroutine thomas(u,d,l,n,x,b)
Integer n,i
Real u(n),d(n),l(n),x(n),b(n)

Do i=1,n-1
   u(i)=u(i)/d(i)
   b(i)=b(i)/d(i)
   d(i)=1.0
   d(i+1)=d(i+1)-l(i+1)*u(i)
   b(i+1)=b(i+1)-l(i+1)*b(i)
   l(i+1)=0.0
End do
b(n)=b(n)/d(n)
Do i=n-1,1,-1
   b(i)=b(i)-u(i)*b(i+1)/d(i)
End do
x=b
End subroutine

Subroutine guarda_datos(x,y,n)
Integer i,n
Real x(n),y(n)

Open(unit=2,file='datos.dat')
Do i=1,n
   Write(2,*) x(i),y(i)
End do
Close(2)
End subroutine

Subroutine grafica(n,cj,bj,dj,x)
Integer i,n
Real cj(n),bj(n),dj(n),h,xo,yo,x(n),y(n)

Write(*,*) 'Elija un paso h para graficar'
Read(*,*) h
Open(unit=4,file='valores.dat')
Do i=1,n-1
   xo=x(i)
   Do while(xo<=x(i+1))
      yo=y(i)+bj(i)*(xo-x(i))+cj(i)*((xo-x(i))**2)+dj(i)*((xo-x(i))**3)
      Write(4,*) xo,yo
      xo=xo+h
   End do
End do
Close(4)
Call system("gnuplot -persist script.p")
End subroutine

Subroutine aproxima (n,bj,cj,dj,x,y)
Integer n,i
Real xo,bj(n),cj(n),dj(n),p,x(n),y(n),pcomprobado

Write(*,*) 'Ingrese xo donde se quiere aproximar la funcion'
Read(*,*) xo
i=0

Do while(xo>x(i+1))
   i=i+1
   p=y(i)+bj(i)*(xo-x(i))+cj(i)*((xo-x(i))**2)+dj(i)*((xo-x(i))**3)
End do
     write(*,*) 'p',i,'(x)=',y(i),'+(',bj(i),')*','(x-',x(i),')+(',cj(i),')*((x-',x(i),')**2)+',dj(i),'*((x-',x(i),')**3)' 
Write(*,*) 'El resultado es ',p


!pcomprobado=   80.00+(8.65836143E-02)*(xo-5.00)+(-1.31673574)*((xo-5.0)**2)+0.207069308*((xo-5.000)**3)

!write(*,*) 'El resultado comprado esssss: ',pcomprobado
End subroutine
END PROGRAM
