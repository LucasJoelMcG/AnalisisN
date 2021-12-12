PROGRAM min_cuadrados

Implicit none
Integer mi,ni,ii
Real,allocatable:: mati(:,:),bi(:),xi(:),yi(:),ai(:),auxmi(:,:),auxbi(:)
character op

!======================================cantidad de puntos===============================
!sigma^2 = VARIANZA ; RMS = Error medio cuadrático
!                SI NO LO CAMBIAS NO ANDAAAAAAAAAAAAAAAAA ===============================
mi=4

!========================================================================================
Allocate(xi(mi),yi(mi))

!Ingrese puntos
xi(1)= 0.0;   yi(1)=58.0;
xi(2)= 20.;   yi(2)=15.0;
xi(3)= 55.;  yi(3)=68.;
xi(4)= 82.;   yi(4)=35.;
!xi(5)= 1.0;   yi(5)=2.7183;
!xi(6)= 150.0;   yi(6)= 9.;
!xi(7)= 300.;   yi(7)= 7.8;
!xi(8)= 400.;  yi(8)= 6.5;
!xi(9)= 550.;  yi(9)= 4.;

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

op='1'
Do while(op=='1')
   Write(*,*) 'Ingrese grado del polinomio aproximante   n='
   Read(*,*) ni
   Allocate(mati(0:ni,0:ni),bi(0:ni),ai(0:ni),auxmi(0:ni,0:ni),auxbi(0:ni))
   call arma_mat(mati,xi,mi,ni)
   call arma_b(bi,yi,xi,mi,ni)
   auxmi=mati; auxbi=bi
   call gauss(auxmi,auxbi,ni)
   call resuelve_gauss(auxmi,auxbi,ni,ai)
   Write(*,*) 'Los coeficientes son: '
   Write(*,*) ' '
   Do ii=0,ni
      Write(*,*) 'a',ii,':  ',ai(ii)
   End do
   Write(*,*) ' '
   call varianza(xi,yi,ni,mi,ai)
   Write(*,*) ' '
   Write(*,*) 'Seleccione una opcion: '
   Write(*,*) '1: Elegir otro n'
   Write(*,*) '2: Graficar esa funcion'
   Read(*,*) op
   If (op=='1') then
      deallocate(mati,bi,ai,auxmi,auxbi)
   End if
   call system("cls")
End do
call grafica(xi,yi,mi,ai,ni)

CONTAINS

Subroutine arma_mat(mat,x,m,n)
Integer i,j,k,m,n
Real x(m),mat(0:n,0:n)

mat=0
Do i=0,n
   Do j=0,n
      Do k=1,m
         mat(i,j)=mat(i,j)+x(k)**(i+j)
      End do
   End do
end do
End subroutine

Subroutine arma_b(b,y,x,m,n)
Integer i,k,m,n
Real b(0:n),y(m),x(m)

b=0
Do i=0,n
   Do k=1,m
      b(i)=b(i)+y(k)*x(k)**i
   End do
End do
End subroutine

Subroutine gauss(aux_m,aux_b,n)
Integer n,pi,k
Real aux_m(0:n,0:n),aux_b(0:n)

Do pi=0,n
   aux_b(pi)=aux_b(pi)/aux_m(pi,pi)
   aux_m(pi,:)=aux_m(pi,:)/aux_m(pi,pi)
   Do k=pi+1,n
      aux_b(k)=aux_b(k)-aux_m(k,pi)*aux_b(pi)
      aux_m(k,:)=aux_m(k,:)-aux_m(k,pi)*aux_m(pi,:)
   End do
End do
End subroutine

Subroutine resuelve_gauss(m,b,n,a)
Integer n,i,k
Real m(0:n,0:n),b(0:n),a(0:n)

a(n)=b(n)/m(n,n)
Do i=n-1,0,-1
   Do k=i+1,n
      a(i)=a(i)+m(i,k)*a(k)
   ENd do
   a(i)=(b(i)-a(i))/m(i,i)
End do
End subroutine

Function p(xo,coef,n)
Integer n,i
Real xo,coef(0:n),p

p=0
Do i=0,n
   p=p+coef(i)*(xo**i)
End do
End function

Subroutine varianza(x,y,n,m,a)
Integer m,k,n
Real y(m),x(m),var,rms,a(0:n),yo

var=0; yo=0
Do k=1,m
   var=var+((y(k)-p(x(k),a,n))**2)
End do
rms=sqrt(var/m)
var=var/(m-n-1)
Write(*,*) 'La varianza es: ',var
Write(*,*) 'RMS es: ',rms
End subroutine

Subroutine grafica(x,y,m,a,n)
Integer m,n,i
Real x(m),y(m),h,a(0:n),xo

Open(unit=2,file='datos.dat')
Do i=1,m
   Write(2,*) x(i),y(i)
End do
Close(2)
Write(*,*) 'Elija un paso h para graficar'
Read(*,*) h
xo=x(1)
Open(unit=4,file='puntos.dat')
Do while(xo<=x(m)) ! PONER xo<=x(m) SI SE QUIERE QUE SEA EL ÚLTIMO PUNTO QUE SE TIENE DATO, SINO PONER EL VALOR QUE UNO NECESITE
   Write(4,*) xo,p(xo,a,n)
   xo=xo+h
End do
Write(4,*) xo,p(xo,a,n)
Close(4)
Call system('gnuplot -persist script1.p')
End subroutine
END PROGRAM
