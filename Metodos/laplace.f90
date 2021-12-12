Program ejercicio1
implicit none
real(8), parameter :: xmax= , ymax= , tol= , 
integer f,c,n
real(8) superior(c),inferior(c),izquierda(f),derecha(f),solucion(n),a(n,n+1)

!cambiar segun nodos    f=cantidad de filas de nodos   c=cantidad de columnas de nodos   n=total de nodos
f=
c=
n=f*c
!---------------------------------------------------------------------------------------------------------

contains

!--------------------------------- CONDICIONES DE CONTORNO ---------------------------------

subroutine bordes(superior, inferior,izquierda,derecha,f,c)
integer f,c
real(8) superior(c),inferior(c),izquierda(f),derecha(f)

superior=
inferior=
izquierda=
derecha=
!si las condiciones de contorno son distintas, cambiar
end subroutine

!------------------------------- ARMA LA MATRIZ DE BANDAS ---------------------------------

subroutine armamatrix

a=0

!diagonal principal
do i=1, n    
   a(i,i)= 4.0 
end do
!diagonal superior larga
do i=1, n-1
  if (mod(n,i)=/0) then
    a(i,i+1)= -1.0
  end if
end do
!diagonal superior corta
do i=1 , n-c
  a(i,c+i)= -1.0
end do
!diagonal inferior larga
do i=1, n-1
  if (mod(n,i)=/0) then
    a(i+1,i)= -1.0
  end if
end do
!diagonal inferior corta
do i=1, n-c
  a(i+c,i)= -1.0
end do
!terminos independientes
do i=1, c
  a(i,n+1)= a(i,n+1) + superior(i)
end do
do i=n, n-c,-1
  a(i,n+1)= a(i,n+1) + inferior(i)
end do
do i=1, f
  a(i*c,n+1)= a(i*c,n+1) + derecha(i)
end do
do i=1, f
  a(f*i-f+i,n+1)= a(i*f-f+i,n+1) + izquierda(i)
end do

end subroutine

!---------------------------------- NORMA VECTOR ---------------------------------------

function normavector(v,n)
integer n,i
real(8) v(n),n,normavector

aux=0
do i=1, n
  if (abs(v(i))>aux) then
    normavector=aux
  end if
end do
end function

!-------------------------------------- w ------------------------------------------------

function factor(f,c)
real(8) w,factor
integer f,c

w= 4.0/(2.0 + SQRT(4.0 - (cos(3.1416/(f-1)) + cos(3.1416/(c-1)))**2))
factor= w
end function

!------------------------ CALCULA LA MATRIZ CON GAUSS SEIDEL -----------------------------

subroutine gaussseidel

iter=0
do while (tol<=error)
  x=xnuevo
  
!-gauss seidel--------                    
  do i=1, n
    xnuevo(i)=a(i,n+1)        !xnuevo=terminos indep
    do j=1, i-1
      xnuevo(i)= xnuevo(i) - a(i,j)*x(j)      
    end do
    do j=i+1, n
      xnuevo(i)= xnuevo(i) - a(i,j)*x(j)
    end do
    xnuevo(i) = xnuevo(i)/a(i,i)    
  end do
!-fin gauss seidel-------
  
  !sobrerelajacion
  w=factor(f,c)
  error= normavector(xnuevo-x,n)/normavector(x,n)  !detencion de proceso iterativo 1
  iter=iter + 1

end subroutine

!------------------------------------- GRAFICA --------------------------------------------

subroutine grafica(n,f,c,izquierda,derecha,superior,inferior,solucion)
integer n,f,c,aux
real(8) x,y,h,k
real(8) izquierda(f),derecha(f),superior(c),inferior(c),solucion(n),esquina(4)

open(unit=2, file='datos.dat', status='replace')

!para los bordes

esquina(1)= (superior(1) + izquierda(1))/2.0       !arriba izquierda
esquina(2)= (superior(c) + derecha(1))/2.0         !arriba derecha
esquina(3)= (izquierda(f) + inferior(1))/2.0       !abajo izquierda
esquina(4)= (derecha(f) + inferior(c))/2.0         !abajo derecha

h= xmax/(c+1)
k= ymax/(f+1)

!fila superior (son las condiciones de contorno)
x=0.0
y=ymax
write(2,'(3f10.5)') x,y,esquina(1)
do i=1, c
  x=x+h
  write(2,'(3f10.5)') x,y,superior(i)
end do
x=x+h
write(2,'(3f10.5)') x,y,esquina(2)
write(2,*)

!centro
aux=0
Do i=1, f
  x=0.0
  y=y - k
  write(2, '(3f10.5)') x,y,izquierda(i)
  do j=1, c
    x=x + h
    write(2, '(3f10.5)') x,y,solucion(aux+j)
  end do
  aux=aux+c
  x=x + h
  write(2, '(3f10.5)') x,y,derecha(i)
  write(2,*)
end do

!fila inferior (son las condiciones de contorno)
x=0.0
y=0.0
write(2,'(3f10.5)') x,y,esquina(3)
do i=1, c
  x=x+h
  write(2,'(3f10.5)') x,y,inferior(i) 
end do
x=x+h
write(2,'(3f10.5)') x,y,esquina(4)
write(2,*)
close(2, status='keep')
call system("C:\gnuplot\bin\gnuplot.exe -persist temperatura.p")

End subroutine

End program
