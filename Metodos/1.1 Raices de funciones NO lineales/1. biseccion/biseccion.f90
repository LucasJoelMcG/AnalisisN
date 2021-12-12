PROGRAM biseccion 

!Lo unico a modificar es la funcion (tambien en el script!) y la tolerancia!

Implicit none
Real a,b,tol,m,e
Integer n,i,iter
character op

m=0
call grafica
Write(*,*) 'Ingrese extremos del intervalo'
Read(*,*) a,b
Write(*,*) 'Seleccione: 1.error en x o 2.error en y'
Read(*,*) op

!Ingreso tolerancia de error
tol=0.00001

n=floor((log((b-a)/tol))/(log(2.0)) + 0.5)
iter=1
call system("clear")
Select case(op)
   case('1') 
      Write(*,*) 'Orden: iter, a, b, f(b)-f(a), m, f(m), error'
      Write(*,'(I3)',advance='no') iter
      Write(*,'(6F9.5)') a,b,abs(f(b)-f(a)),(a+b)/2.0,f((a+b)/2.0),(b-a)/2.0
      Do i=1,n
        iter=iter+1
        call raiz(a,b,m)
        Write(*,'(I3)',advance='no') iter 
        Write(*,'(6F9.5)') a,b,abs(f(b)-f(a)),m,f(m),(b-a)/2.0
      End do 
      Write(*,*) ' '
      Write(*,*) 'La raiz se encuentra en x= ',m
      Write(*,*) 'Con un error en x de ',(b-a)/2.0
      Write(*,*) 'Con un error en y de ',f(m)  
   case('2') 
      e=10
      Write(*,*) 'Orden: iter, a, b, f(b)-f(a), m, error'
      Write(*,'(I3)',advance='no') iter
      Write(*,'(5F9.5)') a,b,abs(f(b)-f(a)),(b+a)/2.0,f((b+a)/2.0)
      Do while(e> tol)
        iter=iter+1
        call raiz(a,b,m)
        e=abs(f(m))
        Write(*,'(I3)',advance='no') iter
        Write(*,'(5F9.5)') a,b,abs(f(b)-f(a)),m,f(m)
      End do   
      Write(*,*) ' '
      Write(*,*) 'La raiz se encuentra en x= ',m
      Write(*,*) 'Con un error en y de ',e       
      Write(*,*) 'Con un error en x de ',(b-a)/2.0
End select
Write(*,*) ' '
Write(*,*) 'El abs de la diferencia DEBE converger a cero. De forma contraria se encontro una singularidad'
CONTAINS

Function f(x)
Real x,f

f=log(x**2+1)-exp(x/2)*cos(3.1416*x)
End function f

Subroutine raiz(aa,bb,mm)
Real aa,bb,mm
  m=(a+b)/2.0
  If (f(a)*f(m)>0) then
    a=m
  Else
    b=m
  End if
  m=(a+b)/2.0
End subroutine raiz

Subroutine grafica
Open(unit=3,file='script.p')
Write(3,*) 'set   autoscale                        # escala los ejes automaticamente'
Write(3,*) 'unset log                              # quita la escala logaritmica (si la hubiera)'
Write(3,*) 'unset label                            # quita los titulos anteriores'
Write(3,*) 'set xtic auto                          # establece automaticamente las divisiones del eje x'
Write(3,*) 'set ytic auto                          # establece automaticamente las divisiones del eje y'
Write(3,*) 'set grid'
Write(3,*) 'set title "GRAFICO"'
Write(3,*) 'set xlabel "x"'
Write(3,*) 'set ylabel "y"'
Write(3,*) 'set xrange[-5:6]'

     !IMPORTANTE: CAMBIAR FUNCION SIN QUITAR LOS PARENTESIS EXTERNOS
Write(3,*) 'plot  (y=log(x**2+1)-exp(x/2)*cos(3.1416*x)) title "f(x)" with lines lt 3 ,\'  
Write(3,*) '(y=0)' 
close(3)
call system ("gnuplot -persist script.p")
End subroutine
END PROGRAM

