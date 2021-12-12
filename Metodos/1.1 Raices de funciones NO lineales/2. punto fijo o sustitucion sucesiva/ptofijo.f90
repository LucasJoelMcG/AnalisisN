PROGRAM puntofijo

!!!!Lo unico a cambiar es f(x), g(x) (tambien en el script!!) y la tolerancia!!!!!!!!!!!!!!!!!!!!

Implicit none
Real tol,error,x, x0
integer iter,maxiter

call grafica
call system ("gnuplot -persist script.dat")
Write(*,*) 'Ingrese valor Xo'
Read(*,*) x

!Ingreso tol en y
tol=0.0005

maxiter=10000
error=10; iter=0
call system ("clear")
Do while((error>tol) .and.(iter<maxiter))  !VER SI LA TOLERANCIA ES EN X O Y --> abs(f(x)) o abs(x-x0)
  x0=x
  iter=iter+1
  x=g(x)
  error=abs(f(x))
 !error=abs(x-x0)
  Write(*,'(I3,3F9.5)') iter,x,f(x),error
End do
Write(*,*) 'Orden: iter, x, f(x), error'
Write(*,*) ' '
Write(*,*) 'La raiz esta en x= ',x
Write(*,*) 'Con un error en y de ',error
Write(*,*) ' '
Write(*,*) 'El x dio cualquier cosa? --> ELEGIR OTRO G(X)'
Write(*,*) ' '
Write(*,*) 'El error DEBE converger SIEMPRE a cero. Sino hay que elegir otro g(x)'
CONTAINS

Function g(x)
Real x,g

g=(x**2.)
End function g

Function f(x)
Real x,f

f= (x**2.)-x
End function

Subroutine grafica

Open(unit=2, file='script.dat')
Write(2,*) 'set   autoscale'                        !# escala los ejes automaticamente
Write(2,*) 'unset log'                              !# quita la escala logaritmica (si la hubiera)
Write(2,*) 'unset label'                            !# quita los titulos anteriores
Write(2,*) 'set xtic auto'                          !# establece automaticamente las divisiones del eje x
Write(2,*) 'set ytic auto'                          !# establece automaticamente las divisiones del eje y
Write(2,*) 'set grid'
Write(2,*) 'set title "Elegir punto semilla"'
Write(2,*) 'set xlabel "x"'
Write(2,*) 'set ylabel "y"'
Write(2,*) 'set xrange[-2:1]'

         !IMPORTANTE: CAMBIAR F Y G SIN QUITAR LOS PARENTESIS EXTERNOS
Write(2,*) 'plot  (y=x) title "y=x" with lines,\'
Write(2,*)      '(y=(x**2)) title "g(x)" with lines lt 8,\'
Write(2,*)      '(y=(x**2)-x) title "f(x)" with lines,\'
Write(2,*)      '(y=0) title "y=0" with lines'
End subroutine

END PROGRAM
