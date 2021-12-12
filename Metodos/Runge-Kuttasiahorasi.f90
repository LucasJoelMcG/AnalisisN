program Rungekutta1
implicit none
  real(8), DIMENSION(0:1)::v_iniciali
  real(8), parameter:: e=2.71828,tolerancia=0.005,kuno=1.0,kdos=1.0
  real(8) hi
  hi=0.1
  v_iniciali(0)=0
  v_iniciali(1)=1
!call Hfijo(v_iniciali,hi)
call hvariable(v_iniciali,hi)

CONTAINS

  FUNCTION v_prima(v)
  ! Definición de la Ecuacion Diferencial
  REAL(8), DIMENSION(0:1) :: v, v_prima
  v_prima(0) = 1.0
  v_prima(1) = (-v(1))/(1+(v(0))**2.0)
  END FUNCTION v_prima
  
  Function v_real(v)
  !Ecuacion exacta
  real(8), dimension(0:1) :: v
  real(8) v_real
  v_real=5*e**(-20*v(0))+(7*(e**(-0.5*v(0))-e**(-20*v(0))))/(19.5)  
  end function
  
function RungeKutta(v,h)
!Metodo de Runge-Kutta (de 4to. Orden)
REAL(8), INTENT(IN) :: h
REAL(8), DIMENSION(0:1) :: v, k1, k2, k3, k4,rungekutta
k1 = h*v_prima(v)
k2 = h*v_prima(v+k1/2.0)
k3 = h*v_prima(v+k2/2.0)
k4 = h*v_prima(v+k3)
v = v + (k1 + 2.0*k2 + 2.0*k3 +k4)/6.0
rungekutta=v
END function
 
  function Hmodificado(v,h)
  logical resultado
  real(8), DIMENSION(0:1):: y,y_ast,v
  real(8) h,error,hmodificado
  integer i,cont
  
  resultado=.true.
  cont=0
  do while (resultado .and. (cont<100))
     y=v
     y=rungekutta(y,h)
     y_ast=v
     do i=1,2
        y_ast=rungekutta(y_ast,h/2.)
     end do
     error=maxval(abs((y - y_ast)))
     if ((error) > tolerancia) then
        h=h/2.
     else if ((error) < (tolerancia/2.)) then
        h=h*2
     else
        resultado=.false.
     end if
     cont=cont+1
  end do
  write(*,*) cont, error
  hmodificado=h
  end function

SUBROUTINE Hfijo(v_inicial,h)
!Rutina con h fijo
REAL(8), INTENT(IN), DIMENSION(0:1) :: v_inicial
REAL(8), INTENT(IN) :: h
real(8) max_valor,min_valor
INTEGER iter,max_iter
REAL(8), DIMENSION(0:1) :: v
open(unit=2,file='datos.dat')
v = v_inicial
 max_valor=11
 min_valor=v(0)
 max_iter=int(((max_valor-min_valor)/h)+0.5)
WRITE (2, '(I5, 5F20.11)') 0, v,h
DO iter = 1, max_iter
v=rungekutta(v,h)
WRITE (2, '(I5, 5F20.11)') iter, v,h
END DO
close(2)
call system ('gnuplot -persist plot.plt')
END SUBROUTINE 

SUBROUTINE Hvariable(v_inicial,h)
!Rutina con h variable
REAL(8), INTENT(IN), DIMENSION(0:1) :: v_inicial
REAL(8), DIMENSION(0:1) :: v,derivada
real(8) h
integer iter
open(unit=2,file='datos.dat',status='replace')
v = v_inicial
WRITE (2, '(I5, 3F20.10)') 0, v,h!,0.0
iter=0
DO while (v(0)<1150)
iter=iter+1
h=hmodificado(v,h)
v=rungekutta(v,h)
!derivada=v_prima(v)
WRITE (2, '(I5, 3F20.10)') iter, v,h !,derivada(1)
END DO
close(2,status='keep')
call system ('gnuplot -persist plot.plt')
END SUBROUTINE 
end program
