PROGRAM nombreProg
  IMPLICIT NONE
! DECLARACION DE VARIABLES
integer yyy, ii,ti,ki,gi,OP
real(8), allocatable, dimension(:) :: ex, yi, xri, coeficientes
real(8) Fm,xs, xxi, yyi, Error,xf,h
real(8), allocatable, dimension(:,:) :: datosi
integer , parameter :: m=4 !cantidad de puntos
integer , parameter :: ni=m-1 ! grado del polinomio

allocate (ex(m),yi(m),coeficientes(0:m-1))
allocate (datosi(0:m-1,2))

ex(1)=1.0 ; yi(1)=0.0
ex(2)=1.3333333 ; yi(2)=0.287682072 
ex(3)=1.6666666 ; yi(3)=0.510825623
ex(4)=2.0  ; yi(4)=0.69314718
!ex(5)=2.0  ; yi(5)=0.4055

do ki=1,m
  datosi(ki-1,1)=ex(ki)
  datosi(ki-1,2)=yi(ki)
enddo

call system ("cls") !limpia la pantalla

yyy=1
DO WHILE (yyy/=4)
write (*,*) ' NECESARIAMENTE HAY QUE PASAR POR LA OPCION 1'
write (*,*)
write (*,*) ' Ingrese 0 si quiere obtener los coeficientes del polinomio'
write (*,*) ' Ingrese 1 si quiere obtener el polinomio evaluado en un X determinado. Tambien te da los Lk evaluados en el punto'
write (*,*) ' Ingrese 2 si quiere plotear'
write (*,*) ' Ingrese 3 si quiere estimar el error de interpolacion'
write (*,*) ' Ingrese 4 si quiere salir'
read (*,*) yyy

select case (yyy)

case(0)

    call lagrange(datosi, coeficientes)
    print*,'El polinomio obtenido es:'
    write(*,'(A2)',advance='NO')'Y='
    do gi=m-1,1,-1
      write(*,'(f10.6,A4,I2,A1)',advance='no') coeficientes(gi),'* x^',gi,'+'
    enddo
    write(*,*) coeficientes(0)

case(1)
!LOS COEFICIENTES DEL POLINOMIO DEBEN SACARSE A MANO USANDO:
!L(k)=Multiplicatoria desde i=0,n de (x-xi)/(xk-xi)
!luego los L(k) se multiplican por el yk dato cada uno
!y luego se suman y resolviendo llego al polinomio de lagrange

!Es la iteracion K, en el punto Xk,Yk se hace la sumatoria con todos los
!xi de todas las iters

call system ("cls")

write (*,*) ' Ingrese el punto donde se evaluara el polinomio'
read (*,*) xxi

!Graba tabla de valores
open (2, file='interp.dat', status='replace')
do ii=1,m
 write (2,'(2F10.6)') ex(ii), yi(ii)
end do
close (2, status='keep')

read(*,*)
call system ("cls")

call PoliLagrange(m,ni,xxi,yyi,ex,yi)
write(*,'(A28,F10.6, A4,F10.6)') 'El valor del polinomio en X=',xxi,' es ',yyi
write(*,'(A1,I1,A1,F10.6,A2,F10.6)') 'P',ni,'(',xxi,')=',yyi
!*************************************************************************
case (2)
open (23, file='poli.dat', status='replace')
print*,'Ingrese un paso para graficar'
read*,h
xf=ex(1)
do while (xf.le.ex(m))
 call PoliLagrange(m,ni,xf,yyi,ex,yi)
 write(23,'(2F10.6)')xf,yyi
 xf=xf+h
end do
close (23, status='keep')
call system ("gnuplot -persist 'script.p'")

!*************************************************************************
case(3) ! ERROR DE INTERPOLACION 
 
  call system ("cls")
 write (*,*) ' Ingrese el numero de puntos a considerar'
 read (*,*) ti
 allocate(xri(0:ti-1))
  call system ("cls")
 write (*,'(A22,I2,A14)') ' Ingrese la derivada F',(ti),' de la funcion'
 read (*,*) Fm
  call system ("cls")
 write (*,*) ' Ingrese los datos en X'
 do ii=0,ti-1
 write (*,'(A17,I2)') ' Ingrese el x', ii
 read (*,*) xri(ii)
 end do
   call system ("cls")
 write (*,*) ' Ingrese el x donde quiere estimar el error'
 read (*,*) xs
   call system ("cls")
   Error=E(xs,ti,xri,Fm)
Print*,'El error es ', Error    !xr son los x(i)
  read (*,*)
end select
END DO

contains

SUBROUTINE lagrange(datos, pol)
! Devuelve el polinomio de lagrange de grado (n-1) 
! como resultado de los productos de todos los (x - xi)
REAL(8), ALLOCATABLE :: polAux(:), xValues(:), polLk(:)
REAL(8) pol(0:), datos(0:, :), denom
INTEGER i, n, g, k

n = SIZE(datos, DIM=1)-1
ALLOCATE(polLk(0:n), polAux(0:n), xValues(0:n-1))

pol = 0
DO k=0, n

! Elimina el elemento k de la lista de valores de x
  xValues(:k-1) = datos(:k-1, 1)
  xValues(k:) = datos(k+1:, 1)

! Arma el primer monomio para comenzar los productos
  polAux(1) = 1.0
  polAux(0) = -xValues(0)
  denom = datos(k, 1) -xValues(0)

  DO g=1, n-1
    polLk(g+1) = 1.0
    DO i = g, 1, -1
      polLk(i) = polAux(i-1) - polAux(i)*xValues(g) !para hallar el coeficiente correspondiente al grado i, sumamos el de i-1 con el producto del i por el termino indep(del q agregamos(x-xk))
    END DO
    polLk(0) = -polAux(0)*xValues(g)
    polAux = polLk
    denom = denom*(datos(k, 1) -xValues(g))
  ENDDO
  polLk = polLk*datos(k, 2)/denom

!Suma cada polinomio polLk
  pol = pol + polLk
ENDDO

END SUBROUTINE

function E(x,t,xr,Fm1)
real(8) E, x, v, Fm1
real(8) xr(0:t-1)
integer i,t,fact,k

fact=1.0
do i=1,t
fact=fact*i
end do


v=1.0
 Do k=0,t-1
  v = v*(x-(xr(k)))
 end do 
print*,'Productoria ',v
print*,'Derivada ', Fm1
print*,'factorial ',fact

E= (v*Fm1)/fact

end function E

Subroutine PoliLagrange(iv,n,xx,yy,X,Y)  
  integer i,iv,k,n
  real(8)  X(iv), Y(iv)
  real(8)  xx,yy, U, S
  ! Se fija si el punto esta dentro del intervalo de interpolacion
  if ((xx<X(1)).or.(xx>X(iv))) then
   print *,' El punto elegido esta fuera del intervalo de interpolacion'
  else
  
  !interpolacion  
  yy=0
  do k = 1, n+1
  S=1  
  U=1
    do i = 1, n+1
      if (i/=k) then
      S=S*(X(k)-X(i))
      U=U*(xx-X(i))
      end if
    end do
    
    yy=yy+(U/S)*Y(k)
    !print*, 'El coeficiente L ',(k-1) ,'es ', (U/S) 
  end do
  !print*,'OJO!No son los coeficientes del polinomio!'
  end if
end subroutine


END PROGRAM


