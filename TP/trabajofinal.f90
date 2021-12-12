program trabajofinal
implicit none
real(8), dimension(63) :: Ui
    integer :: ii
    real(8),parameter:: rho=1130.,cp=3883.,k=0.4324, h=25. ! rho=densidad de los champiñones, cp, calor especifico de los champiñones, k=conductividad termica champiñones, h=coeficiente conveccion aire
    real(8) :: ri, dxi, dti,toli,tfinali,ta,ag
    !tomo como el alto de la lata= 12,4 cm
    dxi = 0.002  ! en metros
    dti=10.    ! en segundos, para t=20 r=0.49273090 estable limite aprox, para t=6.4 aprox para r=1/6
    ri=(k*dti)/(cp*rho*(dxi**2))    ! no tine unidades
    toli=0.01
 
    Ui(1)=118.   ! temperatura que necesito llegar en todo la lata para que la agaritina muera, el autoclave mantiene esa temperatura constante en el fondo, en °C
    do ii=2, 63
      Ui(ii)=20.   !temperatura ambiente, todas las demas zonas que no estan en la base se encuentran a temperatura ambiente, en °C
    end do
  !#### Resolucion ###  
    call ParabolicasExp(Ui,ri,dti,dxi,toli,tfinali)
  !#### Muestra de los valores ###  
 write(*,'(A26)')'Los datos del problema son' 
 write(*,*)
 write(*,'(A30,F6.1,A7)')'Densidad de los champignones= ',rho,' kg/m^3'
 write(*,'(A38,F6.1,A7)')'Calor especifico de los champignones= ',cp,' J/kg*K' 
 write(*,'(A43,F7.5,A5)')'Conductividad termica de los champignones= ',k,' W/m*K'
 write(*,'(A36,F3.0,A9)')'Coeficiente de conveccion del aire= ',h,' W/m^2*K'
 write(*,'(A13,F5.3,A2)')'Avance en x= ',dxi,' m'
 write(*,'(A13,F6.0,A2)')'Avance en t= ',dti,' s'
 Write(*,'(A34,F10.8)')'Con estos valores r queda igual a ',ri
 write(*,'(A57)')'La conserva tiene que llegar a una temperatura de 118 C'
 write(*, '(A88)') 'La conserva llega a la temperatura (con una tolerancia de 0.01) en todo el alto luego de '
 write(*,'(A15,F7.0,A9)')'               ',tfinali,' segundos'
write(*,*)
!#### Destruccion de agaritina ####
!n(t)=n0*e^K/t  pero si la paso a base 10 tengo
!n(t)=n0*10^-t/D donde D=480
ta=0.
dti=20.     
ag=1.!N(t)/n0 en t=0 es igual a 1
open(2,file='nn0.txt')
do while (ag>=1e-12)
  ag=10**(-ta/480.)  !n(t)/n0=10^-t/D
  write(2,'(2F25.15)')ta,ag   
  ta=ta+dti
end do
write(*,'(A30)')'Destruccion de agaritina'
write(*,'(A30)')'N(T)/N0         Tiempo'
write(*,'(E15.10,F25.15)')ag,ta-dti
!#### grafico de n/n0 vs tiempo###
call system ("gnuplot.exe -persist nn0.p")
!grafico de las temperaturas en el tiempo segun la posicion
call system ("gnuplot.exe -persist temp.p")

contains
 subroutine ParabolicasExp(U, r, dt, dx, tol,tfinal)
        implicit none
        real(8), intent(inout) :: U(:)
        real(8), intent(in) :: r, dt, dx, tol      
        real(8) :: t, x, UAnt(size(U))
        integer :: i, n
        real(8),intent (out)::tfinal
 !#### Escritura de los primeros puntos ####        
open(2,file='datos.txt')         
        n = size(U)
        t = 0
        x = 0  
        do i=1, n
        write(2, '(F15.3 F15.3  F25.12)') t,x, U(i)
           x=x+dx
        end do
        write(2, *)
        t=t+dt
 !#### Comienzo del loop ####      
        do while (118.-minval(U)>tol) 
            UAnt = U
            x=0
            do i=2, n-1   !recorre el vector             
                U(i) = r * (UAnt(i+1) + UAnt(i-1)) + (1. - 2. * r) * UAnt(i)
            end do
   !##### Parte de conveccion del aire ######         
            if (U(n)<118.) then                                  
                U(n) = (U(n-1) + h * dx *118./k) / (1 + (h * dx)/k)     
            else 
                U(n) = (U(n-1) - h * dx *118./k) / (1 - (h * dx)/k)
            end if
   !##### escritura de los valores de temperatura #### 
            do i=1, n
            write(2, '(F15.3 F15.3 F25.12)') t,x, U(i)
              x=x+dx
            enddo
            write(2,*)
            t = t + dt
        end do
        tfinal=t-dt  
        close(2)
    end subroutine ParabolicasExp
end program
