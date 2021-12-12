! Entrada= 
!          * n (numero de puntos dato)
!          * Subroutine Ingreso de datos
!
! Salida= 
!          * Matriz de Resolucion
!          * Polinomio por pantalla
!
Program Polinomio_interpolante
 Implicit None
!-----------------------------
 Integer(4), parameter :: n=2 !Numero de puntos dato
!-----------------------------
 !Integer(4), parameter :: gradopol=n-1
 Real(8) coef(0:n-1), a(n,n+1),f(n),x(n)

!-----------------------------
!Secciòn ejecutable
!-----------------------------
 a=0
 Call ingresodatos
 Call creomatriz
 Call solucionsistema
 Call escriberesultado

Contains

 Subroutine ingresodatos

x(1)=0.0 ; f(1)=1.0
x(2)=0.5  ; f(2)=0.87758
!x(3)=0.72 ; f(3)=0.5
!x(4)=1.5  ; f(4)=-0.23
!x(5)=2.7  ; f(5)=4.23
!x(6)=3.8  ; f(6)=4.02
!x(7)=4.3  ; f(7)=1.2
 EndSubroutine
 
 Subroutine creomatriz
 Integer(4) i,j
 
 a(:,1)=1
 
  Do i=1,n
    Do j=2,n
      a(i,j)=x(i)**((j-1)*1.0)
    EndDo
    a(i,n+1)=f(i)
  EndDo
  
  Do i=1,n
    Do j=1,n+1
       write(*,'(F15.6)', advance = 'No')a(i,j)
    EndDo
    write(*,'(F15.6)', advance = 'yes')
  EndDo
  
  write(*,*)

EndSubroutine
 
 Subroutine solucionsistema    !Pivote-Sustitucion Regresiva

 integer(4) fila, i , j , k
 integer(2), parameter :: cantsist=1
 real(8) aux(n+cantsist)

 Do i=1, n
 
 !Analizo fila a usar de pivote----------------------------------------------
  do k=1, (n+cantsist)
    aux(i)=a(i,k)
  enddo
 
  Do j=i, n
    if (abs(a(j,i))>abs(aux(i))) then
       do k=1, (n+cantsist)
          aux(k)=a(j,k)
          a(j,k)=a(i,k)
          a(i,k)=aux(k)    
       enddo
    endif
  enddo

 
 !----------------------Sustitucion Regresiva---------------------------

  a(i,i+1:)=a(i,i+1:)/a(i,i)
  a(i,i)=1
  Do fila=1, i-1
   a(fila,i+1:)=a(fila,i+1:)-a(i,i+1:)*a(fila,i)
   a(fila,i)=0
  EndDo
  Do fila=i+1, n
   a(fila,i+1:)=a(fila,i+1:)-a(i,i+1:)*a(fila,i)
   a(fila,i)=0
  EndDo

 EndDo
 

 
 Do i=1,n
  coef(i-1)=a(i,n+1)
 EndDo

 EndSubroutine 
 
 Subroutine escriberesultado
 
 Integer i
 
 write(*,*) 'Los coeficientes estan escritos de la sig manera: a0,a1,a2...,an-1,an'
 write(*,*) 'Polinomio= a0 x**0 + a1 x**1 + ... + an x**n'
 write(*,*)
 Do i=0 , n-1
  write(*,*) 'A',i,'=',coef(i)
 endDo
 write(*,*) 'Error para 1: ',abs((coef(0)+(coef(1)*0.5))-cos(0.5))
 EndSubroutine
EndProgram











