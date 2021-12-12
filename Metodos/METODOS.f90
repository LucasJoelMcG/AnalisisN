PROGRAM METODOS
IMPLICIT NONE
REAL(8) hi,errtolii,xfinali,Valori
REAL(8), ALLOCATABLE :: vi(:),Aii(:,:),Xi(:),Bi(:),A2i(:,:),Ci(:),AIii(:,:)
INTEGER opcion,ii,ordenii,pasoi,cantecii
 !Filas y columnas de la matriz cuadrada A o tambien uso mi si no es cuadrada:   
INTEGER, PARAMETER :: ni=5, mi=6

WRITE (*,*)
WRITE (*,*)'Tipo de problema a resolver:'
WRITE (*,*)
WRITE (*,*)'1: Sistema de ecuaciones lineales'
WRITE (*,*)'2: Ecuacion diferencial'
WRITE (*,*)
WRITE (*, '(A8)',ADVANCE='no') 'Opcion: '
READ  (*,*) opcion
WRITE (*,*)

IF (opcion==1) THEN

 !VALORES A DEFINIR: #############################################################################



 ALLOCATE(Aii(ni,ni),Xi(ni),Bi(ni),A2i(ni,ni+1),AIii(ni,2*ni))

 !ELEMENTOS DE LA MATRIZ A

 Aii(1,1)=0.0   ; Aii(1,2)=0.0   ; Aii(1,3)=100.0   ; Aii(1,4)=100.0   ; Aii(1,5)=0.0   ;Aii(1,6)=0.0   ;!  Aii(1,7)=0.0    ;  Aii(1,8)=100.0  ! ;  !Aii(1,9)=100.00     ;  !Aii(1,10)=1.3      
 Aii(2,1)=0.0   ; Aii(2,2)=50.0  ; Aii(2,3)=0.0     ; Aii(2,4)=100.0   ; Aii(2,5)=100.0   ;Aii(2,6)=100.0   ;!  Aii(2,7)=0.0    ;  Aii(2,8)=100.0  ! ;  !Aii(2,9)=100.00     ;  !Aii(2,10)=7.01
 Aii(3,1)=0.0   ; Aii(3,2)=0.0   ; Aii(3,3)=0.0    ; Aii(3,4)=0.0     ; Aii(3,5)=0.0   ;Aii(3,6)=100.0   ;!  Aii(3,7)=0.0    ;  Aii(3,8)=100.0  ! ;  !Aii(3,9)=100.00     ;  !Aii(3,10)=2.46
 Aii(4,1)=0.0   ; Aii(4,2)=0.0   ; Aii(4,3)=0.0    ; Aii(4,4)=0.0   ; Aii(4,5)=100.0   ;Aii(4,6)=0.0   ;!  Aii(4,7)=0.0    ;  Aii(4,8)=100.0   !;  !Aii(4,9)=100.00     ;  !Aii(4,10)=0.02
 Aii(5,1)=0.0   ; Aii(5,2)=0.0   ; Aii(5,3)=0.0     ; Aii(5,4)=0.0   ; Aii(5,5)=0.0   ;Aii(5,6)=0.0   ;!  Aii(5,7)=3.10   ;  Aii(5,8)=7.45  ;  Aii(5,9)=100.00       ;  !Aii(5,10)=19.32

 !TERMINOS INDEPENDIENTES

 Bi(1)=6.0
 Bi(2)=25.0
 Bi(3)=-11.0
 Bi(4)=15.0
 !Bi(5)=14.8

 A2i=0.!no modificar
 Xi=0.    !no modificar

 WRITE (*,*)'Seleccionar tipo de operacion a realizar:'
 WRITE (*,*)
 WRITE (*,*)'Metodos DIRECTOS de resolucion del sistema:'
 WRITE (*,*)' 1: Matriz A es diagonal (estuuuuupido)'
 WRITE (*,*)' 2: Matriz A es triangular inferior'
 WRITE (*,*)' 3: Matriz A es triangular superior'
 WRITE (*,*)' 4: Metodo de Gauss'
 WRITE (*,*)' 5: Metodo de Gauss-Jordan (doble Gauss)'
 WRITE (*,*)' 6: Metodo de Thomas'
 WRITE (*,*)' 7: Metodo de Crout, factorizacion de matriz A (Ax = LUx = Lc = B)'
 WRITE (*,*)
 WRITE (*,*)'Metodos INDIRECTOS de resolucion del sistema:'
 WRITE (*,*)' 8: Metodo de Jacobi'
 WRITE (*,*)' 9: Metodo de Gauss-Seidel (Jacobi acelerado)'
 WRITE (*,*)
 WRITE (*,*)'Operaciones:'
 WRITE (*,*)'10: Hallar matriz Inversa de A'
 WRITE (*,*)'11: Hallar numero de condicion de la matriz A'
 WRITE (*,*)'12: Hallar norma Euclidiana del vector C'
 WRITE (*,*)'13: Hallar norma M de la matriz A (maxima suma de abs(filas))'
 WRITE (*,*)'14: Hallar norma L de la matriz A (maxima suma de abs(columnas))'
 WRITE (*,*)'15: Hallar norma de Frobenius de la matriz A'
 WRITE (*,*)
 WRITE (*,*)'16: Gauss-Seidel para elipticas (revisar valores fijos)'
 WRITE (*,*)
 WRITE (*,'(A8)',ADVANCE='no')'Opcion: '
 READ  (*,*)opcion
 WRITE (*,*)

 SELECT CASE(opcion)

  CASE(1)
      CALL Diagonal(Aii,Xi,Bi,ni)

  CASE(2)
      CALL Trianginf(Aii,Xi,Bi,ni)

  CASE(3)
      CALL Triangsup(Aii,Xi,Bi,ni)

  CASE(4)
      WRITE (*,*)'Con mejora iterativa?'
      WRITE (*,*)
      WRITE (*,*)'1: Si'
      WRITE (*,*)'2: No'
      WRITE (*,*)
      WRITE (*,'(A8)',ADVANCE='no')'Opcion: '
      READ  (*,*)opcion
      WRITE (*,*)
      IF (opcion==1) THEN
          CALL GaussMejoraIterativa(Aii,Xi,Bi,ni,A2i)
      ELSE
          CALL Gauss(Aii,Bi,ni,A2i)
          CALL Triangsup(Aii,Xi,Bi,ni)
      END IF

  CASE(5)
      CALL GaussJordan(Aii,Xi,Bi,ni,A2i)

  CASE(6) !Thomas
      CALL Gauss(Aii,Bi,ni,A2i)
      CALL Triangsup(Aii,Xi,Bi,ni)
      WRITE (*,*)

  CASE(7)
      A2i=Aii
      DO ii=1,ni
          A2i(ii,ni+1)=Bi(ii)
      END DO
      CALL Crout(A2i,ni)

  CASE(8)
      CALL Jacobi(Aii,Xi,Bi,ni)

  CASE(9)
      CALL GaussSeidel(Aii,Xi,Bi,ni)

  CASE(10)
      CALL InvertirMatriz(Aii,AIii,ni)
      WRITE (*,*)
      WRITE (*,*)'Corroborar haciendo A x A(-1) ?'
      WRITE (*,*)
      WRITE (*,*)'1: Si'
      WRITE (*,*)'2: No'
      WRITE (*,*)
      WRITE (*,'(A8)',ADVANCE='no')'Opcion: '
      READ  (*,*)opcion
      WRITE (*,*)
      IF (opcion==1) THEN
          CALL CorroborarInversa(Aii,AIii,ni)
      END IF

  CASE(11)
      Valori=Cond(Aii,ni)

  CASE(12)
      ordenii=4                                       !DATOS PARA INGRESAR
      ALLOCATE (Ci(ordenii))
      Ci(1)=1
      Ci(2)=1
      Ci(3)=1
      Ci(4)=1
      Valori=NormaEuclidiana(Ci,ordenii)

  CASE(13)
      Valori=NormaM(Aii,ni)

  CASE(14)
      Valori=NormaL(Aii,ni)

  CASE(15)
      Valori=NormaFrobenius(Aii,ni)

  CASE(16)
    CALL ElipticaGS(Aii)

 END SELECT


end if

!VALORES A DEFINIR: #############################################################################

!ERROR TOLERABLE:  #Solo se usa si el i es variable                  REAL
!errtolii=0.0001

!PASO INICIAL:                                                          REAL
!hi=0.05

!ORDEN DE LA ECUACION DIFERENCIAL                                       INTEGER
!cantecii=2

!ALLOCATE (vi(cantecii+1))

!EL VECTOR "V" ES LAS COMPONENTES X, Y, Y', Y"... ETC HASTA EL ORDEN DE LA E.D.O.

!X INICIAL:                                                             REAL
!vi(1)=0.

!X FINAL:                                                               REAL
!xfinali=3.5

!Y INICIAL:                                                             REAL
!vi(2)=0.

!Y PRIMA INICIAL:                                                       REAL
!vi(3)=0.

!Y SEGUNDA INICIAL:                                                     REAL
!vi(4)=0.4

 !OPEN  (UNIT=2,FILE='Valores.dat',STATUS='REPLACE')
 !WRITE (2,*)'Aca imprime el vector "v" por columnas: x, y, yprima, etc' 
 !WRITE (2,*)
 !WRITE (2,'(F10.6,2X)',ADVANCE='no') vi(1)
 !WRITE (2,'(F10.6,2X)',ADVANCE='no') vi(2)
 !WRITE (2,'(F10.6,2X)',ADVANCE='no') vi(3)
 !WRITE (2,'(F10.6,2X)') vi(4)
 !WRITE (*,*)'Seleccionar metodo de resolucion:'
 !WRITE (*,*)
 !WRITE (*,*)'1: Euler Simple'
 !WRITE (*,*)'2: Euler Modificado'
 !WRITE (*,*)'3: Runge-Kutta de cuarto orden'
 !WRITE (*,*)'4: Runge-Kutta de sexto orden'
 !WRITE (*,*)
 !WRITE (*,'(A8)',ADVANCE='no')'Opcion: '
 !READ  (*,*)opcion
 !WRITE (*,*)
 !WRITE (*,*)'Seleccionar tipo de paso:'
 !WRITE (*,*)
 !WRITE (*,*)'1: Paso constante'
 !WRITE (*,*)'2: Paso variable estrategia "cambio de h"'
 !IF (opcion==4) THEN
 !    WRITE (*,*)'3: Paso variable estrategia "error estimado por RKF"'
 !END IF
 !WRITE (*,*)
 !WRITE (*,'(A8)',ADVANCE='no')'Opcion: '
 !READ  (*,*)pasoi
 !WRITE (*,*)
!
 !SELECT CASE (opcion)
!
 ! CASE(1)
 !     CALL EUSI(hi,vi,pasoi,errtolii,xfinali,cantecii)
!
 ! CASE(2)
 !     CALL EUMO(hi,vi,pasoi,errtolii,xfinali,cantecii)
!
 ! CASE(3)
 !     CALL RUKU4(hi,vi,pasoi,errtolii,xfinali,cantecii)
!
 ! CASE(4)
 !     CALL RUKUFE(hi,vi,pasoi,errtolii,xfinali,cantecii)
!
 !END SELECT
!
 !CLOSE (2,STATUS='KEEP')
 !WRITE (*,*)
 !WRITE (*,*)'Graficar?'
 !WRITE (*,*)
 !WRITE (*,*)'1: Si'
 !WRITE (*,*)'2: No'
 !WRITE (*,*)
 !WRITE (*,'(A8)',ADVANCE='no')'Opcion: '
 !READ  (*,*)opcion
!
 !IF (opcion==1) THEN
 !    CALL SYSTEM("C:\gnuplot\bin\gnuplot.exe -persist scriptMetodos.p ")
!
 !END IF
!
!END IF
!
!WRITE (*,*)
!WRITE (*,*)'\             /\/\/\/\/\/\/\/\/\/\/\            /'
!WRITE (*,*)' >============ Aguante Newton loco  ===========<'
!WRITE (*,*)'/             \/\/\/\/\/\/\/\/\/\/\/            \'


CONTAINS


!============================================================================================
!============================================================================================
!============================================================================================


!RUTINAS PARA PROBLEMAS DE SISTEMA DE ECUACIONES


SUBROUTINE Diagonal(A,X,B,n)
 REAL(8) A(n,n),X(n),B(n)
 INTEGER n,i

 WRITE (*,*)'El vector solucion X es:'
 WRITE (*,*)
 DO i=1,n
     X(i)=B(i)/A(i,i)
     WRITE (*,'(A2)',ADVANCE='no')'X('
     WRITE (*,'(I1)',ADVANCE='no')i
     WRITE (*,'(A3)',ADVANCE='no')')= '
     WRITE(*,'(F16.12)')X(i)

 END DO


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE Trianginf(A,X,B,n)
 REAL(8) A(n,n),X(n),B(n)
 REAL(8) aux
 INTEGER n,i,j

 WRITE (*,*)'El vector solucion X es:'
 WRITE (*,*)
 DO i=1,n
     aux=0
     IF (i>1) THEN
         DO j=1,i-1
             aux=aux+X(j)*A(i,j)
         END DO
     END IF
     X(i)=(B(i)-aux)/A(i,i)
     WRITE (*,'(A2)',ADVANCE='no')'X('
     WRITE (*,'(I1)',ADVANCE='no')i
     WRITE (*,'(A3)',ADVANCE='no')')= '
     WRITE (*,'(F16.12)')X(i)

 END DO

END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE Triangsup(A,X,B,n)
 REAL(8) A(n,n),X(n),B(n)
 REAL(8) aux
 INTEGER n,i,j

 DO i=n,1,-1
     aux=0
     IF (i<n) THEN
         DO j=i,n-1,1
             aux=aux+X(j+1)*A(i,j+1)
         END DO
     END IF
     X(i)=(B(i)-aux)/A(i,i)

 END DO

 WRITE (*,*)'El vector solucion X es:'
 WRITE (*,*)
 DO i=1,n
     WRITE (*,'(A2)',ADVANCE='no')'X('
     WRITE (*,'(I1)',ADVANCE='no')i
     WRITE (*,'(A3)',ADVANCE='no')')= '
     WRITE (*,'(F16.12)')X(i)

 END DO


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE GaussMejoraIterativa(A,X,B,n,A2)
REAL(8) A(n,n),X(n),B(n),A2(n,n+1),r(n),dx(n),normar,tol,suma
INTEGER n,i,j,iter,maxiter

 maxiter=10                                                             !DEFINIR  ==============
 tol=0.001                                                              !DEFINIR  ==============
 normar=10
 iter=0
 CALL Gauss(A,B,n,A2)
 X(n)=A2(n,n+1)

 DO i=n-1,1,-1
     suma=0
     DO j=i+1,n
         suma=suma+A2(i,j)*x(j)
     END DO
     X(i)=A2(i,n+1)-suma
 END DO

 DO WHILE ((normar>tol).AND.(iter<=maxiter))
     r=MATMUL(A,X)-B
     DO i=1,n
         WRITE (*,'(A2,I1,A3,F33.30)')'r(',i,')= ',r(i)
     END DO
     WRITE (*,*)
     normar=MAXVAL(ABS(r))
     A2=A
     DO i=1,n
         A2(i,n+1)=r(i)
     END DO
     iter=iter+1
     IF (normar>tol) THEN
         CALL Gauss(A,B,n,A2)
         dx=0
         dx(n)=A2(n,n+1)
         DO i=n-1,1,-1
             suma=0
             DO j=i+1,n
                 suma=suma+A2(i,j)*dx(j)
             END DO
             dx(i)=A2(i,n+1)-suma
         END DO
         X=X+dx
     END IF
 END DO

 WRITE (*,*)'El vector solucion X es:'
 WRITE (*,*)
 DO i=1,n
     WRITE (*,'(A2,I1,A3,F13.10)')'X(',i,')= ',X(i)
 END DO
 WRITE (*,*)
 WRITE (*,'(A24,I3)')'Iteraciones realizadas: ',iter


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE Gauss(A,B,n,A2)
REAL(8) A(n,n),B(n),A2(n,n+1),aux
INTEGER n,i,j,k

 A2=A
 DO i=1,n
     A2(i,n+1)=B(i)
 END DO

 DO i=1,n
     CALL PivoteoParcialGauss(A2,n,i)
     DO j=i+1,n+1
         A2(i,j)=A2(i,j)/A2(i,i)
     END DO
     A2(i,i)=1
     DO k=i+1,n
         aux=A2(k,i)
         DO j=i+1,n+1
             A2(k,j)=A2(k,j)-(aux*A2(i,j))
         END DO
         A2(k,i)=0
     END DO

 END DO

 DO i=1,n
     DO j=1,n
         A(i,j)=A2(i,j)
     END DO
 END DO

 DO i=1,n
     B(i)=A2(i,n+1)
 END DO


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE GaussJordan(A,X,B,n,A2)
REAL(8) A(n,n),X(n),B(n),A2(n,n+1),aux
INTEGER n,i,j,k

 A2=A
 DO i=1,n
     A2(i,n+1)=B(i)
 END DO

 DO i=1,n
     CALL PivoteoParcialGauss(A2,n,i)
     DO j=i+1,n+1
         A2(i,j)=A2(i,j)/A2(i,i)
     END DO
     A2(i,i)=1
     DO k=i+1,n
         aux=A2(k,i)
         DO j=i+1,n+1
             A2(k,j)=A2(k,j)-(aux*A2(i,j))
         END DO
         A2(k,i)=0
     END DO
 END DO

 DO j=n+1,1,-1
     DO i=j-2,1,-1
         A2(i,n+1)=A2(i,n+1)-(A2(i,j-1)*A2(j-1,n+1))
     END DO
     DO i=j-2,1,-1
         A2(i,j-1)=0
     END DO
 END DO

 DO i=1,n
     X(i)=A2(i,n+1)
 END DO

 WRITE (*,*)'El vector solucion X es:'
 WRITE (*,*)
 DO i=1,n
     WRITE (*,'(A2)',ADVANCE='no')'X('
     WRITE (*,'(I1)',ADVANCE='no')i
     WRITE (*,'(A3)',ADVANCE='no')')= '
     WRITE (*,'(F16.12)')X(i)

 END DO

END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE InvertirMatriz(A,AI,n)
REAL(8) A(n,n),aux,AI(n,2*n),T(n,n),S(n,n),Q(n,n)
INTEGER n,i,j,k

 AI=0
 AI=A
 DO i=1,n
     AI(i,n+i)=1.
 END DO

 DO i=1,n
     CALL PivoteoParcialInvertir(AI,n,i)
     DO j=i+1,2*n
         AI(i,j)=AI(i,j)/AI(i,i)
     END DO
     AI(i,i)=1.
     DO k=i+1,n
         aux=AI(k,i)
         DO j=i+1,2*n
             AI(k,j)=AI(k,j)-(aux*AI(i,j))
         END DO
         AI(k,i)=0.
     END DO

 END DO

 DO i=n,2,-1
     DO k=i-1,1,-1
         AI(k,:)=AI(k,:)-(AI(i,:)*AI(k,i))
     END DO
 END DO

 WRITE (*,*)'Matriz Inversa de A:'
 WRITE (*,*)
 DO i=1,n
     DO j=n+1,2*n
         WRITE (*,'(F16.12)', advance='no')AI(i,j)
     END DO
     WRITE (*,*)
 END DO

 T=0
 S=0
 DO i=1,n
     DO j=1,n
         T(i,j)=A(i,j)
     END DO
 END DO
 DO i=1,n
     DO j=n+1,2*n
         S(i,j-n)=AI(i,j)
     END DO
 END DO

 Q=MATMUL(T,S)


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE CorroborarInversa(A,AI,n)
REAL(8) A(n,n),T(n,n),S(n,n),Q(n,n),AI(n,2*n)
INTEGER n,i,j

 T=0
 S=0
 DO i=1,n
     DO j=1,n
         T(i,j)=A(i,j)
     END DO
 END DO
 DO i=1,n
     DO j=n+1,2*n
         S(i,j-n)=AI(i,j)
     END DO
 END DO

 Q=MATMUL(T,S)

 WRITE (*,*)'El producto vectorial entre A y su inversa hallada es:'
 WRITE (*,*)
 DO i=1,n
     DO j=1,n
         WRITE (*,'(F10.7)',ADVANCE='no')Q(i,j)
     END DO
     WRITE (*,*)
 END DO

END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


Subroutine Crout(A2,n)
REAL(8) A2(n,n+1),L(n,n),U(n,n),X(n),c(n)
INTEGER i,j,k,n
REAL(8) Su,suma

 U=0
 L=0
 DO i=1,n
     L(i,1)=A2(i,1)
 END DO

 DO j=1,n
     U(1,j)=A2(1,j)/L(1,1)
 END DO

 DO i=2,n
     DO j=2, n
         Su=0
         IF (j<=i) THEN
             DO k=1, j-1
                 Su= Su+(L(i,k)*U(k,j))
             END DO
             L(i,j)=A2(i,j)-Su
         ELSE
             DO k=1, i-1
                 Su= Su+(L(i,k)*U(k,j))
             END DO
             U(i,j)=(A2(i,j)-Su)/L(i,i)
         END IF
     END DO
 END DO

 DO i=1, n
     U(i,i)=1
 END DO

 WRITE (*,*)'La matriz L es:'
 DO i=1, n
     DO j=1, n
         WRITE (*,'(F6.2)',ADVANCE='no')L(i,j)
     END DO
     WRITE (*,*)
 END DO

 WRITE (*,*)
 WRITE (*,*)'La matriz U es:'
 DO i=1, n
     DO j=1, n
         WRITE (*,'(F6.2)',ADVANCE='no')U(i,j)
     END DO
     WRITE (*,*)
 END DO
 WRITE (*,*)

 c(1)=A2(1,n+1)/L(1,1)

 DO i=2 , n
     Suma=0
     DO k=1, i-1
         Suma=Suma+L(i,k)*c(k)
     END DO
     c(i)=(A2(i,n+1)-Suma)/(L(i,i))
 END DO

 X(n)=c(n)

 DO i=n-1 , 1 , -1
     Suma=c(i)
     DO k=n-i, 1, -1
         Suma=Suma-U(i,k+i)*x(k+i)
     END DO
     X(i)= Suma
 END DO

 WRITE (*,'(A24)')'El vector solucion X es:'

 DO j=1, n
     WRITE (*,'(A2)',advance='no')'X('
     WRITE (*,'(I1)',advance='no')j
     WRITE (*,'(A3)',advance='no')')= '
     WRITE (*,'(F16.12)')X(j)
 END DO


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE Jacobi(A,X,B,n)
REAL(8) A(n,n),X(n),B(n),Xnueva(n),S
REAL(8) r(n)        !   PARA CRITERIO 2
REAL(8) tol         !   PARA CRITERIOS 2 Y 3
!INTEGER k           !   PARA CRITERIO 1
INTEGER n,i,j,iter

 X=1
 Xnueva=1000

!CRITERIO 1:
!iter=200
!DO k=1,iter    !                                 ORDENAR LAS FILAS ANTES DE INGRESAR LOS DATOS

!CRITERIO 2:
 tol=0.0001 !tolerancia del residuo
 iter=0       !no tocar iter                                       ORDENAR LAS FILAS ANTES DE INGRESAR LOS DATOS
 r=1000      !no tocar r
DO WHILE ((MAXVAL(ABS(r)))>=tol)                               !r=//r//


!CRITERIO 3:
! tol=0.001
! iter=0    !                                      ORDENAR LAS FILAS ANTES DE INGRESAR LOS DATOS
! DO WHILE ((MAXVAL(ABS(Xnueva-X)))>=tol)
!     X=Xnueva


     DO i=1,n
         S=0
         DO j=1,i-1
             S=S+(A(i,j)/A(i,i))*X(j)
         END DO
         DO j=i+1,n
             S=S+(A(i,j)/A(i,i))*X(j)
         END DO
         Xnueva(i)=(B(i)/A(i,i))-S
     END DO
     iter=iter+1                                  !PARA CRITERIOS     2 Y 3
     X=Xnueva                                           !PARA CRITERIOS 1 Y 2
     r=Residuo(A,X,B,n)                               !PARA CRITERIO      2

 END DO

 WRITE (*,'(A24)')'El vector solucion X es:'

 DO j=1, n
     WRITE (*,'(A2,I1,A3,F16.12)')'X(',j,')= ',X(j)
 END DO
 WRITE (*,*)
 WRITE (*,'(A24,I3)')'Iteraciones realizadas: ',iter

END SUBROUTINE

!============================================================================================
!============================================================================================
!============================================================================================

subroutine Datosgrafico(T)
    real(8) T(ni,mi)
    integer i,j
    OPEN(UNIT=2,FILE='datos1.dat',STATUS='REPLACE')
    do i=1, ni
        write(2,'(6F10.4)') T(i,:)
    enddo
    CLOSE(2,STATUS='KEEP')
    CALL SYSTEM ('gnuplot -persist "plotMat.plt"')
end subroutine

!============================================================================================
!============================================================================================
!============================================================================================


subroutine ElipticaGS(u)
    real(8) error, u(ni,mi), uant(ni,mi),tol
    integer i, j, iter, maxiter
    maxiter=2000
    tol=0.0001
    iter=0
    error=2*tol
    do while ((error>tol).and.(iter<maxiter)) !Comienza el ciclo iterativo de Gauss seidel
        uant=u !Guardo u en iteración anterior
        do j=2,ni-1
            do i=2,mi-1
                u(j,i) = (u(j-1,i) + u(j+1,i) + u(j,i-1) + u(j,i+1))/4.
                if(j==2)then
                    u(2,2)=50.; u(2,4)=100. ;u(2,5)=100. !Corrijo nodos interiores conocidos
                end if
                if((j==4).and.(i==5))then
                    u(4,5)=100. !Corrijo nodos interiores conocidos
                end if
            end do
        end do
        error=NormaM(uant-u,ni) !Recalculo el error
        iter=iter+1
        !write(*,*) iter, error
    end do  
    call Datosgrafico(u)
END SUBROUTINE
    
    
    !============================================================================================
    !============================================================================================
    !============================================================================================
    

SUBROUTINE GaussSeidel(A,X,B,n)
REAL(8) A(n,n),X(n),B(n),S
REAL(8) r(n)        !   PARA CRITERIO     2
REAL(8) tol         !   PARA CRITERIOS    2 Y 3
!REAL(8) Xaux(n)     !   PARA CRITERIO        3
!INTEGER k           !   PARA CRITERIO 1
INTEGER n,i,j,iter

 X=0

!CRITERIO 1:
! iter=200
! DO k=1,iter    !                                  ORDENAR LAS FILAS ANTES DE INGRESAR LOS DATOS


!CRITERIO 2:
 tol=0.0001    !tolerancia del resudio tocar
 iter=0         ! ORDENAR LAS FILAS ANTES DE INGRESAR LOS DATOS
 r=1000
 DO WHILE ((MAXVAL(abs(r)))>=tol)                               !r=//r//


!CRITERIO 3:
! tol=0.001
! iter=0                             ORDENAR LAS FILAS ANTES DE INGRESAR LOS DATOS
! Xaux=1000
! DO WHILE ((maxval(abs(Xaux-X)))>=tol)
!     Xaux=X

     iter=iter+1                      !PARA CRITERIOS 2 Y 3
     DO i=1,n
         S=0
         DO j=1,i-1
             S=S+(A(i,j)/A(i,i))*X(j)
         END DO
         DO j=i+1,n
             S=S+(A(i,j)/A(i,i))*X(j)
         END DO
         X(i)=(B(i)/A(i,i))-S
     END DO
     WRITE (*,'(A7,I3)')'Iter = ',iter
     DO j=1, n
         WRITE (*,'(A2,I1,A3,F16.12)')'X(',j,')= ',X(j)
     END DO
     WRITE (*,*)
     r=Residuo(A,X,B,n)               !PARA CRITERIO 2
     READ (*,*)
 END DO
 WRITE (*,'(A24)')'El vector solucion X es:'
 DO j=1, n
     WRITE (*,'(A2,I1,A3,F16.12)')'X(',j,')= ',X(j)
 END DO
 WRITE (*,*)
 WRITE (*,'(A24,I3)')'Iteraciones realiazdas: ',iter
END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE PivoteoParcialGauss(A,n,i)
INTEGER n,h,nuevafila,i
REAL(8) A(n,n+1)
REAL(8) maximo,aux2

 maximo=0
 DO h=i,n
     IF (abs(A(i,h))>maximo) THEN
         maximo=abs(A(i,h))
         nuevafila=h
     END IF
 END DO

 DO h=i,n+1
     aux2=A(i,h)
     A(i,h)=A(nuevafila,h)
     A(nuevafila,h)=aux2
 END DO

END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


SUBROUTINE PivoteoParcialInvertir(A2,n,X)
INTEGER n,h,l,X,nuevafila
REAL(8) A2(n,2*n)
REAL(8) maximo,aux2

 maximo=0
 nuevafila=X
 Aux2=0

 DO h=x,n
     IF (abs(A2(h,x))>maximo) THEN
         maximo=abs(A2(h,x))
         nuevafila=h
     END IF
 END DO

 DO l=1,2*n
     aux2=A2(x,l)
     A2(x,l)=A2(nuevafila,l)
     A2(nuevafila,l)=aux2
 END DO


END SUBROUTINE


!============================================================================================
!============================================================================================
!============================================================================================


FUNCTION Residuo(A,X,B,n)
REAL(8) A(n,n),X(n),B(n),r(n),Residuo(n)
INTEGER n,i

 DO i=1,n
     r=(MATMUL(A,X))-B
 END DO

 WRITE (*,*)'Residuo:'
 DO i=1,n
     WRITE (*,'(F16.12)')r(i)
 END DO
 WRITE (*,*)

 Residuo=r


END FUNCTION


!============================================================================================
!============================================================================================
!============================================================================================


FUNCTION Cond(A,n)
REAL(8) A(n,n),pivote,aux,cte,AI(n,2*n),AI2(n,n),Cond
INTEGER n,i,j,k,l

 AI=0
 AI=A

 DO i=1,n
     AI(i,n+i)=1
 END DO

 DO i=1,n
     CALL PivoteoParcialInvertir(AI,n,i)
     pivote = AI(i,i)
     DO j=1,2*n
         AI(i,j)=AI(i,j)/pivote
     END DO
     DO k=i,n-1
         aux=AI(k+1,i)
         DO j=1,2*n
             AI(k+1,j)=AI(k+1,j)-(aux*AI(i,j))
         END DO
     END DO
 END DO

 DO j=2,n
     DO i=j-1,1,-1
         cte= AI(i,j)
         DO l=1,2*n
             AI(i,l)=AI(i,l)-(cte*AI(j,l))
         END DO
     END DO
 END DO

 DO i=1,n
     DO j=1,n
         AI2(i,j)=AI(i,n+j)
     END DO
 END DO

 WRITE (*,'(A16,F16.12)')'Norma Matriz A: ',maxval(abs(A))
 WRITE (*,*)
 WRITE (*,'(A27,F16.12)')'Norma Matriz inversa de A: ',maxval(abs(AI2))
 WRITE (*,*)
 Cond=(maxval(abs(A)))*(maxval(abs(AI2)))
 WRITE (*,'(A36,F16.12)')'Numero de Condicion de la Matriz A: ',Cond


END FUNCTION


!============================================================================================
!============================================================================================
!============================================================================================


FUNCTION NormaEuclidiana(C,orden)
REAL(8) C(orden),NormaEuclidiana,aux
INTEGER orden,i

 aux=0

 DO i=1,orden
     aux=aux+C(i)**2
 END DO

 NormaEuclidiana=sqrt(aux)
 WRITE (*,'(A19,F16.12)')'Norma Euclidiana = ',NormaEuclidiana


END FUNCTION


!============================================================================================
!============================================================================================
!============================================================================================


FUNCTION NormaM(A,n)
REAL(8) A(n,n),NormaM,aux,auxtot
INTEGER n,i,j

 aux=0
 auxtot=0
 DO i=1,n
     aux=0
     DO j=1,n
         aux=aux+ABS(A(i,j))
     END DO
     IF (aux>auxtot) THEN
         auxtot=aux
     END IF
 END DO

 NormaM=auxtot
 WRITE (*,'(A10,F16.12)')'Norma M = ',NormaM


END FUNCTION


!============================================================================================
!============================================================================================
!============================================================================================


FUNCTION NormaL(A,n)
REAL(8) A(n,n),NormaL,aux,auxtot
INTEGER n,i,j

 aux=0
 auxtot=0
 DO j=1,n
     aux=0
     DO i=1,n
         aux=aux+ABS(A(i,j))
     END DO
     IF (aux>auxtot) THEN
         auxtot=aux
     END IF
 END DO

 NormaL=auxtot
 WRITE (*,'(A10,F16.12)')'Norma L = ',NormaL


END FUNCTION


!============================================================================================
!============================================================================================
!============================================================================================


FUNCTION NormaFrobenius(A,n)
REAL(8) A(n,n),NormaFrobenius,suma
INTEGER n,i,j

 suma=0

 DO i=1,n
     DO j=1,n
         suma=suma+A(i,j)**2
     END DO
 END DO

 NormaFrobenius=SQRT(suma)
 WRITE (*,'(A21,F16.10)')'Norma de Frobenius = ',NormaFrobenius


END FUNCTION


!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================
!============================================================================================


!RUTINAS PARA PROBLEMAS DE ECUACIONES DIFERENCIALES

!
!SUBROUTINE EUSI(h,v,paso,errtol,xfinal,cantec)
!REAL(8) h,v(cantec+1),vaux1(cantec+1),vaux2(cantec+1),errtol,xfinal,hmin
!INTEGER paso,iter,cantec
!CHARACTER(13) formato
!
! WRITE (formato,'(A1,I1,A11)')'(',cantec+1,'(F10.6,2X))'
! iter=0
! SELECT CASE (paso)
!
!  CASE(1)                                                !Paso constante
!     DO WHILE ((v(1)-h)<xfinal)
!         v=v+h*VPRIMA(v,cantec)
!         WRITE(2,formato)v
!         iter=iter+1
!
!     END DO
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!
!  CASE(2)                                                !Paso variable estrategia "cambio de h"
!     vaux1=0
!     vaux2=0
!     DO WHILE ((v(1)-h)<xfinal)
!         vaux1=v+(h/2.)*VPRIMA(v,cantec)
!         vaux1=vaux1+(h/2.)*VPRIMA(vaux1,cantec)
!         vaux2=v+h*VPRIMA(v,cantec)
!         IF (abs((vaux1(2)-vaux2(2)))>errtol) THEN
!             h=(h/2.)
!             WRITE (*,*)'DIVIDE A h POR 2'
!         ELSE
!             IF (abs((vaux1(2)-vaux2(2)))<(errtol/50.)) THEN
!                 h=h*2
!                 WRITE (*,*)'MULTIPLICA A h POR 2'
!             ELSE
!                 v=vaux2
!                 WRITE(2,formato)v
!             END IF
!         END IF
!
!         !SEGUIMIENTO DE x y h
!         WRITE (*,*)
!         WRITE (*,'(A3,F6.3)')'X= ',v(1)
!         WRITE (*,'(A3,F10.8)')'h= ',h
!!         READ (*,*)                    !HABILITAR SI SE QUIERE VER CADA PASO
!         iter=iter+1
!         IF (h<hmin) THEN
!             hmin=h
!         END IF
!
!     END DO
!     WRITE (*,*)
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!     WRITE (*,'(A22,F10.8)')'Menor paso utilizado: ',h
!
! END SELECT
!
!END SUBROUTINE
!
!
!!============================================================================================
!!============================================================================================
!!============================================================================================
!
!
!SUBROUTINE EUMO(h,v,paso,errtol,xfinal,cantec)
!REAL(8) h,v(cantec+1),vaux(cantec+1),vfinal1(cantec+1),vfinal2(cantec+1),errtol,xfinal,hmin
!INTEGER paso,iter,cantec
!CHARACTER(13) formato
!
! WRITE (formato,'(A1,I1,A11)')'(',cantec+1,'(F10.6,2X))'
! iter=0
! SELECT CASE (paso)
!
!  CASE(1)                                                    !Paso constante
!     DO WHILE ((v(1)-h)<xfinal)
!         vaux=v+h*VPRIMA(v,cantec)
!         v=v+(h/2.)*(VPRIMA(v,cantec)+VPRIMA(vaux,cantec))
!         WRITE(2,formato)v
!         iter=iter+1
!
!     END DO
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!
!  CASE(2)                                                    !Paso variable estrategia "cambio de h"
!     DO WHILE ((v(1)-h)<xfinal)
!         vaux=v+(h/2.)*VPRIMA(v,cantec)
!         vaux=v+(h/4.)*(VPRIMA(v,cantec)+VPRIMA(vaux,cantec))
!         vfinal1=vaux
!
!         vaux=vaux+(h/2.)*VPRIMA(vaux,cantec)
!         vfinal1=vaux+(h/4.)*(VPRIMA(vaux,cantec)+VPRIMA(vfinal1,cantec))  !Aca define el 1er valor de "y" a comparar
!
!         vaux=v+h*VPRIMA(v,cantec)
!         vfinal2=v+(h/2.)*(VPRIMA(v,cantec)+VPRIMA(vaux,cantec))           !Aca define el 2do valor de "y" a comparar
!         IF (abs((vfinal1(2)-vfinal2(2)))>errtol) THEN
!             h=(h/2.)
!             WRITE(*,*)'DIVIDE A h POR 2'
!         ELSE
!             IF (abs((vfinal1(2)-vfinal2(2)))<(errtol/100.)) THEN
!                 h=h*2
!                 WRITE(*,*)'MULTIPLICA A h POR 2'
!             ELSE
!                 v=vfinal2
!                 WRITE(2,formato)v
!             END IF
!         END IF
!
!         !SEGUIMIENTO DE x y h
!         WRITE (*,*)
!         WRITE (*,'(A3,F6.3)')'X= ',v(1)
!         WRITE (*,'(A3,F10.8)')'h= ',h
!!         READ (*,*)                    !HABILITAR SI SE QUIERE VER CADA PASO
!         iter=iter+1
!         IF (h<hmin) THEN
!             hmin=h
!         END IF
!
!     END DO
!     WRITE (*,*)
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!     WRITE (*,'(A22,F10.8)')'Menor paso utilizado: ',h
!
! END SELECT
!
!END SUBROUTINE
!
!
!!============================================================================================
!!============================================================================================
!!============================================================================================
!
!
!SUBROUTINE RUKU4(h,v,paso,errtol,xfinal,cantec)
!REAL(8) h,v(cantec+1),k1(cantec+1),k2(cantec+1),k3(cantec+1),k4(cantec+1)
!REAL(8) vfinal1(cantec+1),vfinal2(cantec+1),errtol,xfinal,hmin
!INTEGER paso,iter,cantec
!CHARACTER(13) formato
!
! WRITE (formato,'(A1,I1,A11)')'(',cantec+1,'(F10.6,2X))'
! iter=0
! SELECT CASE (paso)
!
!  CASE(1)                                                !Paso constante
!     DO WHILE ((v(1)-h)<xfinal)
!
!         k1=h*VPRIMA(v)
!         k2=h*VPRIMA(v+k1/2.)
!         k3=h*VPRIMA(v+k2/2.)
!         k4=h*VPRIMA(v+k3)
!
!         v=v+(1./6)*(k1+2*k2+2*k3+k4)
!         WRITE(2,formato)v
!         iter=iter+1
!
!     END DO
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!
!  CASE(2)                                                       !Paso variable estrategia "cambio de h"
!     DO WHILE ((v(1)-h)<xfinal)
!
!         k1=(h/2.)*VPRIMA(v)                           !Usa h=h/2.
!         k2=(h/2.)*VPRIMA(v+k1/2.)
!         k3=(h/2.)*VPRIMA(v+k2/2.)
!         k4=(h/2.)*VPRIMA(v+k3)
!
!         vfinal1=v+(1./6)*(k1+2*k2+2*k3+k4)
!
!         k1=(h/2.)*VPRIMA(vfinal1)                     !Usa h=h/2. de nuevo
!         k2=(h/2.)*VPRIMA(vfinal1+k1/2.)
!         k3=(h/2.)*VPRIMA(vfinal1+k2/2.)
!         k4=(h/2.)*VPRIMA(vfinal1+k3)
!
!         vfinal1=vfinal1+(1./6)*(k1+2*k2+2*k3+k4)      !Aca define el 1er valor de "y" a comparar
!
!         k1=h*VPRIMA(v)                                !Usa h=h
!         k2=h*VPRIMA(v+k1/2.)
!         k3=h*VPRIMA(v+k2/2.)
!         k4=h*VPRIMA(v+k3)
!
!         vfinal2=v+(1./6)*(k1+2*k2+2*k3+k4)            !Aca define el 2do valor de "y" a comparar
!
!         WRITE (*,'(A11,F6.3)')'X final 1= ',vfinal1(2)
!         WRITE (*,'(A11,F6.3)')'X final 2= ',vfinal2(2)
!         READ (*,*)
!         IF (abs((vfinal1(2)-vfinal2(2)))>errtol) THEN
!             h=(h/2.)
!             WRITE (*,*)'DIVIDE A h POR 2'
!         ELSE
!             IF (abs((vfinal1(2)-vfinal2(2)))<(errtol/50.)) THEN
!                 h=h*2
!                 WRITE (*,*)'MULTIPLICA A h POR 2'
!             ELSE
!                 v=vfinal2
!                 WRITE(2,formato)v
!                 iter=iter+1
!             END IF
!
!         END IF
!
!         !SEGUIMIENTO DE x y h
!         WRITE (*,*)
!         WRITE (*,'(A3,F6.3)')'X= ',v(1)
!         WRITE (*,'(A3,F10.8)')'h= ',h
!!         READ (*,*)                    !HABILITAR SI SE QUIERE VER CADA PASO
!         IF (h<hmin) THEN
!             hmin=h
!         END IF
!
!     END DO
!     WRITE (*,*)
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!     WRITE (*,'(A22,F10.8)')'Menor paso utilizado: ',h
!
! END SELECT
!
!END SUBROUTINE
!
!
!!============================================================================================
!!============================================================================================
!!============================================================================================
!
!
!SUBROUTINE RUKUFE(h,v,paso,errtol,xfinal,cantec)
!REAL(8) h,v(cantec+1),k1(cantec+1),k2(cantec+1),k3(cantec+1),k4(cantec+1),k5(cantec+1)
!REAL(8) k6(cantec+1),vfinal1(cantec+1),vfinal2(cantec+1),errtol,xfinal,eRKF(cantec+1),hmin
!INTEGER paso,iter,cantec
!CHARACTER(13) formato
!
! WRITE (formato,'(A1,I1,A11)')'(',cantec+1,'(F10.6,2X))'
! iter=0
! hmin=10
! SELECT CASE (paso)
!
!  CASE(1)                                                !Paso constante
!     DO WHILE ((v(1)-h)<xfinal)
!
!         k1 = h*VPRIMA(v)
!         k2 = h*VPRIMA(v+k1/4.)
!         k3 = h*VPRIMA(v+(3*k1/32.)+(9*k2/32.))
!         k4 = h*VPRIMA(v+(1932*k1/2197.)-(7200*k2/2197.)+(7296*k3/2197.))
!         k5 = h*VPRIMA(v+(439*k1/216.)-(8.*k2)+(3680*k3/513.)-(845*k4/4104.))
!         k6 = h*VPRIMA(v-(8*k1/27.)+(2.*k2)-(3544*k3/2565.)+(1859*k4/4104.)-(11*k5/40.))
!
!         v = v + ((25*k1/216.)+(1408*k3/2565.)+(2197*k4/4104.)-(k5/5.))
!
!         WRITE(2,formato)v
!         iter=iter+1
!
!     END DO
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!
!  CASE(2)                                                !Paso variable estrategia "cambio de h"
!     DO WHILE ((v(1)-h)<xfinal)
!
!         k1=(h/2.)*VPRIMA(v)                                          !Usa h/2
!         k2=(h/2.)*VPRIMA(v+k1/4.)
!         k3=(h/2.)*VPRIMA(v+(3*k1/32.)+(9*k2/32.))
!         k4=(h/2.)*VPRIMA(v+(1932*k1/2197.)-(7200*k2/2197.)+(7296*k3/2197.))
!         k5=(h/2.)*VPRIMA(v+(439*k1/216.)-(8.*k2)+(3680*k3/513.)-(845*k4/4104.))
!         k6=(h/2.)*VPRIMA(v-(8*k1/27.)+(2.*k2)-(3544*k3/2565.)+(1859*k4/4104.)-(11*k5/40.))
!         vfinal1=v+((25*k1/216.)+(1408*k3/2565.)+(2197*k4/4104.)-(k5/5.))
!
!         k1=(h/2.)*VPRIMA(vfinal1)   !HACE OTRO PASO MAS PARA HACER UN h ENTERO
!         k2=(h/2.)*VPRIMA(vfinal1+k1/4.)
!         k3=(h/2.)*VPRIMA(vfinal1+(3*k1/32.)+(9*k2/32.))
!         k4=(h/2.)*VPRIMA(vfinal1+(1932*k1/2197.)-(7200*k2/2197.)+(7296*k3/2197.))
!         k5=(h/2.)*VPRIMA(vfinal1+(439*k1/216.)-(8.*k2)+(3680*k3/513.)-(845*k4/4104.))
!         k6=(h/2.)*VPRIMA(vfinal1-(8*k1/27.)+(2.*k2)-(3544*k3/2565.)+(1859*k4/4104.)-(11*k5/40.))
!
!         vfinal1=vfinal1+((25*k1/216.)+(1408*k3/2565.)+(2197*k4/4104.)-(k5/5.)) !Aca define el 2do valor de "y" a comparar
!
!
!         k1=h*VPRIMA(v)                                               !Usa h
!         k2=h*VPRIMA(v+k1/4.)
!         k3=h*VPRIMA(v+(3*k1/32.)+(9*k2/32.))
!         k4=h*VPRIMA(v+(1932*k1/2197.)-(7200*k2/2197.)+(7296*k3/2197.))
!         k5=h*VPRIMA(v+(439*k1/216.)-(8.0*k2)+(3680*k3/513.)-(845*k4/4104.))
!         k6=h*VPRIMA(v-(8*k1/27.)+(2.*k2)-(3544*k3/2565.)+(1859*k4/4104.)-(11*k5/40.))
!
!         vfinal2=v+((25*k1/216.)+(1408*k3/2565.)+(2197*k4/4104.)-(k5/5.)) !Aca define el 2do valor de "y" a comparar
!
!         IF (ABS(vfinal1(2)-vfinal2(2))>errtol) THEN
!             h=(h/2.)
!             WRITE (*,*)'DIVIDE A h POR 2'
!         ELSE
!             IF (ABS(vfinal1(2)-vfinal2(2))<(errtol/4.)) THEN
!                 h=h*2
!                 WRITE (*,*)'MULTIPLICA A h POR 2'
!             ELSE
!                 v=vfinal2
!                 WRITE(2,formato)v
!                 iter=iter+1
!            END IF
!
!         END IF
!
!         !SEGUIMIENTO DE x y h
!         WRITE (*,*)
!         WRITE (*,'(A3,F6.3)')'X= ',v(1)
!         WRITE (*,'(A3,F10.8)')'h= ',h
!!         READ (*,*)                    !HABILITAR SI SE QUIERE VER CADA PASO
!         IF (h<hmin) THEN
!             hmin=h
!         END IF
!
!     END DO
!     WRITE (*,*)
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!     WRITE (*,'(A22,F10.8)')'Menor paso utilizado: ',h
!
!  CASE(3)                                                !Paso variable estrategia "error estimado por RKF"
!     DO WHILE ((v(1)-h)<xfinal)
!
!         k1=h*VPRIMA(v)
!         k2=h*VPRIMA(v+k1/4.)
!         k3=h*VPRIMA(v+(3*k1/32.)+(9*k2/32.))
!         k4=h*VPRIMA(v+(1932*k1/2197.)-(7200*k2/2197.)+(7296*k3/2197.))
!         k5=h*VPRIMA(v+(439*k1/216.)-(8.*k2)+(3680*k3/513.)-(845*k4/4104.))
!         k6=h*VPRIMA(v-(8*k1/27.)+(2.*k2)-(3544*k3/2565.)+(1859*k4/4104.)-(11*k5/40.))
!
!         v=v+((25*k1/216.)+(1408*k3/2565.)+(2197*k4/4104.)-(k5/5.))
!
!         eRKF=k1/360.-128*(k3/4275.)-2197*(k4/75240.)+k5/50.+2*(k6/55.)
!
!         WRITE(2,formato)v
!
!         IF (errtol>=MAXVAL(ABS(eRKF))) THEN
!             h=h*((errtol/MAXVAl(ABS(eRKF)))**0.2)
!             WRITE (*,*)'INCREMENTA h'
!         ELSE
!             h=h*((errtol/MAXVAl(ABS(eRKF)))**0.22)
!             WRITE (*,*)'REDUCE h'
!         END IF
!         WRITE (*,'(F26.24)')h
!         WRITE (*,*)
!!         READ (*,*)                    !HABILITAR SI SE QUIERE VER CADA PASO
!         iter=iter+1
!         IF (h<hmin) THEN
!             hmin=h
!         END IF
!
!     END DO
!     WRITE (*,'(A24,I4)')'Iteraciones realizadas: ',iter
!     WRITE (*,'(A22,F10.8)')'Menor paso utilizado: ',h
!
! END SELECT
!
!END SUBROUTINE
!
!
!!============================================================================================
!!============================================================================================
!!============================================================================================
!
!
!FUNCTION VPRIMA(v)
!REAL(8) v+1),VPRIMA+1)
!
! VPRIMA(1)=1.          !SIEMPRE VALE 1
!
! VPRIMA(2)=v(3)        !CADA "VPRIMA" VALE EL SIGUIENTE "v", EXCEPTO LA ECUACION DIFERENCIAL
!
! !VPRIMA(3)=v(4)
!
!!MODIFICAR ECUACION DIFERENCIAL ################################################################################
! VPRIMA(3)=(-2963.5-121.7*cos(5*v(1))-750.0*v(3)-29000.0*v(2))/296.35       !ACORDATE: v(1) es x. v(2) es y, v(3) es y' etc en las ec dif
!
!END FUNCTION
!
!
!!============================================================================================
!!============================================================================================
!!============================================================================================
!
!
!
END PROGRAM
