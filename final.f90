!############### MODIFICACIONES DEL PROGRAMA #################
!   -Cambiar los valores largo, deltax, cp, p, cond y tmax
!   -Ingresar condiciones iniciales y de contorno (pueden ser puntos o funciones)

PROGRAM tp_final_intento_2
IMPLICIT NONE

integer(4) orden,k,c,j,i,conjsol
real(8) deltat,r,orden1
real(8),ALLOCATABLE :: M(:,:),A(:,:)
real(8),ALLOCATABLE :: L(:),D(:),U(:),B(:),X(:),vecsol(:)
REAL(8),PARAMETER :: largo=0.124              !LARGO DE LA BARRA
REAL(8),PARAMETER :: deltax=0.002             !DELTA X (CANTIDAD INTERVALOS ENTERO) EN GRAL LO INVENTO YO
REAL(8),PARAMETER :: Cp=3883.                 !CALOR ESPECIFICO
REAL(8),PARAMETER :: p=1130.                  !DENSIDAD
INTEGER,PARAMETER :: tmax=100000              !CANTIDAD DE ITERACIONES TEMPORALES
REAL(8),PARAMETER :: pi=3.1416
REAL(8),PARAMETER :: cond=0.4324              !CONDUCTIVIDAD TERMINA
REAL(8),PARAMETER :: tol=0.01
REAL(8),PARAMETER :: h=25.0                   !coef de convecc del aire
!TMAX ESTÁ RELACIONADO CON EL DELTAX. LO SETEO SEGÚN EN QUÉ TIEMPO QUIERO LA SOLUCIÓN
orden1=largo/deltax 
orden1=orden1+0.5
orden=int(orden1)
WRITE(*,*)'La cantidad de divisiones es: ',orden
        deltat=100.                    !DELTA T
        r=(deltat*cond)/(p*Cp*(deltax**2.))
        WRITE(*,*)'El radio es ',r

ALLOCATE(M(0:tmax,0:orden))
ALLOCATE(A(orden-1,orden-1))   
ALLOCATE(L(orden-1),B(orden-1),vecsol(orden-1),U(orden-1),X(orden-1))
ALLOCATE(D(orden))
A(:,:)=0

DO i=1,orden-1
    DO j=1,orden-1
        IF (j==i) THEN
            A(i,j)=2.+2.*r
        ELSE IF ((j==i+1).or.(j==i-1)) THEN
            A(i,j)=-r
        END IF
    END DO
END DO

M=matriz(orden,tmax)

OPEN(2,FILE='crank.dat',STATUS='REPLACE')
OPEN(3,FILE='prueba.dat',STATUS='REPLACE')
!DO k=1,tmax
k=1
WRITE(3,*) M(k,:)
do while (118.0-minval(M(k-1,:))>=tol .AND. k<=tmax)
  if (M(k,orden)<118.) then                                  
    M(k,orden) = (M(k-1,orden-1) + h * deltax *118./cond) / (1 + (h * deltax)/cond)     
  else 
    M(k,orden) = (M(k-1,orden-1) - h * deltax *118./cond) / (1 - (h * deltax)/cond)
  end if
    B=vectorB(M,orden,tmax,k,r)
    conjsol=1 
    CALL thomas(A,orden-1,conjsol,B,vecsol)
    DO c=1,orden-1
        M(k,c)=vecsol(c)
    END DO
    WRITE(3,*) M(k,:)
    k=k+1
END DO
!CALL titulo(orden)

CALL escrmatriz(tmax+1,orden+1,M,deltat,k)
CLOSE(2,STATUS='KEEP')   
print*,'LA SOLUCION ESTA EN EL ARCHIVO'
call system ("gnuplot.exe -persist crank.p")

CONTAINS
 !####################################################################
FUNCTION matriz(orden,tmax)

 integer(4) i,orden,tmax
 real(8), DIMENSION(0:tmax,0:orden) :: M, matriz
 real(8) x
 M(:,:)=0

           
 !CONDICIONES INICIALES
 !#### O SETEO UNA FUNCIÓN O SETEO VALORES PUNTUALES PARA LOS NODOS INTERNOS#####
M=20

Do j=0,orden
    M(0,j) = 20
end do

!~ Do j=0,100
!~    M(0,j)=0.   !MULTIPLICAR J POR EL VALOR DE DELTAX
!~ End do
!~ do j=100,orden
!~    M(0,j)=1.
!~ end do
 
 !CONDICIONES DE CONTORNO (ACÁ REESCRIBE LAS COND. INIC. DE LS BORDES SETEADAS ARRIBA
 Do i=0,tmax                            
   M(i,0)=118.
   !M(0,orden)=20.
 end do
 
 matriz=M
 
 END FUNCTION matriz
 !####################################################################
 FUNCTION vectorB(M,orden,tmax,t,r1)
 
 integer(4) j,orden,t,tmax
 real(8) r1
 real(8), DIMENSION(0:tmax,0:orden) :: M
 real(8), DIMENSION(orden-1) :: vectorB,B
 
   Do j=2,orden-2
     B(j)=r1*M(t-1,j-1)+(2-2*r1)*M(t-1,j)+r1*M(t-1,j+1)
   end do
     B(1)=r1*M(t,0)+r1*M(t-1,0)+(2-2*r1)*M(t-1,1)+r1*M(t-1,2)
     B(orden-1)=r1*M(t-1,orden-2)+r1*M(t-1,orden)+(2-2*r1)*M(t-1,orden-1)+r1*M(t,orden)
   vectorB=B
 
 END FUNCTION vectorB    
  !####################################################################
 SUBROUTINE Thomas(A,orden,conjsol,B,vecsol)

 integer(4) orden,i,j,conjsol
 real(8), DIMENSION (orden) :: D
 real(8), DIMENSION (orden-1) :: L,U
 real(8), DIMENSION (orden,conjsol) :: B,vecsol
 real(8), DIMENSION (orden,orden+conjsol) :: A
  
 Do i=1,orden
   D(i)=A(i,i)
 end do

 Do i=2,orden
   L(i-1)=A(i,i-1)
 end do  

 Do i=1,orden-1
   U(i)=A(i,i+1)
 end do 

 Do i=1,orden-1
   If (d(i).ne.0) then
     u(i)=u(i)/d(i)
     B(i,:)=B(i,:)/d(i)
     d(i)=1
     d(i+1)=d(i+1)-l(i)*u(i)
     do j=1,conjsol
     B(i+1,j)=B(i+1,j)-l(i)*B(i,j)
     end do
     l(i)=0
   end if  
 end do
 Do j=1,conjsol
   vecsol(orden,j)=B(orden,j)/D(orden)
   i=orden-1
   Do While (i.ge.1)
     vecsol(i,j)=B(i,j)-U(i)*vecsol(i+1,j)
     i=i-1
   end do
 end do
 end subroutine thomas    
 !####################################################################
Subroutine titulo(nodo)
integer(4) nodo,j
Write(2,'(A10)', advance='no') 'tiempo='
Do j=1, nodo  
  Write(2,'(A10)', advance='no') 'nodo=', j 
end do 
end subroutine
 !####################################################################
subroutine escrmatriz(filas,columnas,mat,deltat,k)
integer(4) i, j, filas, columnas,k 
real(8), dimension(filas,columnas) :: mat
real(8) t,deltat,x

t=0
x=0
Do i=1,k
  do j=1,columnas
    write(2,'(F15.4)', advance='no') t
    write(2, '( F15.4  F25.12)') x, mat(i,j)
    x= x+deltax
  end do
  x=0
  write(2,*)
  t=t+deltat
end do
end subroutine escrmatriz
 !####################################################################
END PROGRAM
