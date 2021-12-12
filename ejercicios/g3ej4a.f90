program g3ej2
    implicit none
    integer ii
    integer, parameter :: ord=4 , indep=1    !n:orden de la matriz, indep=orden de vector terminos independientes
    real(8) Mi(ord,ord), bi(ord), mat_ampl(ord,ord+indep) !para gauss,regresiva y crout
    !real(8) vecU(n),vecD(n),vecL(n),tIndep(n,indep)!para thomas
    Mi(1,1)=4.5556;    Mi(1,2)=0.0345;     Mi(1,3)=12.0009;     Mi(1,4)=-1.8887;     bi(1)=17.102
    Mi(2,1)=0.567;    Mi(2,2)=1.3456;     Mi(2,3)=1234.9;     Mi(2,4)=34.5;     bi(2)=-6.1593
    Mi(3,1)=4.56;    Mi(3,2)=0.1;    Mi(3,3)=12.0;     Mi(3,4)=-2.0;     bi(3)=3.0004
    Mi(4,1)=38.01;    Mi(4,2)=-3.67;    Mi(4,3)=2.10;     Mi(4,4)=-3.0;     bi(4)=0.0
    !do ii=1,n
    !    vecU(ii)=1
    !    vecL(ii)=1
    !    vecD(ii)=4
    !    read *, tIndep(ii,indep) 
    !end do
    !vecU(n)=0
    !vecL(1)=0
    !call matrizampliada(Mi,bi,mat_ampl)
    !call gauss(mat_ampl)
    !call GaussJordan(mat_ampl,n)
    !call regresiva(mat_ampl) 
    !call CroutReduc(Mi,n)
    !call CroutSol(Mi,n,bi)
    !call Thomas(vecU,vecD,vecL,tIndep)
    call Inversa(Mi,ord)
    

contains

subroutine pivoteo(M,col)
    real(8) M(ord,ord+indep)
    integer i, col,fila
    real(8) maximo,aux !analiza tambien el primer valor
    write(*,*)
    maximo=0.0
    do i=col, ord
      if (abs(M(col,i))>maximo) then
        maximo=abs(M(col,i))
        fila=i 
      end if
    end do
    do i=col, ord+1
      aux=M(col,i)
      M(col,i)=M(fila,i)
      M(fila,i)=aux 
    end do
    write(*, '(5f10.3)') M 
    write(*,*)
end subroutine
    
subroutine gauss(M)
    real(8) M(ord,ord+1)
    integer t,fila
    Do t=1, ord !t=columna
    call pivoteo(M,t)
      do fila=t+1, ord
        M(fila, t+1:) = M(fila, t+1:) - (M(t, t+1:)*M(fila, t)/M(t, t))
        M(fila,t)=0.0
      end do
    end do
end subroutine
    
Subroutine regresiva(M)
    real(8) M(ord,ord+1),x(ord)
    integer i,j,aux
    real(8) suma
    aux=ord
    x(ord)=M(ord,ord+1)/M(ord,ord)
    do i=aux-1,1,-1
      suma=0.0
      do j=i+1, ord
        suma=suma + (M(i,j)*x(j))
      end do
      x(i)=(M(i,ord+1) - suma)/M(i,i)
    end do
    print *, 'El vector solución es: '
    write (*,'(4f15.5)') X !cambiar segun el orden de la matriz
end subroutine

subroutine MostrarMat(matriz,m)
    integer i,j,m
    real(8),DIMENSION(ord,m) :: matriz
    do i=1,m
        WRITE(*,*) matriz(i,:)
    end do
end subroutine

subroutine matrizampliada(M,b,mat_amp)
    real(8) M(ord,ord), b(ord)
    real(8) mat_amp(ord,ord+indep)
    integer i,j
    mat_amp=0.0
    !print*, 'la matriz ampliada es:'
    do i=1, ord
      do j=1, ord
        mat_amp(i,j)=M(i,j)
      end do
    end do
    do i=1,ord
      mat_amp(i,ord+1)=b(i)
    end do
    !do i=1, n
    !  write(*, '(5f10.3)', advance='no') mat_amp(i,:)
    !  write(*,*)
    !end do
end subroutine

SUBROUTINE GaussJordan(mat_amp,orden)! Metodo de Gauss-Jordan
    !REAL(8), DIMENSION(:,:), INTENT(IN) :: matriz, term_indep
    INTEGER t, fila, orden
    REAL(8), DIMENSION(orden,orden+1) :: mat_amp
    !orden = SIZE(matriz, DIM=1)
    !CALL creaMatrizAmpliada(matriz, term_indep, mat_amp)
    DO t=1, orden
        DO fila=1, t-1
            mat_amp(fila,t+1:) = mat_amp(fila,t+1:) - mat_amp(t,t+1:)*mat_amp(fila,t) / mat_amp(t,t) 
            mat_amp(fila,t) = 0.0
        END DO
        DO fila=t+1, orden
            mat_amp(fila,t+1:)=mat_amp(fila,t+1:)-mat_amp(t,t+1:)*mat_amp(fila,t) / mat_amp(t,t)
            mat_amp(fila,t) = 0.0
        END DO
    END DO
    !CALL imprimeMatriz(mat_amp(:,orden+1:))
    !DEALLOCATE(mat_amp) 
END SUBROUTINE

SUBROUTINE creaMatrizIdentidad(N,MatID)
    REAL(8), INTENT(INOUT), ALLOCATABLE :: MatID(:,:)
    INTEGER N,i
    ALLOCATE (MatID(N,N))
    MatID=0.0
    DO i=1,N
        MatID(i,i)=1.0		
    END DO
    !CALL muestramat(MatID,N,N)
END SUBROUTINE creaMatrizIdentidad

SUBROUTINE Pivoteo_Inversa(A, MatID, k, N, pivote)
    !Variables
    REAL(8), ALLOCATABLE :: A(:,:), MatID(:,:)
    REAL(8), ALLOCATABLE :: aux(:), aux2(:)
    REAL(8) pivote
    INTEGER N, i, k, filap
    ALLOCATE(aux(N))
    ALLOCATE(aux2(N))
    pivote = A(k,k)
    filap = k
    !write(*,'(F10.4,I5)')pivote,filap
    DO i=k+1, N
        IF (abs(A(i,k)) > abs(pivote)) THEN
            pivote = A(i,k)
            filap = i
        END IF
    END DO
    aux(:) = A(k,:)
    aux2(:) = MatID(k,:)
    A(k,:) = A(filap,:)
    MatID(k,:) = MatID(filap,:)
    A(filap,:) = aux(:)
    MatID(filap,:) = aux2(:)
    DEALLOCATE(aux)
    DEALLOCATE(aux2)
END SUBROUTINE Pivoteo_Inversa

SUBROUTINE Gauss_Jordan_Inversa(Mat_Inv,MatID,N)
    !Variable
    REAL(8), INTENT(INOUT), ALLOCATABLE :: MatID(:,:),Mat_Inv(:,:)
    !REAL(8), INTENT(IN), ALLOCATABLE :: A(:,:)
    REAL(8) pivote
    INTEGER N, i, j, k
    !Triangulación Inferior
    DO k=1, N
        !write(*,*) 'prueba'
        CALL Pivoteo_Inversa(Mat_Inv, MatID, k, N, pivote)
        Mat_Inv(k,:)= Mat_Inv(k,:) / pivote
        MatID(k,:) = MatID(k,:) / pivote
        DO i = k+1, N
            DO j=k+1, N
                Mat_Inv(i,j) = Mat_Inv(i,j) - Mat_Inv(k,j) * Mat_Inv(i,k)
            END DO
            MatID(i,:) = MatID(i,:) - MatID(k,:) * Mat_Inv(i,k)
            Mat_Inv(i,k) = 0.
        END DO
    END DO
    !Triangulación Superior
    pivote = 1.
    Mat_Inv(N,:) = Mat_Inv(N,:) / pivote
    MatID(N,:) = MatID(N,:) / pivote
    DO k=N, 1, -1
        DO i=k-1, 1, -1
            DO j=i+1, k-1
                Mat_Inv(i,j) = Mat_Inv(i,j) - Mat_Inv(k,j) * Mat_Inv(i,k)
            END DO
            MatID(i,:) = MatID(i,:) - MatID(k,:) * Mat_Inv(i,k)
            Mat_Inv(i,k) = 0.
        END DO
    END DO
    Mat_Inv = MatID !A la matriz A le asignamos su inversa
END SUBROUTINE Gauss_Jordan_Inversa

SUBROUTINE Inversa(Mat_Inv,N)
!Calcula la inversa de una matriz, por Gauss-Jordan
    REAL(8), INTENT(INOUT) :: Mat_Inv(:,:)
    REAL(8), ALLOCATABLE :: MatID(:,:)
    INTEGER N
    !creo la matriz identidad
    N = ord
    CALL creaMatrizIdentidad(N, MatID)
    CALL Gauss_Jordan_Inversa(Mat_Inv, MatID, N)
    !DEALLOCATE(MatID)
END SUBROUTINE

SUBROUTINE Thomas(u_orig, d_orig, l_orig, term_indep)! Metodo de Thomas para matrices Tri-Diagonales
    REAL(8), DIMENSION(:,:), INTENT(IN) :: term_indep
    REAL(8), DIMENSION(:), INTENT(IN) :: u_orig, d_orig, l_orig
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: b
    REAL(8), DIMENSION(:), ALLOCATABLE :: u, d, l  
    INTEGER orden, cant_vec, i
    orden = SIZE(term_indep, DIM=1)
    cant_vec = SIZE(term_indep, DIM=2)! Realiza copias del original
    ALLOCATE(u(orden), d(orden), l(orden))
    u = u_orig
    d = d_orig
    l = l_orig
    ALLOCATE(b(orden, cant_vec))
    b = term_indep
    DO i=1, orden-1
        u(i) = u(i) / d(i)
        b(i,:) = b(i,:) / d(i)
        d(i) = 1.0
        d(i+1) = d(i+1) - l(i+1)*u(i)
        b(i+1,:) = b(i+1,:) - l(i+1)*b(i,:)
        l(i+1) = 0.0
    END DO
    ! Obtencion de la solucion por Sustitucion Inversa
    b(orden,:) = b(orden,:)/d(orden)
    DO i=orden-1, 1, -1
        b(i,:) = b(i,:) - u(i)*b(i+1,:) / d(i)
    END DO
    do i=1,orden
        write (*,*) b(i,1)
    end do
    DEALLOCATE(u, d, l, b)
END SUBROUTINE

SUBROUTINE Residuos(A2,vecsol,B,orden,conjsol)
    real(8) A2(:,:), B(:,:),vecsol(:,:)
    integer(4) conjsol,j,orden
    real(8), dimension(orden,conjsol) :: residuo
    do j=1,conjsol
        residuo=(matmul(A2,vecsol)-B)
    end do
    print*,'El vector residuo es'
    write(*,*) residuo
end subroutine Residuos

subroutine CroutReduc(A,orden)
    REAL(8) A(ord,ord)
    INTEGER fila, col, i, k, orden
    A(1,2:)=A(1,2:) / A(1,1)
    DO i=2, orden ! Calcula la columna i
        DO fila=i, orden
            DO k=1,i-1
                A(fila,i)=A(fila,i) - A(fila,k)*A(k,i)
            END DO 
        END DO! Calcula fila i
        DO col=i+1, orden
            DO k=1,i-1
                A(i, col)=A(i,col) - A(i,k)*A(k,col)
            END DO
            A(i,col) = A(i,col) / A(i,i) 
        END DO
    END DO
end subroutine

subroutine CroutSol(A,orden,b)
    REAL(8) A(orden,orden), c(orden),b(orden),x(orden)
    INTEGER fila, col, i, k, orden! Calcula el vector c
    c(1) = b(1) / A(1,1)
    DO i=2, orden
        c(i) = b(i)
        DO k=1, i-1
            c(i) = c(i) - A(i,k)*c(k)
        END DO
        c(i) = c(i) / A(i,i)
    END DO! Calcula la Solución
    x(orden) = c(orden)
    DO i=(orden-1), 1, -1
        x(i) = c(i)
        DO k=i+1,orden
            x(i) = x(i) - A(i,k)*x(k)
        END DO
    END DO
    write(*,*) x
end subroutine

END PROGRAM