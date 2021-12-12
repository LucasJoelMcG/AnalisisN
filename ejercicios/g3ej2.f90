program g3ej2
    implicit none
    integer, parameter :: n=4 , indep=1    !n:orden de la matriz, indep=orden de vector terminos independientes
    real(8) Mi(n,n), bi(n), mat_ampl(n,n+indep)

    Mi(1,1)=2.1756;    Mi(1,2)=4.0231;     Mi(1,3)=-2.1732;     Mi(1,4)=5.1967;     bi(1)=17.102
    Mi(2,1)=-4.0231;    Mi(2,2)=6.0;     Mi(2,3)=0.0;     Mi(2,4)=1.1973;     bi(2)=-6.1593
    Mi(3,1)=-1.0;    Mi(3,2)=-5.2107;    Mi(3,3)=1.1111;     Mi(3,4)=0.0;     bi(3)=3.0004
    Mi(4,1)=6.0235;    Mi(4,2)=7.0;    Mi(4,3)=0.0;     Mi(4,4)=-4.1561;     bi(4)=0.0

    call matrizampliada(Mi,bi,mat_ampl)
    !call gauss(mat_ampl)
    !call GaussJordan(mat_ampl,n)
    !call regresiva(mat_ampl) 
    call CroutReduc(Mi,n)
    call CroutSol(Mi,n,bi)

contains

subroutine pivoteo(M,col)
    real(8) M(n,n+indep)
    integer i, col,fila
    real(8) maximo,aux !analiza tambien el primer valor
    write(*,*)
    maximo=0.0
    do i=col, n
      if (abs(M(col,i))>maximo) then
        maximo=abs(M(col,i))
        fila=i 
      end if
    end do
    do i=col, n+1
      aux=M(col,i)
      M(col,i)=M(fila,i)
      M(fila,i)=aux 
    end do
    write(*, '(5f10.3)') M 
    write(*,*)
end subroutine
    
subroutine gauss(M)
    real(8) M(n,n+1)
    integer t,fila
    Do t=1, n !t=columna
    call pivoteo(M,t)
      do fila=t+1, n
        M(fila, t+1:) = M(fila, t+1:) - (M(t, t+1:)*M(fila, t)/M(t, t))
        M(fila,t)=0.0
      end do
    end do
end subroutine
    
Subroutine regresiva(M)
    real(8) M(n,n+1),x(n)
    integer i,j,aux
    real(8) suma
    aux=n
    x(n)=M(n,n+1)/M(n,n)
    do i=aux-1,1,-1
      suma=0.0
      do j=i+1, n
        suma=suma + (M(i,j)*x(j))
      end do
      x(i)=(M(i,n+1) - suma)/M(i,i)
    end do
    print *, 'El vector solución es: '
    write (*,'(4f15.5)') X !cambiar segun el orden de la matriz
end subroutine

subroutine MostrarMat(matriz,m)
    integer i,j,m
    real(8),DIMENSION(n,m) :: matriz
    do i=1,m
        WRITE(*,*) matriz(i,:)
    end do
end subroutine

subroutine matrizampliada(M,b,mat_amp)
    real(8) M(n,n), b(n)
    real(8) mat_amp(n,n+indep)
    integer i,j
    mat_amp=0.0
    !print*, 'la matriz ampliada es:'
    do i=1, n
      do j=1, n
        mat_amp(i,j)=M(i,j)
      end do
    end do
    do i=1,n
      mat_amp(i,n+1)=b(i)
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

!SUBROUTINE Inversa(matriz)! Calcula la inversa de una matriz, por Gauss-Jordan
!    REAL(8), DIMENSION(:,:), INTENT(IN) :: matriz
!    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: matrizIdentidad
!    INTEGER orden
!    orden = SIZE(matriz, DIM=1)
!    CALL creaMatrizIdentidad(orden, matrizIdentidad)
!    CALL GaussJordan(matriz, matrizIdentidad)
!    DEALLOCATE(matrizIdentidad)
!END SUBROUTINE

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
    !CALL imprimeMatriz(b)
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
    write(3,*) residuo
end subroutine Residuos

subroutine CroutReduc(A,orden)
    REAL(8) A(n,n)
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