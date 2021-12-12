SUBROUTINE GaussJordan(A,X,B,n,A2) !A es la matris, B vec de t.ind. A2 Matriz ampliada
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
end subroutine