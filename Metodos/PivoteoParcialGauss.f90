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