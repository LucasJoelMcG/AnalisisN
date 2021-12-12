program g2ej3a
    implicit none
    !cant escu = cant de y's , 
    !vii[0]=0 (viene de y(0)) y vii[1]=5 
    INTEGER :: cant_ecu=1, max_itera=50
    REAL(8), DIMENSION(0:1) :: vii
    REAL(8) :: h1=0.01
    vii(1)=1.0 !esto es y(0)=1
    vii(0)=0.0 !este es el 0 de aca arriba
    call EulerSimple(vii,cant_ecu,max_itera,h1)
    call execute_command_line ('gnuplot -p plot.plt')
    contains

    FUNCTION v_prima(v)
        !Definición de la Ecuacion Diferencial
        !v(0)=x v(1)=y
        REAL(8), DIMENSION(0:1) :: v, v_prima
        v_prima(0) = 1.0
        v_prima(1) = -3*v(1)+exp(-v(0))
    END FUNCTION v_prima

    SUBROUTINE EulerSimple(vi, cant_ec, max_iter, h)
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8), INTENT(IN) :: h
        INTEGER iter
        REAL(8), DIMENSION(0:cant_ec) :: v
        open (2,file='datos.txt',status='replace')
        WRITE (*, '(I5, 2F10.6)') 0, vi
        WRITE (2, '(I5, 2F10.6)') 0, vi
        v = vi
        DO iter = 1, max_iter
            v = v + h*v_prima(v)
            WRITE (*, '(I5, 2F10.6)') iter, v
            WRITE (2, '(I5, 2F10.6)') iter, v
        END DO
        close(2,status='keep')
    END SUBROUTINE EulerSimple
end program g2ej3a