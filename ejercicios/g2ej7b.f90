program g2ej3b
    implicit none
    !cant escu = cant de y's , 
    !vii[0]=0 (viene de y(0)) y vii[1]=5 
    INTEGER :: cant_ecu=1, max_itera=40
    REAL(8), DIMENSION(0:1) :: vii
    REAL(8) :: h1=0.1
    vii(1)=1.0 !esto es y
    vii(0)=0.0 !este es x
    call RungeKuttaFehlberg(vii,cant_ecu,max_itera,h1)
    call execute_command_line ('gnuplot -p plot.plt')
    contains

    FUNCTION v_prima(v)
        !Definición de la Ecuacion Diferencial
        !v(0)=x v(1)=y
        REAL(8), DIMENSION(0:1) :: v, v_prima
        v_prima(0) = 1.0
        v_prima(1) = -v(0)*v(1)+(4*v(0)/v(1))
    END FUNCTION v_prima

    SUBROUTINE RungeKuttaFehlberg(vi, cant_ec, max_iter, h)
        !Metodo de Runge-Kutta-Fehlberg (de 6to. orden)
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8), INTENT(IN) :: h
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        INTEGER iter
        REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4, k5, k6, e
        open (2,file='datos.txt',status='replace')
        WRITE (*, '(I5, 2F10.6, 2F15.7)') 0, vi, 0.0, 0.0
        WRITE (2, '(I5, 2F10.6, 2F15.7)') 0, vi, 0.0, 0.0
        v = vi
        DO iter = 1, max_iter
            k1 = h*v_prima(v)
            k2 = h*v_prima(v + k1/4.0)
            k3 = h*v_prima(v + (3.0*k1 + 9.0*k2)/32.0)
            k4 = h*v_prima(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
            k5 = h*v_prima(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0-845.0*k4/4104.0)
            k6 = h*v_prima(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
            v = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
            e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
            WRITE (*, '(I5, 2F10.6, 2F22.19)') iter, v, e
            WRITE (2, '(I5, 2F10.6, 2F22.19)') iter, v, e
        END DO
        CLOSE(2,status='keep')
    END SUBROUTINE RungeKuttaFehlberg
end program g2ej3b