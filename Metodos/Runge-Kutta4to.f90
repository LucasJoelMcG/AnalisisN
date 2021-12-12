SUBROUTINE RungeKutta(vi, cant_ec, max_iter, h)
    !Metodo de Runge-Kutta (de 4to. Orden)
    REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
    REAL(8), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: cant_ec, max_iter
    INTEGER iter
    REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4
    WRITE (*, '(I5, 2F10.6)') 0, vi
    v = vi
    DO iter = 1, max_iter
        k1 = h*v_prima(v)
        k2 = h*v_prima(v+k1/2.0)
        k3 = h*v_prima(v+k2/2.0)
        k4 = h*v_prima(v+k3)
        v = v + (k1 + 2.0*k2 + 2.0*k3 +k4)/6.0
    WRITE (*, '(I5, 2F10.6)') iter, v
    END DO
END SUBROUTINE RungeKutta