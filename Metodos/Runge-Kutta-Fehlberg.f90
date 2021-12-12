SUBROUTINE RungeKuttaFehlberg(vi, cant_ec, max_iter, h)
    !Metodo de Runge-Kutta-Fehlberg (de 6to. orden)
    REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
    REAL(8), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: cant_ec, max_iter
    INTEGER iter
    REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4, k5, k6, e
    WRITE (*, '(I5, 3F10.6)') 0, vi, 0
    v = vi
    DO iter = 1, max_iter
        k1 = h*v_prima(v)
        k2 = h*v_prima(v + k1/4.0)
        k3 = h*v_prima(v + (3.0*k1 + 9.0*k2)/32.0)
        k4 = h*v_prima(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
        k5 = h*v_prima(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 – 845.0*k4/4104.0)
        k6 = h*v_prima(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 – 11.0*k5/40.0)
        v = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
        e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
        WRITE (*, '(I5, 3F10.6)') iter, v, e
    END DO
END SUBROUTINE RungeKuttaFehlberg