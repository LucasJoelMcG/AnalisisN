program g2ej1
    implicit none
    !y'=x.y^1/3
    REAL(8), DIMENSION(0:cant_ecu) :: vii
    INTEGER :: cant_ecu, max_itera
    REAL(8) :: h1

    contains

    FUNCTION v_prima(v)
        ! Definición de la Ecuacion Diferencial
        REAL(8), DIMENSION(0:cant_ec) :: v, v_prima
        v_prima(0) = 1.0
        v_prima(1) = v(0)*v(1)**(1/3.)
    END FUNCTION v_prima

    SUBROUTINE EulerSimple(vi, cant_ec, max_iter, h)
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8), INTENT(IN) :: h
        INTEGER iter
        REAL(8), DIMENSION(0:cant_ec) :: v
        WRITE (*, '(I5, 2F10.6)') 0, vi
        v = vi
        DO iter = 1, max_iter
            v = v + h*v_prima(v)
            WRITE (*, '(I5, 2F10.6)') iter, v
        END DO
    END SUBROUTINE EulerSimple

end program g2ej1