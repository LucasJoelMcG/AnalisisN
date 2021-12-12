MODULE derivada
    ! Cantidad de Ecuaciones
    INTEGER, PARAMETER :: cant_ec = 1
    CONTAINS
    FUNCTION v_prima(v)
    ! Definición de la Ecuacion Diferencial
    REAL(8), DIMENSION(0:cant_ec) :: v, v_prima
    v_prima(0) = 1.0
    v_prima(1) = v(0)*v(1)**(1/3.)
    END FUNCTION v_prima
END MODULE