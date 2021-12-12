SUBROUTINE EulerModificado(vi, cant_ec, max_iter, h)
    INTEGER, INTENT(IN) :: cant_ec, max_iter
    REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
    REAL(8), INTENT(IN) :: h
    INTEGER iter
    REAL(8), DIMENSION(0:cant_ec) :: v, vp
    WRITE (*, '(I5, 2F10.6)') 0, vi
    v = vi
    DO iter = 1, max_iter
        vp = v_prima(v)
        v = v + h*(vp + v_prima(v + h*vp))/2.0
    WRITE (*, '(I5, 2F10.6)') iter, v
    END DO
END SUBROUTINE EulerModificado
