program g2ej3b
    implicit none
    !cant escu = cant de y's , 
    !vii[0]=0 (viene de y(0)) y vii[1]=5 
    INTEGER, PARAMETER :: cant_ecu=2, max_itera=50
    REAL(8), DIMENSION(0:cant_ecu) :: vii
    REAL(8) :: h1=0.1
    real(8) :: tolerancia=0.0001
    INTEGER j
    !vii(6)=0.0
    !vii(5)=0.0
    !vii(4)=0.0
    !vii(3)=0.0
    vii(2)=0.0
    vii(1)=0.0 !esto es y
    vii(0)=0.0 !este es x
    call RungeKuttaFehlberg(vii,cant_ecu,max_itera,h1)
    !call RungeKuttaP(vii,cant_ecu,max_itera,h1,tolerancia)
    !call RungeKuttaFehlbergH(vii,cant_ecu,max_itera,h1,tolerancia)
    !open (2,file='datos.txt',status='replace')
    !WRITE (2, '(I5, 5F10.6)') 0, vii
    !do j=1,max_itera
    !    vii=RungeKutta(vii,h1,cant_ecu)
    !    WRITE (2, '(I5, 5F16.10)') j, vii
    !end do
    call execute_command_line ('gnuplot -p plot.plt')
    contains

    FUNCTION v_prima(v)
        !Definición de la Ecuacion Diferencial
        !v(0)=x v(1)=y
        REAL(8), DIMENSION(0:cant_ecu) :: v, v_prima
        v_prima(0) = 1.0
        v_prima(1) = v(2)
        v_prima(2) = (1-v(1)/0.001-v(2)*300)/200
        !v_prima(3) = (1-(v(2)*0)-(v(1)/0.001))/200.0
        !v_prima(4) = (1-v(2)*50-v(1)/0.001)/200
        !v_prima(5) = (1-v(2)*100-v(1)/0.001)/200
        !v_prima(6) = (1-v(2)*300-v(1)/0.001)/200
    END FUNCTION v_prima

    SUBROUTINE RungeKuttaFehlberg(vi, cant_ec, max_iter, h)
        !Metodo de Runge-Kutta-Fehlberg (de 6to. orden)
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8), INTENT(IN) :: h
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        INTEGER iter, i
        REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4, k5, k6, e
        open(2,file='datos3.txt',status='replace')
        WRITE (*, '(I5, 7F10.6)') 0, vi
        WRITE (2, '(I5, 7F10.6)') 0, vi
        v = vi
        DO i=1, max_iter
            k1 = h*v_prima(v)
            k2 = h*v_prima(v + k1/4.0)
            k3 = h*v_prima(v + (3.0*k1 + 9.0*k2)/32.0)
            k4 = h*v_prima(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
            k5 = h*v_prima(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0)
            k6 = h*v_prima(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
            v = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
            e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
            WRITE (*, '(I5, 7F16.10)') i, v
            WRITE (2, '(I5, 7F16.10)') i, v
        END DO
        close(2,status='keep')
    END SUBROUTINE RungeKuttaFehlberg

    SUBROUTINE RungeKuttaFehlbergH(vi, cant_ec, max_iter, h, tol)
        !Metodo de Runge-Kutta-Fehlberg (de 6to. orden)
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8) :: h
        real(8) :: tol, alfa
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        INTEGER iter
        REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4, k5, k6, e
        open (2,file='datos.txt',status='replace')
        WRITE (*, '(I5, 3F10.6, 2F15.7)') 0, vi, 0.0, h
        WRITE (2, '(I5, 3F10.6, 2F15.7)') 0, vi, 0.0, h
        v = vi
        alfa=0.2
        DO while (v(0)<2)
            k1 = h*v_prima(v)
            k2 = h*v_prima(v + k1/4.0)
            k3 = h*v_prima(v + (3.0*k1 + 9.0*k2)/32.0)
            k4 = h*v_prima(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
            k5 = h*v_prima(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0-845.0*k4/4104.0)
            k6 = h*v_prima(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
            v = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
            e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
            if (tol>=maxval(abs(e))) then
                alfa=0.2
            else
                alfa=0.22
            end if
            h=h*((abs(tol/maxval(abs(e))))**alfa)
            WRITE (*, '(I5, 3F10.6, 2F22.19)') iter, v, maxval(abs(e)), h
            WRITE (2, '(I5, 3F10.6, 2F22.19)') iter, v, maxval(abs(e)), h
        END DO
        CLOSE(2,status='keep')
    END SUBROUTINE RungeKuttaFehlbergH

    SUBROUTINE RungeKuttaP(vi, cant_ec, max_iter, h, tol)
        !Metodo de Runge-Kutta (de 4to. Orden)
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8) :: h
        real (8) :: tol
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        INTEGER iter, iter2, cont
        REAL(8), DIMENSION(0:cant_ec) :: v, auxv, auxvp, k1, k2, k3, k4
        open (2,file='datos.txt',status='replace')
        v=vi
        !WRITE (*, '(I5, 3F10.6)') 0, vi, h
        WRITE (2, '(I5, 6F10.6)') 0, vi, h
        DO iter = 1, max_iter
            cont=0
            auxv=v
            auxvp=v
            auxv=RungeKutta(auxv,h,cant_ec)
            do iter2=1,2
                auxvp=RungeKutta(auxvp,h/2.0,cant_ec)
            end do
            do while (((maxval(abs(auxv-auxvp)).LT.tol/2.0) .OR. (maxval(abs(auxv-auxvp)).GT.tol)) .AND. cont.LT.22) 
                if (maxval(abs(auxv-auxvp)).LT.tol/2.0) then
                    h=h*2.0
                else if (maxval(abs(auxv-auxvp)).GT.tol) then
                    h=h/2.0
                end if
            auxv=v
            auxvp=v
            auxv=RungeKutta(auxv,h,cant_ec)
            do iter2=1,2
                auxvp=RungeKutta(auxvp,h/2.0,cant_ec)
            end do
                cont=cont+1
            end do
            WRITE (*, '(I5, F16.6)') cont, maxval(abs(auxv-auxvp))
            v=RungeKutta(v,h,cant_ec)
            !WRITE (*, '(I5, 3F16.6)') iter, v, h
            WRITE (2, '(I5, 6F16.6)') iter, v, h
        END DO
        CLOSE(2,status='keep')
    END SUBROUTINE RungeKuttaP

    function RungeKutta(v,h,cant_ec)
        !Metodo de Runge-Kutta (de 4to. Orden)
        REAL(8), INTENT(IN) :: h
        INTEGER, INTENT(IN) :: cant_ec
        REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4, RungeKutta
        k1 = h*v_prima(v)
        k2 = h*v_prima(v+k1/2.0)
        k3 = h*v_prima(v+k2/2.0)
        k4 = h*v_prima(v+k3)
        v = v + (k1 + 2.0*k2 + 2.0*k3 +k4)/6.0
        RungeKutta=v
    END function
end program g2ej3b