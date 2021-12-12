program g2ej3b
    implicit none
    !cant escu = cant de y's , 
    !vii[0]=0 (viene de y(0)) y vii[1]=5 
    INTEGER, PARAMETER :: cant_ecu=2, max_itera=2500
    REAL(8), DIMENSION(0:cant_ecu) :: vii
    REAL(8) :: h1=0.001
    real(8) :: tolerancia=0.005
    INTEGER j
    vii(2)=0.0
    vii(1)=0.0 !esto es y
    vii(0)=0.0 !este es x
    !call RungeKuttaP(vii,cant_ecu,max_itera,h1,tolerancia)
    !call RungeKuttaFehlbergH(vii,cant_ecu,max_itera,h1,tolerancia)
    open (2,file='datos.txt',status='replace')
    WRITE (2, '(I5, 3F10.6)') 0, vii
    do j=1,max_itera
        vii=RungeKutta(vii,h1,cant_ecu)
        WRITE (2, '(I5, 3F16.10)') j, vii
    end do
    call execute_command_line ('gnuplot -p plot.plt')
    contains

    FUNCTION v_prima(v)
        !Definición de la Ecuacion Diferencial
        !v(0)=x v(1)=y
        REAL(8), DIMENSION(0:2) :: v, v_prima
        v_prima(0) = 1.0
        v_prima(1) = v(2)
        v_prima(2) = (24.0*sin(10.0*v(0))-50.0*v(1)-6.0*v(2))/0.5
    END FUNCTION v_prima

    SUBROUTINE RungeKuttaFehlberg(vi, cant_ec, max_iter, h)
        !Metodo de Runge-Kutta-Fehlberg (de 6to. orden)
        REAL(8), INTENT(IN), DIMENSION(0:cant_ec) :: vi
        REAL(8), INTENT(IN) :: h
        INTEGER, INTENT(IN) :: cant_ec, max_iter
        INTEGER iter
        REAL(8), DIMENSION(0:cant_ec) :: v, k1, k2, k3, k4, k5, k6, e
        open(2,file='datos.txt',status='replace')
        WRITE (*, '(I5, 3F10.6)') 0, vi, 0.0
        WRITE (2, '(I5, 3F10.6)') 0, vi, 0.0
        v = vi
        DO while (v(1)>4)
            k1 = h*v_prima(v)
            k2 = h*v_prima(v + k1/4.0)
            k3 = h*v_prima(v + (3.0*k1 + 9.0*k2)/32.0)
            k4 = h*v_prima(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
            k5 = h*v_prima(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0)
            k6 = h*v_prima(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
            v = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
            e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
            WRITE (*, '(I5, 4F16.10)') iter, v, e
            WRITE (2, '(I5, 4F16.10)') iter, v, e
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
        WRITE (2, '(I5, 3F10.6)') 0, vi, h
        DO iter = 1, max_iter
            cont=0
            auxv=v
            auxvp=v
            do while (((maxval(abs(auxv-auxvp)).LT.tol/2.0) .OR. (maxval(abs(auxv-auxvp)).GT.tol)) .AND. cont.LT.22)    
                auxv=v
                auxvp=v
                auxv=RungeKutta(auxv,h,cant_ec)
                do iter2=1,2
                    auxvp=RungeKutta(auxvp,h/2.0,cant_ec)
                end do
                if (maxval(abs(auxv-auxvp)).LT.tol/2.0) then
                    h=h*2.0
                else if (maxval(abs(auxv-auxvp)).GT.tol) then
                    h=h/2.0
                end if
                cont=cont+1
            end do
            WRITE (*, '(I5, F16.6)') cont, maxval(abs(auxv-auxvp))
            v=RungeKutta(v,h,cant_ec)
            !WRITE (*, '(I5, 3F16.6)') iter, v, h
            WRITE (2, '(I5, 3F16.6)') iter, v, h
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