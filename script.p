 set   autoscale                        # escala los ejes automaticamente
 unset log                              # quita la escala logaritmica (si la hubiera)
 unset label                            # quita los titulos anteriores
 set xtic auto                          # establece automaticamente las divisiones del eje x
 set ytic auto                          # establece automaticamente las divisiones del eje y
 set grid
 set title "GRAFICO"
 set xlabel "x"
 set ylabel "y"
 set xrange[-5:6]
 plot  (y=log(x**2+1)-exp(x/2)*cos(3.1416*x)) title "f(x)" with lines lt 3 ,\
 (y=0)
