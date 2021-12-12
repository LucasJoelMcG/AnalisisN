 set   autoscale                        # escala los ejes automaticamente
 unset log                              # quita la escala logaritmica (si la hubiera)
 unset label                            # quita los titulos anteriores
 set xtic auto                          # establece automaticamente las divisiones del eje x
 set ytic auto                          # establece automaticamente las divisiones del eje y
 set grid
 set title "GRAFICO"
 set xlabel "x"
 set ylabel "y"
 set xrange[0:4500]
 plot  'nn0.txt' using 1:2 title 'N/N0' with lines
