set   autoscale                        # Escala los ejes automáticamente
unset log                              # Quita la escala logarítmica (si la hubiera)
unset label                            # Quita los títulos anteriores
set xtic auto                          # Establece automáticamente las divisiones del eje x
set ytic auto                          # Establece automáticamente las divisiones del eje y
set grid
set xlabel "x"
set ylabel "Valores de y"
# set zlabel "z"
								# Se puede también establecer el rango del eje x
							
								# Para más información: Terminal --> GNUPLOT --> test

set yrange [-10:15]	
set xrange [-15:15]

#linecolor 15
#dt dashtype
#lw linewidth
	
set title "Evolución de las Concentraciones"
plot  (y=0.5*exp(x/3)-sin(x)) title 'f(x)' with lines

	  
     
