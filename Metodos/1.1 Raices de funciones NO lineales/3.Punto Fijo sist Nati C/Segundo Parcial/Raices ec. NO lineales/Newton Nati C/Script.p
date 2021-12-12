set   autoscale                        # Escala los ejes automáticamente
unset log                              # Quita la escala logarítmica (si la hubiera)
unset label                            # Quita los títulos anteriores
set xtic auto                          # Establece automáticamente las divisiones del eje x
set ytic auto                          # Establece automáticamente las divisiones del eje y
set grid
set xlabel "x"
set ylabel "y"
# set zlabel "z"
								# Se puede también establecer el rango del eje x
							
								# Para más información: Terminal --> GNUPLOT --> test

set yrange [-1:1]	
set xrange [-5:5]

#linecolor 15
#dt dashtype
#lw linewidth
	
set title "Raices de f(x)"
plot  (y=x**4. + x**3. - 4.*x**2. -3.*x +3.) title 'f(x)' with lines,\
      (y=0) with lines

