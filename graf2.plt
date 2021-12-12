# Ejemplo de script Gnuplot, para graficar los datos 
unset   angles                      # escala los ejes automaticamente
unset log                              # quita la escala logaritmica (si la hubiera)
unset label                            # quita los titulos anteriores
set xtic auto                          # establece automaticamente las divisiones del eje x
set ytic auto                          # establece automaticamente las divisiones del eje y
set grid
set title "Función"
set xlabel "x"
set ylabel "y"

plot "graficoptofijo.dat" using 1:2 title 'Función f(x)' with lines,\
     "graficoptofijo.dat" using 1:3 title 'recta y = x' with lines,\
     "graficoptofijo.dat" using 1:4 title 'Funcion g(x)' with lines,\
     "graficoptofijo.dat" using 1:5 title 'Funcion g(x)derivada ' with lines