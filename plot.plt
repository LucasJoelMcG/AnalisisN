set  autoscale                  #escala los ejes automaticamente
unset log                       #quita la escala logaritmica
unset label                     #quita los titulos anteriores
set xtic auto                   #establece las divisiones del eje x
set ytic auto                   #establece las divisines del eje y

set title "Concentraciones - RK 4to"
set xlabel "x - V(0)"
set ylabel "f(x)"

plot "datos.txt" using 1:2 title "C1" with lines ,\
"datos.txt" using 1:3 title "C2" with lines ,\
"datos.txt" using 1:4 title "C3" with lines ,\
"datos.txt" using 1:5 title "C4" with lines ,\