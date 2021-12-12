 set   autoscale
 unset log
 unset label
 set xtic auto
 set ytic auto
 set grid
 set title "Funcion aproximada"
 set xlabel "x"
 set ylabel "y"
  
 plot  'datos.dat' using 1:2 title 'datos' with points lw 3,\
       'puntos.dat' using 1:2 title 'aproximacion' with lines lt 3
