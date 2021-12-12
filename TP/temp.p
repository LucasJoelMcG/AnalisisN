set autoscale
unset label
unset log
set xtic auto
set ytic auto
set pm3d map
set nokey 
set xlabel 'Tiempo [s]'
set ylabel 'Altura de las latas [m]'
set title ' Variacion de la temperatura segun la altura en las latas [K]'
set yrange[0:0.124]
set xrange[0:144000]
splot 'datos.txt'
