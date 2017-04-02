set terminal wxt
set output 'wavefunc0.png'
set autoscale 
set xtic auto                      
set ytic auto
set xlabel "r" 
set ylabel "w_f(r)"
set key right top	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2 pt 7 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 2 pt 11 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 13 ps 1
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 2 ps 1
set style line 5 lt 1 lc rgb "#191970" lw 2 pt 1 ps 1
set style line 6 lt 1 lc rgb "#006400" lw 2 pt 64 ps 1
set style line 7 lt 0 lc rgb "#f055f0" lw 2 pt 65 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 2 pt 73 ps 1

plot 'hydrogen_wavefunc.txt' u 1:2 t "wave function", \
'hydrogen_wavefunc.txt' u 1:3 t "derivative of wave function"