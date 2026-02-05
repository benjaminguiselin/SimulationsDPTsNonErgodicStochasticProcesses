set term qt
set xlabel 'q'
set ylabel '{/Symbol f}(q)'
set g
set key top center
set samples 10000
p './Data/ldf_end_position.txt' w p title 'simu', \
[-2:2] abs(x) w l lc rgb "#000000" title 'asymptotic', \
[2:] 1 + 0.25 * x ** 2 w l lc rgb "#000000" notitle, \
[:-2] 1 + 0.25 * x ** 2 w l lc rgb "#000000" notitle, \
'./Data/ldf_end_position.txt' every 12::1 w p title 'simu (short)'
pause -1

set term cairolatex pdf font ',8' linewidth 2 size 8.5cm, 5.3cm standalone header "\\usepackage{amsmath,bm}"
set output 'ldf_end_position_mortal.tex'
set xlabel '$q$'
set ylabel '$\phi(q)$' offset char -0.1, 0
set xr [-4:4]
unset key
unset g
set samples 10000
set ytics 1
set mytics 2
set xtics 1
set mxtics 2
set arrow 1 from 2, 0 to 2, 2 lw 3 dt 2 nohead
set arrow 2 from -2, 0 to -2, 2 lw 3 dt 2 nohead
set label 1 at 2.12, 0.3 '$q_c$'
set label 2 at -2.55, 0.3 '$-q_c$'
p './Data/ldf_end_position.txt' every 12::1 w p pt 1 ps 0.2 lc rgb "#ffffff", \
[-2:2] abs(x) w l lw 3 lc rgb "#1D60FF", \
[2:] 1 + 0.25 * x ** 2 w l lw 3 lc rgb "#FF0380", \
[:-2] 1 + 0.25 * x ** 2 w l lw 3 lc rgb "#FF0380", \
'./Data/ldf_end_position.txt' every 12::1 w p pt 7 ps 0.3 lc rgb "#000000"