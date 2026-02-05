set term qt
set xlabel 'q'
set ylabel '{/Symbol f}(q)'
set g
set key top center
set samples 10000
p './Data/ldf_area.txt' w p title 'simu', \
[-2./3:2./3] 2. * sqrt(2. * abs(x) / 3.) w l lc rgb "#000000" title 'asymptotic', \
[2./3:] 1 + 0.75 * x ** 2 w l lc rgb "#000000" notitle, \
[:-2./3] 1 + 0.75 * x ** 2 w l lc rgb "#000000" notitle, \
'./Data/ldf_area.txt' every 40::32 w p title 'simu (short)', \
'./Analytics/ldf_area_theory.txt' u 1:($2+0.06) w l lw 2 dt 2 lc rgb "#28b463"
pause -1

set term cairolatex pdf font ',8' linewidth 2 size 8.5cm, 5.3cm standalone header "\\usepackage{amsmath,bm}"
set output 'ldf_area_mortal.tex'
set xlabel '$q$'
set ylabel '$\phi(q)$' offset char -0.1, 0
set xr [-1.5:1.5]
set yr [0:3]
unset key
unset g
set samples 10000
set ytics 1
set mytics 2
set xtics 0.5
set mxtics 2
set arrow 1 from 2./3, 0 to 2./3, 4./3 lw 3 dt 2 nohead
set arrow 2 from -2./3, 0 to -2./3, 4./3 lw 3 dt 2 nohead
set label 1 at 0.72, 0.2 '$q_c$'
set label 2 at -0.88, 0.2 '$-q_c$'
p './Data/ldf_area.txt' every 40::32 w p pt 1 ps 0.2 lc rgb "#ffffff", \
[-2./3:2./3] 2. * sqrt(2. * abs(x) / 3.) w l lw 3 lc rgb "#1D60FF", \
[2./3:] 1 + 0.75 * x ** 2 w l lw 3 lc rgb "#FF0380", \
[:-2./3] 1 + 0.75 * x ** 2 w l lw 3 lc rgb "#FF0380", \
'./Data/ldf_area.txt' every 40::32 w p pt 7 ps 0.3 lc rgb "#000000", \
'./Analytics/ldf_area_theory.txt' u 1:($2+0.06) w l lw 2 dt 2 lc rgb "#28b463"
