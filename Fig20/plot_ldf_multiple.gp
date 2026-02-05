set term qt
set xlabel 'q'
set ylabel '{/Symbol f}(q)'
set g
set key top center
set samples 10000
set xr [0:16]
p './Data/ldf_3walker_full.txt' w p title 'simu', \
(x < 2)?x:((x<2*sqrt(3))?(0.25*x**2+1):((x<4*sqrt(3))?(x*sqrt(3)-2):((x<4*sqrt(5))?(0.125*x**2+4):((x<6*sqrt(5))?(x*sqrt(5)-6):(x**2/12.+9.))))) w l lc rgb "#000000" title 'asymptotic', \
'./Data/ldf_3walker_full.txt' every 55::10 w p title 'simu (short)'
pause -1

set term cairolatex pdf font ',8' linewidth 2 size 8.5cm, 5.3cm standalone header "\\usepackage{amsmath,bm}"
set output 'ldf_multiple_mortal.tex'
set xlabel '$q$'
set ylabel '$\phi(q)$' offset char -.1, 0
set xr [0:16]
set yr [0:32]
unset key
unset g
set samples 1000
set ytics 8
set mytics 2
set xtics 2
set mxtics 2
set arrow 1 from 2, 0 to 2, 2 lw 3 dt 2 nohead
set arrow 2 from 2*sqrt(3), 0 to 2*sqrt(3), 4 lw 3 dt 2 nohead
set arrow 3 from 4*sqrt(3), 0 to 4*sqrt(3), 10 lw 3 dt 2 nohead
set arrow 4 from 4*sqrt(5), 0 to 4*sqrt(5), 14 lw 3 dt 2 nohead
set arrow 5 from 6*sqrt(5), 0 to 6*sqrt(5), 24 lw 3 dt 2 nohead
p './Data/ldf_3walker_full.txt' every 55::10 w p pt 1 ps 0.2 lc rgb "#ffffff", \
[0:2] x w l lw 3 lc rgb "#1D60FF", \
[2:2*sqrt(3)] 0.25*x**2+1 w l lw 3 lc rgb "#FF0380", \
[2*sqrt(3):4*sqrt(3)]x*sqrt(3)-2 w l lw 3 lc rgb "#1D60FF", \
[4*sqrt(3):4*sqrt(5)]0.125*x**2+4 w l lw 3 lc rgb "#FF0380", \
[4*sqrt(5):6*sqrt(5)] x*sqrt(5)-6 w l lw 3 lc rgb "#1D60FF", \
[6*sqrt(5):] x**2/12+9 w l lw 3 lc rgb "#FF0380", \
'./Data/ldf_3walker_full.txt' every 55::10 w p pt 7 ps 0.3 lc rgb "#000000"