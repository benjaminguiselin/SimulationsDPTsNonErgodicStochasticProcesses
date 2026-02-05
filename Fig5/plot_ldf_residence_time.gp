set term qt
set xlabel 'q'
set ylabel '{/Symbol f}(q)'
set g
set key top center
set samples 10000
p './Data/ldf_residence_time.txt' w p title 'simu', \
'./Analytics/ldf_residence_time_theory.txt' w l lc rgb "#000000" title 'asymptotic', \
'./Data/ldf_residence_time.txt' every 39::1 w p title 'simu (short)'
pause -1

set term cairolatex pdf font ',8' linewidth 2 size 8.5cm, 5.3cm standalone header "\\usepackage{amsmath,bm}"
set output 'ldf_residence_time_mortal.tex'
set xlabel '$q$'
set ylabel '$\phi(q)$' offset char -0.1, 0
set xr [0:1]
set yr [0:3.]
unset key
unset g
set ytics 1
set mytics 2
set xtics 0.2
set mxtics 2
set arrow 1 from 0.7873276081682232, 0 to 0.7873276081682232, 1.3700869421974955 lw 3 dt 2 nohead
set label 1 at 0.8, 0.25 '$q_c$'
p './Data/ldf_residence_time.txt' every 39::1 w p pt 1 ps 0.2 lc rgb "#ffffff", \
'./Analytics/ldf_residence_time_theory.txt' index 0 w l lw 3 lc rgb "#1D60FF", \
'./Analytics/ldf_residence_time_theory.txt' index 1 w l lw 3 lc rgb "#FF0380", \
'./Data/ldf_residence_time.txt' every 39::1 w p pt 7 ps 0.3 lc rgb "#000000"
