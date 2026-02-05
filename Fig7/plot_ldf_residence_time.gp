set term qt
set xlabel 'q'
set ylabel '{/Symbol f}(q)'
set g
set key top center
set samples 10000
p './Data/ldf_residence_time.txt' w p title 'simu', \
'./Analytics/ldf_residence_time_theory.txt' w l lc rgb "#000000" title 'asymptotic', \
'./Data/ldf_residence_time.txt' every 25::1 w p title 'simu (short)'
pause -1

set term cairolatex pdf font ',8' linewidth 2 size 8.5cm, 5.3cm standalone header "\\usepackage{amsmath,bm}"
set output 'ldf_residence_time_absorbing.tex'
set multiplot
set xlabel '$q$'
set ylabel '$\phi(q)$' offset char -0.1, 0
set xr [0:1]
set yr [0:20.]
unset key
unset g
set ytics 4
set mytics 2
set xtics 0.2
set mxtics 2
set arrow 1 from 0.4081301395778893, 0 to 0.4081301395778893, 2.7187854367500845 lw 3 dt 2 nohead
set label 1 at 0.42, 1. '$q_{c}$'
p './Data/ldf_residence_time.txt' every 25::1 w p pt 1 ps 0.2 lc rgb "#ffffff", \
'./Analytics/ldf_residence_time_theory.txt' index 0 w l lw 3 lc rgb "#1D60FF", \
'./Analytics/ldf_residence_time_theory.txt' index 1 w l lw 3 lc rgb "#FF0380", \
'./Data/ldf_residence_time.txt' every 25::1 w p pt 7 ps 0.3 lc rgb "#000000"
set size 0.5, 0.55
set origin 0.12, 0.4
set xlabel '$p$'
set ylabel '$\mu(p)$'
set xr [0:14]
set yr [-1:3]
set xtics 2
set ytics 1
set mxtics 2
set mytics 2
unset arrow
unset label
set label 2 at 6.661564959554329, -0.3 '$p_c$'
p './Analytics/scgf_residence_time_theory.txt' w l lw 2 lc rgb "#000000"
unset multiplot
