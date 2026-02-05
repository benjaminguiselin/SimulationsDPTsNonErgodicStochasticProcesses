set term qt
set xlabel 'q/N'
set ylabel '{/Symbol f}(q)/N^2'
set g
set key top center
set samples 10000
set yr [0:1]
set xr [0:1]
p './Data/One_walker/ldf_1walker.txt' w p title '1 walker', \
'./Data/Two_walkers/ldf_2walker.txt' u ($1/2):($2/4) w p title '2 walkers', \
'./Data/Three_walkers/ldf_3walker.txt' u ($1/3):($2/9) w p title '3 walkers', \
x w l lc rgb "#000000" title 'asympotic (1 walker)', \
(x<0.5)?(x/2):((3*x-1)/2) w l lc rgb "#000000" title 'asympotic (2 walkers)', \
(x<1./3.)?(x/3):((x<2./3.)?(x-2./9.):((5*x-2)/3)) w l lc rgb "#000000" title 'asympotic (3 walkers)'
pause -1

set term cairolatex pdf font ',8' linewidth 2 size 8.5cm, 5.3cm standalone header "\\usepackage{amsmath,bm}"
set output 'ldf_epidemiology.tex'
set xlabel '$q/N$'
set ylabel '$\phi(q)/N^2$' offset char -0.1,0
set xr [0:1]
set yr [0:1]
unset key
unset g
set ytics 0.2
set mytics 2
set xtics 0.2
set mxtics 2
set arrow 1 from .5, 0 to .5, .25 lw 3 dt 2 lc rgb "#28b463" nohead 
set arrow 2 from 1./3., 0 to 1./3., 1./9. lw 3 dt 2 lc rgb "#FF0380" nohead
set arrow 3 from 2./3., 0 to 2./3., 4./9. lw 3 dt 2 lc rgb "#FF0380" nohead
set samples 10000
p './Data/One_walker/ldf_1walker.txt' every 20 w p pt 1 ps 0.2 lc rgb "#ffffff", \
x w l lw 2 lc rgb "#1D60FF", \
(x<0.5)?(x/2):((3*x-1)/2) w l lw 2 lc rgb "#28b463", \
(x<1./3.)?(x/3):((x<2./3.)?(x-2./9.):((5*x-2)/3)) w l lw 2 lc rgb "#FF0380", \
x**2 w l lw 2 dt (5,8) lc rgb "#707b7c", \
'./Data/One_walker/ldf_1walker.txt' every 4 w p pt 7 ps 0.3 lc rgb "#000000", \
'./Data/Two_walkers/ldf_2walker.txt' every 8 u ($1/2):($2/4) w p pt 7 ps 0.3 lc rgb "#000000", \
'./Data/Three_walkers/ldf_3walker.txt' every 12 u ($1/3):($2/9) w p pt 7 ps 0.3 lc rgb "#000000"
