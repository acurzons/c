gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 1 lc rgb 'blue'
set ytics format '%.1e'
set yrange [1E-17:1E-9]
set xrange [1E-5:1E5]
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' 
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'kelnervskafexhiu.ps'
plot 'OUTPUT/kafexhiu_ppgamma.txt' u 1:6 w lines title 'kafexhiu', 'OUTPUT/kelner_ppgamma.txt' u 1:(\$6) w lines title 'kelner'

EOF