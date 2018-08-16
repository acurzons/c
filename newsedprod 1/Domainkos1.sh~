gnuplot<<EOF
set log xy
set ytics format '%.1e'
set yrange [1E-15:1E-10]
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})'
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'AharonianspecIC.ps'
set style line 1 lt 1 lw 1 lc rgb 'blue'
plot 'OUTPUT/Aharonianspec_ICprimary.txt' u 1:6 w lines title 'IC CMB', 'OUTPUT/Aharonianspec_ICRprimary.txt' u 1:6 w lines title ' IC FIR', 'Ahaspec24ICCMB.dat' u 1:2 notitle ls 1 , 'Ahaspec24ICFIR.dat' u 1:2 title 'Aharonian 1997' ls 1
EOF