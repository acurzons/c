gnuplot<<EOF
set log xy
set ytics format '%.1e'
set yrange [3E-14:1E-9]
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})'
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'Domainkolep1.ps'
set style line 1 lt 1 lw 1 lc rgb 'blue'
plot  'OUTPUT/Domainkoelectron_syncprimary.txt' u 1:6 w lines title 'synchrotron', 'OUTPUT/Domainkoelectron_ICprimary.txt' u 1:6 w lines title 'IC CMB',  'Domainko1.dat' u (\$1*1E-12):2  title 'Domainko lepton {/Symbol a}=2' ls 1 
EOF