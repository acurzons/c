gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 1 lc rgb 'blue'
set style line 2 lt 1 lw 1 lc rgb  'black'

set ytics format '%.1e'
set yrange [1E-16:1E-8]
set xrange [1E-15:1E5]
set label 1 at 1E-4,1E-10 'W_{p}=8.0x10^{50} erg'
set label 2 at 1E-4,1E-11 'n_{H}=50 cm^{-3}'                            
set label 3 at 1E-4,3E-12 'proton spectra {/Symbol a}=2.2'             
set label 4 at 1E-4,3E-11 'E_{cut}=200 TeV''
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' 
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'SanoDorcelec2.ps'
f(x)=(x>0.5E-9 && x<9E-9)? 1.56E-17*(x/1.6)**(-0.5):1/0
g(x)=(x>0.5E-9 && x<9E-9)? 8.96E-16*(x/1.6)**(-0.3):1/0
h(x)=(x>0.5E-9 && x<9E-9)? 2.66E-19*(x/1.6)**(-0.7):1/0


plot 'OUTPUT/SanoDorcelec_syncprimary.txt' u 1:6 w lines title 'Synchrotron',  'OUTPUT/SanoDorcelec_ICprimary.txt' u 1:6 w lines title 'IC CMB',  'OUTPUT/SanoDorcelec_ICRprimary.txt' u 1:6 w lines title 'IC IR', '< paste OUTPUT/SanoDorcelec_ICprimary.txt OUTPUT/SanoDorcelec_ICRprimary.txt' u 1:(\$6 +\$12) w lines title 'IC total', 'DorC.txt' u 1:(\$2*\$1*\$1):(\$3*\$1*\$1):(\$4*\$1*\$1) with yerrorbars title 'DorC HESS' ls 1,  'DorC.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1):(0):(-0.8*\$2*\$1*\$1) w vector title '' ls 1, "<echo '3E-3 7.64E-13 1E-3 1E-2 0' " u 1:2:3:4:(\$2):(\$2) w xyerrorbars notitle ls 2, '' u 1:(\$2):(0):(-0.65*\$2) w vectors notitle ls 2, '+' using 1:(g(\$1)):(h(\$1)) w filledcurves title 'Suzaku'


EOF