gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 1 lc rgb 'blue'
set style line 2 lt 1 lw 1 lc rgb  'black'

set ytics format '%.1e'
set yrange [1E-16:1E-8]
set xrange [1E-15:1E5]
set label 1 at 1E-4,1E-10 'W_{e}=7.0x10^{48} erg'
set label 2 at 1E-4,1E-11 'B=13 {/Symbol m}G '                            
set label 3 at 1E-4,3E-12 'electron spectra {/Symbol a}=2.3'             
set label 4 at 1E-4,3E-11 'E_{cut}=800 TeV'
set label 5 at 1E-4,1E-12 't_{age}=800 yr'
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' 
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'SanoregCDorcelectime2.ps'
f(x)=(x>0.5E-9 && x<9E-9)? 3.0E-16*(x)**(-0.3):1/0
g(x)=(x>0.5E-9 && x<9E-9)? 4.9E-17*(x)**(-0.4):1/0
h(x)=(x>0.5E-9 && x<9E-9)? 2.0E-14*(x)**(-0.1):1/0


plot 'OUTPUT/regCDorcelectime_syncprimary.txt' u 1:6 w lines title 'Synchrotron',  'OUTPUT/regCDorcelectime_ICprimary.txt' u 1:6 w lines title 'IC CMB',  'OUTPUT/regCDorcelectime_ICRprimary.txt' u 1:6 w lines title 'IC IR', 'OUTPUT/regCDorcelectime_ICOprimary.txt' u 1:6 w lines title 'IC IR2', '< paste OUTPUT/regCDorcelectime_ICprimary.txt OUTPUT/regCDorcelectime_ICRprimary.txt OUTPUT/regCDorcelectime_ICOprimary.txt' u 1:(\$6 +\$12 +\$18 ) w lines title 'IC total', 'DorC.txt' u 1:(\$2*\$1*\$1*0.14*1.6):(\$3*\$1*\$1*0.14*1.6):(\$4*\$1*\$1*0.14*1.6) with yerrorbars title 'DorC HESS' ls 1,  'DorC.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1*0.14*1.6):(0):(-0.8*\$2*\$1*\$1*0.14*1.6) w vector title '' ls 1, '+' using 1:(g(\$1)):(h(\$1)) w filledcurves title 'Suzaku'


EOF