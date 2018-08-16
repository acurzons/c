gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 1 lc rgb 'blue'
set ytics format '%.1e'
set yrange [1E-16:1E-8]
set xrange [1E-15:1E5]
set label 1 at 1E-4,1E-10 'W_{p}=1.4x10^{50} erg'
set label 2 at 1E-4,1E-11 'n_{H}=50 cm^{-3}'                            
set label 3 at 1E-4,3E-12 'proton spectra {/Symbol a}=2.0'             
set label 4 at 1E-4,3E-11 'E_{cut}=100 TeV''
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' 
set xlabel 'E_{/Symbol g} (TeV)'
#f(x)= 6.69E-20*x**(-0.86)
f(x)=(x>0.3E-9 && x<8.0E-9) ? 6.69E-20*x**(-0.86) : 1/0
g(x)=(x>0.3E-9 && x<8.0E-9) ? 1.10E-6*x**(0.67) : 1/0
set term postscript enhanced color
set output 'SanoN132Delec.ps'
plot 'OUTPUT/SanoN131Delec_syncprimary.txt' u 1:6 w lines title 'Synchrotron', 'OUTPUT/SanoN131Delec_ICprimary.txt' u 1:6 w lines title 'IC CMB',  'OUTPUT/SanoN131Delec_ICRprimary.txt' u 1:6 w lines title 'IC IR', '< paste OUTPUT/SanoN131Delec_ICprimary.txt OUTPUT/SanoN131Delec_ICRprimary.txt' u 1:(\$6 +\$12) w lines title 'IC total' , 'N132D.txt' u 1:(\$2*\$1*\$1):(\$3*\$1*\$1):(\$4*\$1*\$1) with yerrorbars title 'N132D HESS' ls 1,  'N132D.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1):(0):(-0.8*\$2*\$1*\$1) w vector title '' ls 1, '+' using 1:(f(\$1)):(g(\$1)) w filledcurves closed  title 'Chandra'


EOF