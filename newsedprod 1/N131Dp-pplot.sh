gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 1 lc rgb 'blue'
set ytics format '%.1e'
set yrange [1E-17:1E-9]
set xrange [1E-5:1E5]
set label 1 at 1E-4,1E-10 'W_{p}=1.4x10^{50} erg'
set label 2 at 1E-4,1E-11 'n_{H}=50 cm^{-3}'                            
set label 3 at 1E-4,3E-12 'proton spectra {/Symbol a}=2.0'             
set label 4 at 1E-4,3E-11 'E_{cut}=100 TeV''
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' 
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'SanoN132D.ps'
a(x)=((x>0.1E-3)&&(x<100E-3))? 3.6E-12*x**(0.3):1/0
b(x)=((x>0.1E-3)&&(x<100E-3))? 3.8E-11*x**(0.9):1/0

c(x)=((x>0.1E-3)&&(x<100E-3))? 6.0E-13*x**(0.3):1/0
d(x)=((x>0.1E-3)&&(x<100E-3))? 6.30E-12*x**(0.9):1/0

max1(x,y)=(x>y) ? x:y
min1(x,y)=(x<y) ? x:y

plot '+' using 1:(max1(a(\$1),b(\$1))):(min1(c(\$1),d(\$1))) w filledcurves closed  title 'Fermi', 'OUTPUT/SanoN131D_ppgamma.txt' u 1:6 w lines title 'hadronic scenario', 'N132D.txt' u 1:(\$2*\$1*\$1):(\$3*\$1*\$1):(\$4*\$1*\$1) with yerrorbars title 'N132D HESS' ls 1,  'N132D.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1):(0):(-0.8*\$2*\$1*\$1) w vector title '' ls 1


EOF