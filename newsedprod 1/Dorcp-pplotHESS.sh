gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 2 lc rgb 'blue'
 set style line 2 lt 1 lw 2 lc rgb  'black'
set style line 3 lt 4 lw 4 lc rgb 'purple'
set style line 4 lt 2 lw 2 lc rgb 'blue'
set style line 5 lt 1 lw 2 lc rgb 'red'
set ytics format '%.1e'
set yrange [1E-17:1E-9]
set xrange [1E-5:1E5]
set label 1 at 1E-1,3E-10 'W_{p}=1.2x10^{50} erg' font ',15'
set label 2 at 1E-1,1E-11 'n_{H}=60 cm^{-3}' font ',15'                            
set label 3 at 1E-1,1E-10 'proton spectra {/Symbol a}=2.0' font',15'             
set label 4 at 1E-1,3E-11 'E_{cut}=100 TeV' font ',15'
set label 5 at 6E-5,6E-11 '30 DorC Hadronic scenario' font ',15'
set label 6 at 200,2E-12 'CTA south 100 h' font ',15' textcolor rgb 'purple'
set arrow 1 nohead from 0.8E-1,3E-10 to 0.8E-1,1E-11 lc rgb 'black' 
#set label 5 at 3E-4,1E-13 'Fermi-LAT'
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' font ',15' 
set xlabel 'E_{/Symbol g} (TeV)'font ',15'
set xtics font ',15'
set ytics font ',15'
set key font ',15'
set term postscript enhanced color
set output 'SanoDorchess.ps'
i(x)=(x>0.2E-3 && x<200E-3)? 1.6*9.53E-13*x**(-0.05):1/0
j(x)=(x>0.2E-3 && x<200E-3)? 1.6*7.62E-13*x**(-0.05):1/0
k(x)=(x>0.2E-3 && x<200E-3)? 1.6*4.37E-13*x**(-0.15):1/0
l(x)=(x>0.2E-3 && x<200E-3)? 1.6*5.46E-13*x**(-0.15):1/0

ismin(x1,x2,x3,x4)=(x1<x2 && x1<x3 && x1<x4  ? x1: x2<x1 && x2<x3 && x2<x4  ? x2: x3<x1 && x3<x2 && x3<x4 ? x3: x4<x1 && x4<x2 && x4<x3 ? x4:1/0)

ismax(x1,x2,x3,x4)=(x1>x2 && x1>x3 && x1>x4  ? x1: x2>x1 && x2>x3 && x2>x4  ? x2: x3>x1 && x3>x2 && x3>x4 ? x3: x4>x1 && x4>x2 && x4>x3 ? x4:1/0)
#####------FERMI functions-----#######

plot 'OUTPUT/Sanopaperhadron_ppgamma.txt' u 1:6 w lines title '30 DorC','OUTPUT/LMCcoreMC3_ppgamma.txt' u 1:6  w lines ls 4 title 'core MC3', 'DorC.txt' u 1:(\$2*\$1*\$1*1.6):(\$3*\$1*\$1*1.6):(\$4*\$1*\$1*1.6) with yerrorbars title 'DorC HESS' ls 2,  'DorC.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1):(0):(-0.8*\$2*\$1*\$1) w vector title '' ls 2,  '+' using 1:(ismin(i(\$1),j(\$1),k(\$1),l(\$1))):(ismax(i(\$1),j(\$1),k(\$1),l(\$1))) w filledcurves ls 2 title 'Fermi E1+E3', 'CTAsouthsens50h.dat' u 1:(\$2/1.4) w lines ls 3 title '' 

#,"<echo '3E-3 7.64E-13 1E-3 1E-2 0' " u 1:2:3:4:(\$2):(\$2) w xyerrorbars notitle ls 2, '' u 1:(\$2):(0):(-0.65*\$2) w vectors notitle ls 2, 'DorC.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1):(0):(-0.8*\$2*\$1*\$1) w vector title '' ls 1

EOF