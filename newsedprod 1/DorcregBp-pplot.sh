gnuplot<<EOF
set log xy
 set style line 1 lt 1 lw 1 lc rgb 'blue'
 set style line 2 lt 1 lw 1 lc rgb  'black'
set ytics format '%.1e'
set yrange [1E-17:1E-9]
set xrange [1E-5:1E5]
#set label 1 at 1E-4,1E-10 'W_{p}=7.0x10^{47} erg'
#set label 2 at 1E-4,1E-11 'n_{H}=50 cm^{-3}'                            
#set label 3 at 1E-4,3E-12 'proton spectra {/Symbol a}=2.2'             
#set label 4 at 1E-4,3E-11 'E_{cut}=200 TeV''
#set label 5 at 3E-4,1E-13 'Fermi-LAT'
set ylabel 'E_{/Symbol g}^{2} dN/dE_{/Symbol g} (erg cm^{-2} s^{-1})' 
set xlabel 'E_{/Symbol g} (TeV)'
set term postscript enhanced color
set output 'regBDorc.ps'

#####------FERMI functions-----#######

plot 'OUTPUT/DorcregB_ppgamma.txt' u 1:6 w lines title 'hadronic scenario', 'DorC.txt' u 1:(\$2*\$1*\$1*0.14):(\$3*\$1*\$1*0.14):(\$4*\$1*\$1*0.14) with yerrorbars title 'DorC HESS' ls 1,  'DorC.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1*0.14):(0):(-0.8*\$2*\$1*\$1*0.14) w vector title '' ls 1,"<echo '3E-3 7.64E-13 1E-3 1E-2 0' " u 1:(\$2*0.14):3:4:(\$2*0.14):(\$2*0.14) w xyerrorbars notitle ls 2, '' u 1:(\$2*0.14):(0):(-0.65*\$2*0.14) w vectors notitle ls 2, 'DorC.txt' u 1:(valid(3)? 1/0: \$2*\$1*\$1*0.14):(0):(-0.8*\$2*\$1*\$1*0.14) w vector title '' ls 1



EOF