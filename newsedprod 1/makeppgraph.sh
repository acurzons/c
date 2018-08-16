infile=$1
outfile=$2
gnuplot<<EOF
set log xy
set xlabel 'E_{/Symbol g} (TeV)'
set ylabel 'E^2dN/dE (erg cm^{-2} s^{-1})
set ytics format '10^{%T}'
set xrange [0.0001:1E5]
set yrange [10E-16:10E-8]
set term postscript enhanced
set output '${outfile}
plot '${infile}' u 1:6 


EOF
