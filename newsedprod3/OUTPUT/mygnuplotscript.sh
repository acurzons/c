#!/bin/sh

outfile=$1
infile1=$2
infile2=$3
infile3=$4
gnuplot <<EOF

set log xy
set yrange [1.0e-60:1]
set ytics format "%.1e"
set xtics format "%.1e"
set term postscript enhanced color 
set output '${outfile}'

plot '${infile1}_ppgamma.txt' u 1:6, '${infile1}_syncprimary.txt' u 1:6, '${infile2}_ppgamma.txt' u 1:6, '${infile2}_syncprimary.txt' u 1:6, '${infile3}_ppgamma.txt' u 1:6, '${infile3}_syncprimary.txt' u 1:6      

#plot '${infile}_ppgamma.txt' u 1:6, '${infile}_syncprimary.txt' u 1:6, '${infile}_ICprimary.txt' u 1:6, '${infile}_ICRprimary.txt' u 1:6, '${infile}_bremssprimary.txt' u 1:6 

EOF