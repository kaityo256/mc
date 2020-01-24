set term png
set out "m2.png"
set xlabel "Temperature"
set ylabel "m^2"

p "L064.dat" u 1:2:4 w e t "Normal"\
, "L064.dat" u 1:5:7 w e t "Improved Estimator"\

set out "m4.png"
set ylabel "m^4"

p "L064.dat" u 1:8:10 w e t "Normal"\
, "L064.dat" u 1:11:13 w e t "Improved Estimator"\

set out "binder.png"
set ylabel "U"

p "L064.dat" u 1:14:16 w e t "Normal"\
, "L064.dat" u 1:17:19 w e t "Improved Estimator"\
