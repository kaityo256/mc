set term pngcairo
set out "q3.png"

Q=3
bc = log(1.0+sqrt(Q))

set xla "Beta"
set yla "U"
set yra [1:2.2]
set style data line
set arrow 1 from bc,1.25 to bc, 1
p "Q3_L16.dat" u 1:3 t "L=16"\
, "Q3_L32.dat" u 1:3 t "L=32"\
