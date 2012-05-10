set term png 
set out "acf.png"
set xla "MCs"
set yla "ACF"

set style data line
set xra[:100]
set yra[0:1]
p "single.acf" t "Single-Flip" lw 3,\
  "sw.acf" t "Swendsen-Wang" lw 3,\
  "wolff.acf" t "Wolff" lw 3

