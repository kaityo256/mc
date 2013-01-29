set term png
set out "nucleation.png"

set log y
set xla "Inverse Temperature"
set yla "Nucleation Time [MCs]"
unset key
p [:51] [50:1e4] "nucleation.dat" w e pt 6

