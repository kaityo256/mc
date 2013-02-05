set term png 
set out "index.png"

set yra [-1:10]
set style data linespoints
set pointsize 1.5
unset ytics
unset key

set xla "MCs"

set label "9" at 9.2, 0
set label "8" at 9.2, 1
set label "7" at 9.2, 2
set label "6" at 9.2, 3
set label "5" at 9.2, 4
set label "4" at 9.2, 5
set label "3" at 9.2, 6
set label "2" at 9.2, 7
set label "1" at 9.2, 8
set label "0" at 9.2, 9

p "index.dat" u 1:2 w l\
 ,"index.dat" u 1:3 w l\
 ,"index.dat" u 1:4 w l\
 ,"index.dat" u 1:5 w l\
 ,"index.dat" u 1:6 w l\
 ,"index.dat" u 1:7 w l\
 ,"index.dat" u 1:8 w l\
 ,"index.dat" u 1:9 w l\
 ,"index.dat" u 1:10 w l\
 ,"index.dat" u 1:11 w l\
 ,"index.dat" u 1:2 pt 6 lc 1\
 ,"index.dat" u 1:3 pt 6 lc 2\
 ,"index.dat" u 1:4 pt 6 lc 3\
 ,"index.dat" u 1:5 pt 6 lc 4\
 ,"index.dat" u 1:6 pt 6 lc 5\
 ,"index.dat" u 1:7 pt 6 lc 6\
 ,"index.dat" u 1:8 pt 6 lc 7\
 ,"index.dat" u 1:9 pt 6 lc 8\
 ,"index.dat" u 1:10 pt 6 lc 9\
 ,"index.dat" u 1:11 pt 6 lc 10
