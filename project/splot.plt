set term pngcairo size 1000, 750
#set termoption dash

set style line 1 lt 1 lc rgb "red" lw 2
set style line 2 lt 2 lc rgb "#00BB00" lw 2
set style line 3 lt 3 lc 3 lw 2
set style line 4 lt 4 lc rgb "#FFD700" lw 2
set style line 5 lt 5 lc rgb "#BB00BB" lw 2
set style line 6 lt 2 lc rgb "orange" lw 3
set style line 7 lt 2 lc rgb "red" lw 3
set style line 8 lt 4 lc rgb "greenyellow" lw 2
set style line 9 lt 4 lc rgb "brown" lw 3
set style line 10 lt 2 lc rgb "black" lw 2
set style line 11 lt 3 lc rgb "gold" lw 3


set output "Solution.png"
splot "solution.dat" u 1:2:3 with lines ls 2

