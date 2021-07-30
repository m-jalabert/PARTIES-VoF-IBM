set datafile separator ","
set size 1,1

set terminal pngcairo dashed size 1200,600 

set xrange [-0.5:2]
#set yrange [-0.2:1]
set key top right
set key outside

set output 'st_27.png'
set ylabel "Y"
plot \
"mobile_st_27.dat" u (($1-0.1667)/0.024731):($4/0.006-0.5) w l lc 1 lw 2 title "PARTIES",\
"../../extracted_data/st_27.csv" u ($1):($2) w p pt 6 ps 2 lc "black" lw 2 title "Gondret",\
"../../extracted_data/st_27_num.csv" u ($1):($2) w l lc "black" lw 1 title "Previous numerical data"


set datafile separator ","
set size 1,1

set terminal pngcairo dashed size 1200,600 

set xrange [-0.5:10]
#set yrange [-0.2:1]
set key top right
set key outside

set output 'st_152.png'
set ylabel "Y"
plot \
"mobile_st_152.dat" u (($1/0.017487424)-20.6):(($4/0.003)-0.5) w l lc 1 lw 2 title "PARTIES",\
"../../extracted_data/st_152.csv" u ($1):($2) w p pt 6 ps 2 lc "black" lw 2 title "Gondret",\
"../../extracted_data/st_152_num.csv" u ($1):($2) w l lc "black" lw 1 title "Previous numerical data",\









