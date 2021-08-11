#This is a GNUPLOT script to generate plots of center of the particle over time
#Open the terminal in the window of this file and use command "gnuplot plot.p" 


#Setup common properties
set datafile separator ","
set size 1,1
set terminal pngcairo dashed size 1200,600 
set xlabel "t/t_{ref}"
set ylabel "{/Symbol e}/D_p"
set key top right

#Plot St = 27
set output 'st_27.png'
set xrange [-0.3:2]
plot \
"../Reference_data/mobile_st_27.dat" u (($1-0.1667)/0.024731):($4/0.006-0.5) w l lc 1 lw 3 dt (10,3) title "Reference data",\
"../Reference_data/gondret27_expy.dat" u ($1/0.024731):($2/0.006) w p pt 6 ps 3 lc "black" lw 2 title "Gondret, exp. data",\
"../Reference_data/st_27_prev_num.csv" u ($1):($2) w l lc "black" lw 1 title "JCP(2017) num data"

#Plot St = 152
set output 'st_152.png'
set xrange [-0.5:10]
plot \
"../Reference_data/mobile_st_152.dat" u (($1/0.017487424)-20.6):(($4/0.003)-0.5) w l lc 1 lw 3 dt (10,3) title "Reference data",\
"../Reference_data/gondret152_expy.dat" u ($1/0.017487424):($2/0.003) w p pt 6 ps 2 lc "black" lw 2 title "Gondret, exp. data",\
"../Reference_data/st_152_prev_num.csv" u ($1):($2) w l lc "black" lw 1 title "JCP(2017) num data"









