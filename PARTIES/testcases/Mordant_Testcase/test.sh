#!/bin/bash

#echo 'Hello World'

#g=9.81



#awk 'NR==1{print}' mobile.dat

#awk -F'\t' '{ print $7 }'

#reference=$(awk 'BEGIN{FS=","} {print $8}' mobile.dat)

#for ((i=0; i<10; ++i))

#for i in awk 'BEGIN{FS=","} {print $8}' mobile.dat
#for i in $reference
	#a=$(( 2*i ))
#test = $reference*2

#echo $i

#python3 test.py
#done

awk '{print $2}' mobile.dat