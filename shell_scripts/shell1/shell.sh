#----------------------------------------------------------------------------------------#
# This shell script performs several runs of runge_kutta.dat for various values          #
# of electron energy, barrier width and barrier height. This data is then collected      #
# in a .csv file which can be used to represent the data graphically                     #
#----------------------------------------------------------------------------------------#

#!/bin/bash

#Width in Bohr
width_min=1000
width_max=5000
width_n=10
width_h=$(echo "scale=0;($width_max - $width_min)/$width_n" | bc -l)

#Height in Hartrees
height_min=0.005
height_max=0.05
height_n=10
height_h=$(echo "($height_max - $height_min)/$height_n" | bc -l)

#Electron energy in Hartrees
energy_min=$(echo "0.5*$height_min" | bc -l)
energy_max=$(echo "0.1*$height_max + $height_max" | bc -l)
energy_n=1
energy_h=$(echo "($energy_max - $energy_min)/$energy_n" | bc -l)

#Arrays containing data to be written to a file
energy_array1=()
energy_array2=()

for (( i=1; i<=$width_n; i++))
do
	width=$(echo "$width_min + ($i-1)*$width_h" | bc -l)
	width_start=$(echo "$width - 0.5*$width" | bc -l)
	width_end=$(echo "$width + 0.5*$width" | bc -l)
	boundary=$(echo "scale=0;2*$width_end/1" | bc -l)

	sed -r -i "s/([0-9]+\.*[0-9]*\s+)([0-9]+\.*[0-9]*\s+)([0-9]+\.*[0-9]*\s+)(::\s+height)/\1 $width_start $width_end \4/" input2.dat
	sed -r -i "s/([0-9]+\.*[0-9]*\s+)([0-9]+\.*[0-9]*\s+)(::\s+start\s+finish)/0 $boundary \3/" input1.dat

	for (( x=1; x<=$energy_n; x++))
	do
		energy=$(echo "$energy_min + ($x-1)*$energy_h" | bc -l)
		energy_array1+=("$energy ,")
	done

	echo "$width ${energy_array1[*]}" >> EvH_data_${width}.csv
	energy_array1=()

	for (( j=1; j<=$height_n; j++))
	do
		height=$(echo "$height_min + ($j-1)*$height_h" | bc -l)

		sed -r -i "s/([0-9]*\.[0-9]*\s+)([0-9]*\.[0-9]*\s+)([0-9]*\.[0-9]*\s+)(::\s+height)/$height \2\3\4/" input2.dat

		for ((k=1; k<=$energy_n; k++))
		do
			energy=$(echo "$energy_min + ($k-1)*$energy_h" | bc -l)

			sed -r -i "s/([0-9]*\.*[0-9]*\s+)(::\s+energy)/$energy \2/" input2.dat 

			temp=$(./runge_kutta.out)
			echo "$height $width $temp"
			echo "$height $width $temp" >> data.dat
			energy_array2+=("$temp ,")

		done

		echo "$height ,${energy_array2[*]}" >> EvH_data_${width}.csv
		energy_array2=()

	done
done
