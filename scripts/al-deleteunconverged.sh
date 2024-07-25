##usage: check the convergence from OSZICAR, and delete the unconverged loop from the .data file
filename=$1  #.data file
inputnelw=$2  #nelw number in INCAR
nelw=${inputnelw:-60} # Set default value of nelw to 60 if inputnelw is not provided
echo $nelw
cp $filename $filename-all 

# Extract unconverged loop numbers based on nelw from OSZICAR file, and Sorting them in descending order, to make sure delete the line from bottom to top to aviod the change of line number
#NELW=$(bzgrep NELM OUTCAR.1.bz2 | awk {'print substr($3, 1, length($3)-1)'} )
#Nloops=$(bzgrep "RMM:  $NELW" -A1 OSZICAR.1.bz2 | grep T| awk {'print $1'}|sort -nr | cut -f2) ###number of unconverged loops
Nloops=$(grep -a "RMM:\+[[:space:]]\+$nelw" -A1 OSZICAR | grep F| awk {'print $1'}|sort -nr | cut -f2);echo "make sure your NELW is $nelw"
#Natoms=$1 #number of atoms
Natoms=$(grep -a NIONS ${filename%.data} | awk {' print $12'})

# Iterate over each unconverged loop number
for N in $Nloops;
do
#I=$((14+(9+Natoms)*(N-1)))
#F=$((13+(9+Natoms)*N))
# Find the line numbers containing "begin" and "end" for the current loop number
grepI=`grep -a -n begin $filename|head -$N|tail -1`
grepF=`grep -a -n end $filename|head -$N|tail -1`
I=${grepI%:begin} # to get the starting line number
F=${grepF%:end}  # to get the ending line number
echo "########deleting unconverged loops $N as below"
#sed -n "${I},${F}p" $filename
#if [ $F -gt $(wc $filename | awk {'print $1'}) ];then
#echo "error,;Loops in input.data is less than the unconverged loop number"
#exit
#fi

# Delete the lines between the starting and ending line numbers in the filename
sed -i "${I},${F}d" $filename
done
mv $filename $filename-deleunconverged
