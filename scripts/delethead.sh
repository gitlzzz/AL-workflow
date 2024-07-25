# delete the distractor for active learning
filename=$1
sed -i '/comment vasp./d' $filename
sed -i '/comment VRHFIN/d' $filename
sed -i '/comment TITEL/d' $filename
sed -i '/comment VRHFIN/d' $filename
sed -i '/comment k-points/d' $filename
sed -i '/comment PREC/d' $filename
sed -i '/comment ISPIN/d' $filename
sed -i '/comment LASPH/d' $filename
sed -i '/comment ENCUT/d' $filename
sed -i '/comment ISMEAR/d' $filename
sed -i '/comment GGA/d' $filename
sed -i '/comment IVDW/d' $filename
sed -i '/comment VDW_SR/d' $filename
sed -i '/comment Values for/d' $filename
