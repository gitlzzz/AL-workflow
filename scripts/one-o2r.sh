#one-o2r.sh OUTCAR, convert singlepoint OUTCAR to RuNNer .data file
dir=$(pwd);RuNNerUC.py vasp_outcar runner -i $1 -o $1.data -c "singlepoint frame 1 $dir" -f 1
