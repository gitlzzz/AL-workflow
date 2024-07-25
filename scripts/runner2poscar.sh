#convert the structures need to be run in DFT to POSCAR format
module load gcc/8.1.0 impi/2018.1 mkl/2018.1 opencv/4.1.2 python/3.6.4_ML
RuNNerUC.py runner vasp_poscar -i $1 -o POSCAR
