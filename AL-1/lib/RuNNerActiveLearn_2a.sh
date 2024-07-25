####################################################################################################

####################################################################################################

# RuNNerActiveLearn_2a.sh
# written by Marco Eckhoff (26th May 2021)
# modified by Zan 30 Apr 2022
# recalculates selected structures employing RuNNer and the two HDNNPs.
#
# No modifications are required for regular usage.

####################################################################################################

####################################################################################################

#RuNNer_exe=$1
energy_threshold=$1
force_threshold=$2
mkdir mode2
mkdir mode2/nnp-data-1 mode2/nnp-data-2

cp RuNNer/HDNNP_1/input.nn mode2/nnp-data-1/input.nn
cp RuNNer/HDNNP_2/input.nn mode2/nnp-data-2/input.nn
cp RuNNer/HDNNP_1/scaling.data mode2/nnp-data-1
cp RuNNer/HDNNP_1/weights.*.data mode2/nnp-data-1
cp RuNNer/HDNNP_2/scaling.data mode2/nnp-data-2
cp RuNNer/HDNNP_2/weights.*.data mode2/nnp-data-2
cd mode2
ln -s ../input.data-new input.data
#sed s/'.*runner_mode.*'/'runner_mode 1'/g RuNNer/HDNNP_1/input.nn > mode2/HDNNP_1/input.nn
#sed s/'.*runner_mode.*'/'runner_mode 1'/g RuNNer/HDNNP_2/input.nn > mode2/HDNNP_2/input.nn
#sed -i s/'.*test_fraction.*'/'test_fraction 0.0'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*epochs.*'/'epochs 0'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*use_old_scaling'/'use_old_scaling'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*use_old_weights_short'/'use_old_weights_short'/g mode2/HDNNP_*/input.nn
#sed -i s/"#CAUTION: don't forget use_short_forces below (if you want to generate the training files for the forces)"/"#CAUTION: don't forget use-short-forces below (if you want to generate the training files for the forces)"/g mode2/HDNNP_*/input.nn
#sed -i s/'.*use_short_forces'/'use_short_forces'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*remove_atom_energies'/'#remove_atom_energies'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*atom_energy'/'#atom_energy'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*energy_threshold'/'#energy_threshold'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*write_trainpoints'/'write_trainpoints'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*write_trainforces'/'write_trainforces'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*precondition_weights'/'#precondition_weights'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*use_systematic_weights_short'/'#use_systematic_weights_short'/g mode2/HDNNP_*/input.nn
#sed -i s/'.*nguyen_widrow_weights_short'/'#nguyen_widrow_weights_short'/g mode2/HDNNP_*/input.nn

#cd mode2/nnp-data-1
#  ln -s ../../input.data-new input.data
##  ${RuNNer_exe} > mode_1.out &
#cd ../../mode2/nnp-data-2
#  ln -s ../../input.data-new input.data
# ${RuNNer_exe} > mode_1.out &
#cd ../
wait
module load impi mkl intel gsl/2.7 n2p2/2.1.4-lammps-29Sept2021
mpirun -np 48 nnp-comp2 compare
nnp-comp2 select 1 $energy_threshold $force_threshold Cu O
#nnp-comp2 select 1 0.002 0.02 Cu O
#n_atoms=$(head -n 1 mode2/HDNNP_1/function.data)
#pbc=$(head -n 1 mode2/HDNNP_1/trainstruct.data | awk '{print $2}')
#if [ ${pbc} == T ]
#then
#  lattice=3
#else
#  lattice=0
#fi
#
#head -n $[${n_atoms}+2] mode2/HDNNP_1/function.data > mode2/HDNNP_1/testing.data
#head -n $[${n_atoms}+1] mode2/HDNNP_1/trainforces.data > mode2/HDNNP_1/testforces.data
#head -n $[${n_atoms}+1+${lattice}] mode2/HDNNP_1/trainstruct.data > mode2/HDNNP_1/teststruct.data
#head -n $[${n_atoms}+2] mode2/HDNNP_2/function.data > mode2/HDNNP_2/testing.data
#head -n $[${n_atoms}+1] mode2/HDNNP_2/trainforces.data > mode2/HDNNP_2/testforces.data
#head -n $[${n_atoms}+1+${lattice}] mode2/HDNNP_2/trainstruct.data > mode2/HDNNP_2/teststruct.data
#sed -i s/'.*runner_mode.*'/'runner_mode 2'/g mode2/HDNNP_*/input.nn
#
#cd mode2/HDNNP_1
#  ${RuNNer_exe} > mode_2.out &
#cd ../../mode2/HDNNP_2
#  ${RuNNer_exe} > mode_2.out &
#cd ../..
#wait
#
#rm mode2/HDNNP_*/000000.short.*.out mode2/HDNNP_*/debug.out mode2/HDNNP_*/function.data mode2/HDNNP_*/testforces.000000.out mode2/HDNNP_*/testforces.data mode2/HDNNP_*/testing.data mode2/HDNNP_*/testpoints.000000.out mode2/HDNNP_*/teststruct.data mode2/HDNNP_*/trainforces.data mode2/HDNNP_*/trainstruct.data
#
####################################################################################################

####################################################################################################
