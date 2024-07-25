for x in */;do cd $x;grep thermo log.lammps > data.txt;sed -i '/#/d' data.txt;sed -i '/v_thermo/d' data.txt ;cd ..;done
#for x in */;do cd $x;grep thermo log.lammps > data.txt;sed -i '341,347d' data.txt;sed -i '675,681d' data.txt ;cd ..;done
