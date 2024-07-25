#!/usr/bin/env bash

# Set the base run path
runpath=$(pwd)

# Loop through iterations from 1 to 5
for loop in {1..5}; do
    echo "################# Start Loop $loop #################"
    
    # Navigate to the loop directory
    cd "$runpath/AL-$loop"
    
    # Update seed value in the Python script
    sed -i "s/seed =.*/seed = $loop/g" RuNNerActiveLearn.py
    
    # Execute initial setup
    ./initial
    
    # Collect summary of selected structures
    echo $loop >> ../selectsum
    grep "Total number of structures   :" mode2.out >> ../selectsum
    grep "selected structures" mode2.out >> ../selectsum
    grep "threshold exceeded" mode2.out >> ../selectsum
    
    # Convert input data to POSCAR format
    runner2poscar.sh input.data-add
    
    # Check for POSCAR files
    if [ ! -f POSCAR.1 ]; then
        echo "No POSCAR from mode1"
        break
    fi
    
    if [ -f POSCAR.2000 ]; then
        echo "Too many POSCAR from mode1 > 2000"
        break
    fi
    
    # Organize POSCAR files into separate folders
    rm -rf poscar
    mkdir poscar
    for x in POSCAR.*; do
        cp -rf ../templ/ poscar/$x
        mv $x poscar/$x/POSCAR
    done
    
    cd poscar
    numberpos=$(ls | wc -l)
    echo "################# Submitting $numberpos Static #################"
    

    # Submit jobs for POSCAR folders
    for posfolder in *; do
        while [ -f "$posfolder/unsub" ]; do
            checklist=$(squeue | awk '{ print $1 }' | grep "\[" | wc -l)
            while [[ $checklist -ne 0 ]]; do
                echo "There is a list of jobs, waiting"
                sleep 300
                checklist=$(squeue | grep "\[" | wc -l)
            done
            subatparentfolder.sh
            echo "Check and submit jobs every 5 mins"
            sleep 300
        done
    done
    
    # Check and resubmit unfinished jobs
    for posfolder in *; do
        while [ -f "$posfolder/unfinished" ]; do
            cd $posfolder
            checkandresub.sh
            cd ..
            sleep 300
        done
    done
    
    # Combine output data
    cat */OUTCAR.1.data-deleunconverged > input.add$loop
    echo "################# Preparing Next Loop Folder #################"
    
    # Prepare next loop folder
    cp -rf "$runpath/AL-$loop" "$runpath/AL-$((loop+1))"
    rm "$runpath/AL-$((loop+1))/train.sum"
    cp input.add$loop "$runpath/AL-$loop/input.add$loop-DFT"
    cat input.add$loop >> "$runpath/AL-$((loop+1))/RuNNer/input.data"
    
    echo "############# Starting Training in Folder $((loop+1)) ###########"
    
    # Train new neural networks
    cd "$runpath/AL-$((loop+1))/RuNNer/HDNNP_1/train"
    rm *out
    rm input.data
    ln -s ../../input.data
    sbatch job
    
    cd "$runpath/AL-$((loop+1))/RuNNer/HDNNP_2/train"
    rm *out
    rm input.data
    ln -s ../../input.data
    sbatch job
    
    # Wait for training to finish
    checktrain=$(squeue | grep n2p2tra | wc -l)
    while [[ $checktrain -ne 0 ]]; do
        sleep 300
        checktrain=$(squeue | grep n2p2tra | wc -l)
    done
    
    # Copy trained model files to the next loop folder
    cp input.nn scaling.data weights.*data ../
    cd "$runpath/AL-$((loop+1))/RuNNer/HDNNP_1/train"
    cp input.nn scaling.data weights.*data ../
    
    echo "################# Finish Loop $loop #################"
    
    # Cleanup and prepare for the next loop
    cd "$runpath/AL-$((loop+1))"
    rm -rf mode2.out mode1 mode2
    cat train.sum
    echo $((loop+1)) >> "$runpath/train.sum"
    cat train.sum >> "$runpath/train.sum"
    echo "" >> "$runpath/train.sum"
    
    cd $runpath
    sed -i "s/loop=.*/loop=$loop/g" subzip
    sbatch subzip
done

echo "################# End All Loops #################"
