#!/bin/bash
#

thisrun='090-GM_NNB_IR3'
BUILD=Build_Parallel_Optim #Optimmax #Debugmax #Parallel_Debug # Parallel_Optim    
NTHREADS=24  # number of OpenMP threads 
# OMP_STACKSIZE=16M #stack size for each OMP thread: default 4M, recommended 16M, mytest 512M

MYDIR="$( cd "$( dirname "$0" )" && pwd )"
projectname="${MYDIR##*/}"
MY_IP="$(wget -q -O - checkip.dyndns.org|sed -e 's/.*Current IP Address: //' -e 's/<.*$//')"
USERNAME="$(whoami)"

path_to_results=~/work/Research/Eqprem_idiorisk/test_results/Euler
mkdir --parents $path_to_results/$thisrun
cd $MYDIR

cd model_input/last_results
# The next line is crazy: it finds the name of all selected calibrations (not outcommented, and without the .txt) and pipes them to tar. Note it handles also several calibs.
grep -Po '^[^!].*?(?<=/)\K[^/.]*(?=\.)' ../select_calibration_here.txt | xargs tar -cjf last_results.tar.bz2
cd ../..

mv eu-* eu-$thisrun
sed -i "s/epss4\w*/$projectname/g" eu-$thisrun
tar -zcf code.tar.gz *.f90 eu* src*
mv code.tar.gz $path_to_results/$thisrun
rm $BUILD/*.o $BUILD/*.mod $BUILD/src_utilities/*.o $BUILD/src_utilities/*.mod $BUILD/src_classes/*.o $BUILD/src_classes/*.mod

tar -cf - *.f90 eu-$thisrun $BUILD src_classes src_matlab src_utilities model_input/*.* model_input/calib* model_input/data model_input/last_results/last_results.tar.bz2 | ssh -C danielh@euler.ethz.ch "export BUILD=$BUILD OMP_NUM_THREADS=$NTHREADS MY_IP=$MY_IP USERNAME=$USERNAME path_to_results=$path_to_results && mkdir -p $thisrun/model_output && cd $thisrun && tar -xf - ; cd model_input/last_results; tar -xjf last_results.tar.bz2; rm last_results.tar.bz2; cd ../..; bsub -J $thisrun < eu-$thisrun"
rm model_input/last_results/last_results.tar.bz2